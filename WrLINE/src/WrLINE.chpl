/*
 * Copyright (C) 2026 Alexander Dukhan
 * Chapel transcription of the WrLINE Python package developed by
 * Victor Manuel Velasco Berrelleza and Agnes Noy.
 *
 * WrLINE processes C1' atomic coordinates from molecular dynamics
 * trajectories to compute helical axis, twist, register angles,
 * and writhe of closed/circular DNA structures.
 */

module WrLINE {
  use Parms;
  use Functions;
  use Math;
  use IO;
  use FileSystem;

  config var name, top, traj: string = "";
  config var nbp, nstep: int = 0;
  config var doStrip, isCircular: bool = false;

  // -----------------------------------------------------------------------
  // TUPLE-BASED VECTOR MATH (stack-allocated, no heap overhead)
  // -----------------------------------------------------------------------

  type vec3 = 3*real;      // (real, real, real) — stack-allocated
  type mat3 = 3*(3*real);  // 3 row-vectors — stack-allocated

  proc mkVec(x: real, y: real, z: real): vec3 {
    return (x, y, z);
  }

  proc vecAdd(a: vec3, b: vec3): vec3 {
    return (a(0) + b(0), a(1) + b(1), a(2) + b(2));
  }

  proc vecSub(a: vec3, b: vec3): vec3 {
    return (a(0) - b(0), a(1) - b(1), a(2) - b(2));
  }

  proc vecScale(a: vec3, s: real): vec3 {
    return (a(0) * s, a(1) * s, a(2) * s);
  }

  proc vecDot(a: vec3, b: vec3): real {
    return a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
  }

  proc vecNorm(a: vec3): real {
    return sqrt(a(0)**2 + a(1)**2 + a(2)**2);
  }

  proc vecNormalize(a: vec3): vec3 {
    var n = vecNorm(a);
    return (a(0) / n, a(1) / n, a(2) / n);
  }

  proc vecCross(a: vec3, b: vec3): vec3 {
    return (a(1) * b(2) - a(2) * b(1),
            a(2) * b(0) - a(0) * b(2),
            a(0) * b(1) - a(1) * b(0));
  }

  // Matrix-vector multiply: M * v
  proc matVec(m: mat3, v: vec3): vec3 {
    return (vecDot(m(0), v),
            vecDot(m(1), v),
            vecDot(m(2), v));
  }

  // Matrix-matrix multiply: A * B
  proc matMul(a: mat3, b: mat3): mat3 {
    // Columns of B
    var b0 = mkVec(b(0)(0), b(1)(0), b(2)(0)),
        b1 = mkVec(b(0)(1), b(1)(1), b(2)(1)),
        b2 = mkVec(b(0)(2), b(1)(2), b(2)(2));
    return ((vecDot(a(0), b0), vecDot(a(0), b1), vecDot(a(0), b2)),
            (vecDot(a(1), b0), vecDot(a(1), b1), vecDot(a(1), b2)),
            (vecDot(a(2), b0), vecDot(a(2), b1), vecDot(a(2), b2)));
  }

  // -----------------------------------------------------------------------
  // MATHEMATICAL OPERATIONS (from caxislib.py) — using tuples
  // -----------------------------------------------------------------------

  /*
   * Find arctan of y/x and return values between -pi and pi.
   */
  proc arctan360(x: real, y: real): real {
    if x >= 0.0 && y >= 0.0 {
      return atan(y / x);        // Q1
    } else if x < 0.0 && y >= 0.0 {
      return atan(y / x) + pi;   // Q2
    }
    else if x < 0.0 && y < 0.0 {
      return atan(y / x) - pi;   // Q3
    } else {
      return atan(y / x);        // Q4
    }
  }

  /*
   * Returns rotation matrix that rotates a vector to the x-axis.
   */
  proc rotateToX(v: vec3): mat3 {
    var n = vecNormalize(v);
    var c = -atan2(n(1), n(0));
    return ((cos(c), -sin(c), 0.0),
            (sin(c),  cos(c), 0.0),
            (0.0,     0.0,    1.0));
  }

  /*
   * Returns rotation matrix that rotates a vector to the z-axis.
   */
  proc rotateToZ(v: vec3): mat3 {
    var n = vecNormalize(v);
    var a = atan2(n(1), n(2));
    var b = atan(-n(0) / sqrt(n(1)**2 + n(2)**2));

    var Rx = ((1.0, 0.0, 0.0),
              (0.0, cos(a), -sin(a)),
              (0.0, sin(a),  cos(a)));
    var Ry = ((cos(b), 0.0, sin(b)),
              (0.0,    1.0, 0.0),
              (-sin(b), 0.0, cos(b)));

    return matMul(Ry, Rx);
  }

  /*
   * Returns twist angle between two vectors (a1..a2) and (b1..b2)
   * relative to a local z-axis.
   */
  proc twistAngle(a1: vec3, a2: vec3, b1: vec3, b2: vec3, z: vec3): real {
    var vA = vecNormalize(vecSub(a2, a1));
    var vB = vecNormalize(vecSub(b2, b1));
    var vZ = vecNormalize(z);

    var Rz = rotateToZ(vZ);
    var unitA = matVec(Rz, vA);
    var unitB = matVec(Rz, vB);

    var Rx = rotateToX(unitA);
    var rotatedB = matVec(Rx, unitB);

    return atan2(rotatedB(1), rotatedB(0)) * 180.0 / pi;
  }

  // -----------------------------------------------------------------------
  // FILE I/O
  // -----------------------------------------------------------------------

  /*
   * Read 3-column coordinate file.
   */
  proc readThreeCol(path: string, inNbp: int, inNstep: int) {
    var X: [1..inNstep, 1..inNbp, 1..3] real;
    try! {
      var f = open(path, ioMode.r);
      var reader = f.reader();
      for t in 1..inNstep {
        for j in 1..inNbp {
          for c in 1..3 {
            reader.read(X[t, j, c]);
          }
        }
      }
      f.close();
    }
    return X;
  }

  /*
   * Read AMBER mdcrd file and return strandA, strandB and midpoints.
   */
  proc readMdcrd(inNbp: int, inNstep: int, inTraj: string) {
    try! {
      var f = open(inTraj, ioMode.r);
      var reader = f.reader(locking=false);

      // Skip header line
      var headerLine: string;
      if !(try! reader.readLine(headerLine)) then halt("Empty mdcrd file");

      // Total coordinates: nstep * nbp * 2 atoms * 3 components
      var totalCoords = inNstep * inNbp * 2 * 3;
      var coords: [0..#totalCoords] real;
      var idx = 0;

      // Read all coordinate values (8 characters wide, fixed format)
      for line in reader.lines(stripNewline=true) {
        var l = line.numBytes;
        var n = l / 8;
        for i in 0..#n {
          if idx < totalCoords {
            var token = line[i * 8..#8];
            var isBlank = true;
            for b in token.bytes() do
              if b != 32 {
                isBlank = false;
                break;
              }
            if !isBlank {
              coords[idx] = token: real;
              idx += 1;
            }
          }
        }
      }
      f.close();

      // Build strand arrays directly from flat coords
      var rA, rB: [1..inNstep, 1..inNbp, 1..3] real;

      // coords layout: for each atom (nstep * nbp*2 atoms), 3 values
      var atomsPerStep = inNbp * 2;
      for t in 1..inNstep {
        for j in 1..inNbp {
          // strand_a = data_array[0:num_bp]
          var aAtomIdx = (j - 1);
          var aIdx = ((t - 1) * atomsPerStep + aAtomIdx) * 3;
          rA[t, j, 1] = coords[aIdx];
          rA[t, j, 2] = coords[aIdx + 1];
          rA[t, j, 3] = coords[aIdx + 2];

          // strand_b = data_array[(2*num_bp - 1):(num_bp - 1):-1]
          var bAtomIdx = (2 * inNbp - j);
          var bIdx = ((t - 1) * atomsPerStep + bAtomIdx) * 3;
          rB[t, j, 1] = coords[bIdx];
          rB[t, j, 2] = coords[bIdx + 1];
          rB[t, j, 3] = coords[bIdx + 2];
        }
      }

      // r = midpoints of neighboring bp steps
      var r: [1..inNstep, 1..inNbp, 1..3] real;
      for j in 1..inNbp {
        var jNext = (j % inNbp) + 1;
        for t in 1..inNstep {
          for c in 1..3 {
            if !isCircular && j == inNbp {
              r[t, j, c] = 0.5 * (rA[t, j, c] + rB[t, j, c]);
            } else {
              r[t, j, c] = 0.25 * (rA[t, j, c] + rA[t, jNext, c]
                                   + rB[t, j, c] + rB[t, jNext, c]);
            }
          }
        }
      }

      return (rA, rB, r);
    }
  }

  // -----------------------------------------------------------------------
  // MAIN DATA PROCESSING
  // -----------------------------------------------------------------------

  /*
   * 1st-order helical axis
   */
  proc haxis(inNbp: int, inNstep: int,
             ref r: [1..inNstep, 1..inNbp, 1..3] real) {
    var rC: [1..inNstep, 1..inNbp, 1..3] real;
    forall t in 1..inNstep {
      for j in 1..inNbp {
        var sx = r[t, j, 1], sy = r[t, j, 2], sz = r[t, j, 3];
        var k = 0;
        while k < 5 {
          k += 1;
          if !isCircular {
            // Python: midpoints[:, t, j-k] works (wraps), 
            //         midpoints[:, t, j+k] can IndexError
            if j + k > inNbp {
              k -= 1;
              break;
            }
            var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
            var jp = j + k;
            sx += r[t, jm, 1] + r[t, jp, 1];
            sy += r[t, jm, 2] + r[t, jp, 2];
            sz += r[t, jm, 3] + r[t, jp, 3];
          } else {
            var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
            var jp = ((j - 1 + k) % inNbp) + 1;
            sx += r[t, jm, 1] + r[t, jp, 1];
            sy += r[t, jm, 2] + r[t, jp, 2];
            sz += r[t, jm, 3] + r[t, jp, 3];
          }
        }
        var denom = (2.0 * k + 1.0);
        rC[t, j, 1] = sx / denom;
        rC[t, j, 2] = sy / denom;
        rC[t, j, 3] = sz / denom;
      }
    }
    return rC;
  }

  /*
   * Calculate twist
   */
  proc calcTwist(inName: string, inNbp: int, inNstep: int,
                 ref rA: [1..inNstep, 1..inNbp, 1..3] real,
                 ref rB: [1..inNstep, 1..inNbp, 1..3] real,
                 ref rC: [1..inNstep, 1..inNbp, 1..3] real) {
    var tw: [1..inNstep, 1..inNbp] real;
    forall t in 1..inNstep {
      for j in 1..inNbp {
        var Z: vec3;
        if isCircular {
          var jNext = (j % inNbp) + 1;
          var jPrev = ((j - 2 + inNbp) % inNbp) + 1;
          Z = mkVec(rC[t, jNext, 1] - rC[t, jPrev, 1],
                    rC[t, jNext, 2] - rC[t, jPrev, 2],
                    rC[t, jNext, 3] - rC[t, jPrev, 3]);
          tw[t, j] = twistAngle(mkVec(rA[t, j, 1], rA[t, j, 2], rA[t, j, 3]),
                       mkVec(rB[t, j, 1], rB[t, j, 2], rB[t, j, 3]),
                       mkVec(rA[t, jNext, 1], rA[t, jNext, 2], rA[t, jNext, 3]),
                       mkVec(rB[t, jNext, 1], rB[t, jNext, 2], rB[t, jNext, 3]),
                       Z);
        } else {
          if j < inNbp {
            if j > 1 {
              Z = mkVec(rC[t, j+1, 1] - rC[t, j-1, 1],
                        rC[t, j+1, 2] - rC[t, j-1, 2],
                        rC[t, j+1, 3] - rC[t, j-1, 3]);
            } else {
              Z = mkVec(2.0 * (rC[t, j+1, 1] - rC[t, j, 1]),
                        2.0 * (rC[t, j+1, 2] - rC[t, j, 2]),
                        2.0 * (rC[t, j+1, 3] - rC[t, j, 3]));
            }
            tw[t, j] = twistAngle(mkVec(rA[t, j, 1], rA[t, j, 2], rA[t, j, 3]),
                         mkVec(rB[t, j, 1], rB[t, j, 2], rB[t, j, 3]),
                         mkVec(rA[t, j+1, 1], rA[t, j+1, 2], rA[t, j+1, 3]),
                         mkVec(rB[t, j+1, 1], rB[t, j+1, 2], rB[t, j+1, 3]),
                         Z);
          } else {
            tw[t, j] = 0.0;
          }
        }
      }
    }

    // Write tw.ser (array block I/O)
    try! {
      var f = open(inName + "/tw.ser", ioMode.cw);
      var writer = f.writer(locking=false);
      for t in 1..inNstep do
        writer.writeln(tw[t, ..]);
      writer.close();
    }

    return tw;
  }

  /*
   * Calculate central helical axis (CAXIS).
   */
  proc caxis(inNbp: int, inNstep: int,
             ref r: [1..inNstep, 1..inNbp, 1..3] real,
             ref tw: [1..inNstep, 1..inNbp] real) {
    var r1: [1..inNstep, 1..inNbp, 1..3] real;
    forall t in 1..inNstep {
      for j in 1..inNbp {
        var Tw = tw[t, j];
        var sx = r[t, j, 1], sy = r[t, j, 2], sz = r[t, j, 3];
        var k = 0;
        var prev = Tw;
        var triggeredBreak = false;

        while Tw < 360.0 {
          k += 1;
          prev = Tw;
          
          if !isCircular {
            if j + k > inNbp {
              triggeredBreak = true;
              break;
            }
          }
          if isCircular && k > inNbp then break;
          
          var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
          var jp = if !isCircular then j + k else ((j - 1 + k) % inNbp) + 1;
          
          Tw += tw[t, jm] + tw[t, jp];
          sx += r[t, jm, 1] + r[t, jp, 1];
          sy += r[t, jm, 2] + r[t, jp, 2];
          sz += r[t, jm, 3] + r[t, jp, 3];
        }

        var w: real = 0.0;
        if !isCircular && triggeredBreak {
           w = 0.0;
        } else {
          if Tw != prev then
            w = (360.0 - prev) / (Tw - prev);
          else
            w = 0.0;
            
          var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
          var jp = if !isCircular then j + k else ((j - 1 + k) % inNbp) + 1;
          
          sx -= (1.0 - w) * (r[t, jm, 1] + r[t, jp, 1]);
          sy -= (1.0 - w) * (r[t, jm, 2] + r[t, jp, 2]);
          sz -= (1.0 - w) * (r[t, jm, 3] + r[t, jp, 3]);
        }

        var denom = 2.0 * (k: real + w) - 1.0;
        if denom == 0 then denom = 1.0;
        r1[t, j, 1] = sx / denom;
        r1[t, j, 2] = sy / denom;
        r1[t, j, 3] = sz / denom;
      }
    }
    return r1;
  }

  /*
   * Calculate sine of the register angles.
   */
  proc sinreg(inName: string, inNbp: int, inNstep: int,
              ref r: [1..inNstep, 1..inNbp, 1..3] real,
              ref r1: [1..inNstep, 1..inNbp, 1..3] real) {
    var sinRegArr: [1..inNstep, 0..inNbp] real;
    forall t in 1..inNstep {
      for j in 1..inNbp {
        var jPrev = ((j - 2 + inNbp) % inNbp) + 1;
        var jNext = (j % inNbp) + 1;

        var v0 = mkVec(r1[t, j, 1] - r1[t, jPrev, 1],
                       r1[t, j, 2] - r1[t, jPrev, 2],
                       r1[t, j, 3] - r1[t, jPrev, 3]);
        var v1 = mkVec(r1[t, jNext, 1] - r1[t, j, 1],
                       r1[t, jNext, 2] - r1[t, j, 2],
                       r1[t, jNext, 3] - r1[t, j, 3]);

        var planeVec = vecCross(v1, v0);
        var minorGroove = mkVec(r[t, jPrev, 1] - r1[t, jPrev, 1],
                                r[t, jPrev, 2] - r1[t, jPrev, 2],
                                r[t, jPrev, 3] - r1[t, jPrev, 3]);

        var normM = vecNorm(minorGroove);
        var normP = vecNorm(planeVec);

        if normM > 0 && normP > 0 {
          var val = vecNorm(vecCross(minorGroove, planeVec)) / (normM * normP);
          var dotVal = vecDot(vecSub(v1, v0), minorGroove);
          if dotVal < 0 then val = -val;
          sinRegArr[t, j] = val;
        } else {
          sinRegArr[t, j] = 0.0;
        }
      }
    }

    // Time scale column
    for t in 1..inNstep do sinRegArr[t, 0] = 0.01 * t;

    try! {
      var f = open(inName + "/sinreg.ser", ioMode.cw);
      var writer = f.writer(locking=false);
      for t in 1..inNstep do
        writer.writeln(sinRegArr[t, ..]);
      writer.close();
    }
  }

  /*
   * Output files
   */
  proc makeFiles(inName: string, inNbp: int, inNstep: int,
                 ref r: [1..inNstep, 1..inNbp, 1..3] real,
                 ref r1: [1..inNstep, 1..inNbp, 1..3] real) {
    try! {
      {
        var f = open(inName + "/C.xyz", ioMode.cw);
        var w = f.writer(locking=false);
        for i in 1..inNstep {
          w.writeln(inNbp);
          w.writeln();
          for j in 1..inNbp do
            w.writef("H %8.3dr %8.3dr %8.3dr \n",
                    r[i, j, 1], r[i, j, 2], r[i, j, 3]);
        }
      }
      {
        var f = open(inName + "/C1.xyz", ioMode.cw);
        var w = f.writer(locking=false);
        for i in 1..inNstep {
          w.writeln(inNbp);
          w.writeln();
          for j in 1..inNbp do
            w.writef("H %8.3dr %8.3dr %8.3dr \n",
                    r1[i, j, 1], r1[i, j, 2], r1[i, j, 3]);
        }
      }
      {
        var f = open(inName + "/C.3col", ioMode.cw);
        var w = f.writer(locking=false);
        for i in 1..inNstep do for j in 1..inNbp do
          w.writef("%8.3dr %8.3dr %8.3dr \n",
                  r[i, j, 1], r[i, j, 2], r[i, j, 3]);
      }
      {
        var f = open(inName + "/C1.3col", ioMode.cw);
        var w = f.writer(locking=false);
        for i in 1..inNstep do for j in 1..inNbp do
          w.writef("%8.3dr %8.3dr %8.3dr \n",
                 r1[i, j, 1], r1[i, j, 2], r1[i, j, 3]);
      }
    }
  }

  /*
   * Writhe
   */
  proc wr(ref X: [] real, t: int, l: int): real {
    var Wr: real = 0.0;
    var loopLimit = if isCircular then l else l - 1;
    // Removed nested forall loop and local array 'y' allocations
    for j in 1..loopLimit {
      var next_j = if isCircular && j == l then 1 else j + 1;
      var tj = mkVec(X[t, next_j, 1] - X[t, j, 1],
                     X[t, next_j, 2] - X[t, j, 2],
                     X[t, next_j, 3] - X[t, j, 3]);
      for k in 1..j - 1 {
        var next_k = k + 1;
        var tk = mkVec(X[t, next_k, 1] - X[t, k, 1],
                       X[t, next_k, 2] - X[t, k, 2],
                       X[t, next_k, 3] - X[t, k, 3]);
        var rjk = mkVec(X[t, j, 1] - X[t, k, 1],
                       X[t, j, 2] - X[t, k, 2],
                       X[t, j, 3] - X[t, k, 3]);
        var crossJK = vecCross(tj, tk);
        var dotVal = vecDot(rjk, crossJK);
        var dist = vecNorm(rjk);
        if dist > 1e-9 {
          Wr += dotVal / (dist * dist * dist * 2.0 * pi);
        }
      }
    }
    return Wr;
  }

  proc calcWrithe(inName: string, inNbp: int,
                  inNstep: int, ref X: [1..inNstep, 1..inNbp, 1..3] real) {
    var writheArr: [1..inNstep, 1..2] real;
    forall t in 1..inNstep {
      writheArr[t, 1] = t: real;
      writheArr[t, 2] = wr(X, t, inNbp);
    }
    try! {
      var f = open(inName + "/writhe.ser", ioMode.cw);
      var writer = f.writer(locking=false);
      for t in 1..inNstep do
        writer.writeln(writheArr[t, 1]: int, " ", writheArr[t, 2]);
    }
  }

  proc main() {
    try! {
      if name == "" || nbp == 0 || nstep == 0 {
        writeln(
    "Usage: WrLINE --name=<name> --traj=<traj> --nbp=<nbp> --nstep=<nstep> ..."
        );
        return;
      }
      if doStrip {
        use Subprocess;
        var sub = spawn(["bash", "stripC.sh", name, top, traj]);
        sub.wait();
      }
      if !isDir(name) then mkdir(name, parents=true);
      writeln("Reading trajectory...");
      var (rA, rB, r) = readMdcrd(nbp, nstep, traj);
      writeln("Calculating helical axis...");
      var rC = haxis(nbp, nstep, r);
      writeln("Calculating twist...");
      var tw = calcTwist(name, nbp, nstep, rA, rB, rC);
      writeln("Calculating central helical axis...");
      var r1 = caxis(nbp, nstep, r, tw);
      writeln("Calculating register angles...");
      sinreg(name, nbp, nstep, r, r1);
      writeln("Generating coordinate files...");
      makeFiles(name, nbp, nstep, r, r1);
      writeln("Calculating writhe...");
      calcWrithe(name, nbp, nstep, r1);
      writeln("job ", name, " done!!! XD!! :) ");
    }
  }
}
