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

  config var name: string = "";
  config var top: string = "";
  config var traj: string = "";
  config var nbp: int = 0;
  config var nstep: int = 0;
  config var doStrip: bool = false;
  config var isCircular: bool = false;

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
    var b0 = mkVec(b(0)(0), b(1)(0), b(2)(0));
    var b1 = mkVec(b(0)(1), b(1)(1), b(2)(1));
    var b2 = mkVec(b(0)(2), b(1)(2), b(2)(2));
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
    } else 
    if x < 0.0 && y < 0.0 {
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
    var inp = inTraj;
    try! {
      var f = open(inp, ioMode.r);
      var reader = f.reader();

      // Skip header line
      var headerLine: string;
      reader.readLine(headerLine);

      // Total coordinates: nstep * nbp * 2 atoms * 3 components
      var totalCoords = inNstep * inNbp * 2 * 3;
      var coords: [0..#totalCoords] real;
      var idx = 0;

      // Read all coordinate values (8 characters wide, fixed format)
      var line: string;
      while reader.readLine(line) {
        var l = line.size;
        if l > 0 && line[l - 1] == "\n" then l -= 1;
        var n = l / 8;
        for i in 0..#n {
          if idx < totalCoords {
            var token = line[i * 8..#8];
            coords[idx] = token: real;
            idx += 1;
          }
        }
      }
      f.close();

      // Build strand arrays directly from flat coords
      var rA: [1..3, 1..inNstep, 1..inNbp] real;
      var rB: [1..3, 1..inNstep, 1..inNbp] real;

      // coords layout: for each atom (nstep * nbp*2 atoms), 3 values
      var atomsPerStep = inNbp * 2;
      for t in 1..inNstep {
        for j in 1..inNbp {
          var aIdx = ((t - 1) * atomsPerStep + (j - 1)) * 3;
          rA[1, t, j] = coords[aIdx];
          rA[2, t, j] = coords[aIdx + 1];
          rA[3, t, j] = coords[aIdx + 2];

          var bAtom = 2 * inNbp - j; 
          var bIdx = ((t - 1) * atomsPerStep + bAtom) * 3;
          rB[1, t, j] = coords[bIdx];
          rB[2, t, j] = coords[bIdx + 1];
          rB[3, t, j] = coords[bIdx + 2];
        }
      }

      // r = midpoints of neighboring bp steps
      var r: [1..3, 1..inNstep, 1..inNbp] real;
      for j in 1..inNbp {
        var jNext = (j % inNbp) + 1;
        for t in 1..inNstep {
          for c in 1..3 {
            if !isCircular && j == inNbp {
              r[c, t, j] = 0.5 * (rA[c, t, j] + rB[c, t, j]);
            } else {
              r[c, t, j] = 0.25 * (rA[c, t, j] + rA[c, t, jNext]
                                   + rB[c, t, j] + rB[c, t, jNext]);
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
             ref r: [1..3, 1..inNstep, 1..inNbp] real) {
    var rC: [1..3, 1..inNstep, 1..inNbp] real;
    forall j in 1..inNbp {
      for t in 1..inNstep {
        var sx = r[1, t, j], sy = r[2, t, j], sz = r[3, t, j];
        var k_actual = 0;
        for k in 1..5 {
          if isCircular {
            var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
            var jp = ((j - 1 + k) % inNbp) + 1;
            sx += r[1, t, jm] + r[1, t, jp];
            sy += r[2, t, jm] + r[2, t, jp];
            sz += r[3, t, jm] + r[3, t, jp];
            k_actual = k;
          } else {
            // Replicate Python reference behavior: negative index wraps
            // around (Python's negative indexing), but positive out-of-range
            // raises IndexError and breaks.
            if j + k > inNbp {
              break;
            }
            var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
            var jp = j + k;
            sx += r[1, t, jm] + r[1, t, jp];
            sy += r[2, t, jm] + r[2, t, jp];
            sz += r[3, t, jm] + r[3, t, jp];
            k_actual = k;
          }
        }
        var denom = (2.0 * k_actual + 1.0);
        rC[1, t, j] = sx / denom;
        rC[2, t, j] = sy / denom;
        rC[3, t, j] = sz / denom;
      }
    }
    return rC;
  }

  /*
   * Calculate twist
   */
  proc calcTwist(inName: string, inNbp: int, inNstep: int,
                 ref rA: [1..3, 1..inNstep, 1..inNbp] real,
                 ref rB: [1..3, 1..inNstep, 1..inNbp] real,
                 ref rC: [1..3, 1..inNstep, 1..inNbp] real) {
    var tw: [1..inNstep, 1..inNbp] real;
    forall j in 1..inNbp {
      for t in 1..inNstep {
        var Z: vec3;
        if isCircular {
          var jNext = (j % inNbp) + 1;
          var jPrev = ((j - 2 + inNbp) % inNbp) + 1;
          Z = mkVec(rC[1, t, jNext] - rC[1, t, jPrev],
                    rC[2, t, jNext] - rC[2, t, jPrev],
                    rC[3, t, jNext] - rC[3, t, jPrev]);
          tw[t, j] = twistAngle(mkVec(rA[1, t, j], rA[2, t, j], rA[3, t, j]),
                       mkVec(rB[1, t, j], rB[2, t, j], rB[3, t, j]),
                       mkVec(rA[1, t, jNext], rA[2, t, jNext], rA[3, t, jNext]),
                       mkVec(rB[1, t, jNext], rB[2, t, jNext], rB[3, t, jNext]),
                       Z);
        } else {
          if j < inNbp {
            if j > 1 {
              Z = mkVec(rC[1, t, j+1] - rC[1, t, j-1],
                        rC[2, t, j+1] - rC[2, t, j-1],
                        rC[3, t, j+1] - rC[3, t, j-1]);
            } else {
              Z = mkVec(2.0 * (rC[1, t, j+1] - rC[1, t, j]),
                        2.0 * (rC[2, t, j+1] - rC[2, t, j]),
                        2.0 * (rC[3, t, j+1] - rC[3, t, j]));
            }
            tw[t, j] = twistAngle(mkVec(rA[1, t, j], rA[2, t, j], rA[3, t, j]),
                         mkVec(rB[1, t, j], rB[2, t, j], rB[3, t, j]),
                         mkVec(rA[1, t, j+1], rA[2, t, j+1], rA[3, t, j+1]),
                         mkVec(rB[1, t, j+1], rB[2, t, j+1], rB[3, t, j+1]),
                         Z);
          } else {
            tw[t, j] = 0.0;
          }
        }
      }
    }

    // Write tw.ser (serial I/O)
    try! {
      var f = open(inName + "/tw.ser", ioMode.cw);
      var writer = f.writer();
      for t in 1..inNstep {
        for j in 1..inNbp {
          writer.write(tw[t, j], " ");
        }
        writer.writeln();
      }
      writer.close();
    }

    return tw;
  }

  /*
   * Calculate central helical axis (CAXIS).
   */
  proc caxis(inNbp: int, inNstep: int,
             ref r: [1..3, 1..inNstep, 1..inNbp] real,
             ref tw: [1..inNstep, 1..inNbp] real) {
    var r1: [1..3, 1..inNstep, 1..inNbp] real;
    forall j in 1..inNbp {
      for t in 1..inNstep {
        var Tw = tw[t, j];
        var sx = r[1, t, j], sy = r[2, t, j], sz = r[3, t, j];
        var k = 0;
        var prev = Tw;
        while Tw < 360.0 {
          k += 1;
          prev = Tw;
          if !isCircular && (j - k < 1 || j + k > inNbp) {
            break;
          }
          var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
          var jp = ((j - 1 + k) % inNbp) + 1;
          Tw += tw[t, jm] + tw[t, jp];
          sx += r[1, t, jm] + r[1, t, jp];
          sy += r[2, t, jm] + r[2, t, jp];
          sz += r[3, t, jm] + r[3, t, jp];
        }

        var w: real = 0.0;
        if isCircular || (j - k >= 1 && j + k <= inNbp) {
          if Tw != prev then
            w = (360.0 - prev) / (Tw - prev);
          else
            w = 0.0;
            
          var jmW = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
          var jpW = ((j - 1 + k) % inNbp) + 1;
          sx -= (1.0 - w) * (r[1, t, jmW] + r[1, t, jpW]);
          sy -= (1.0 - w) * (r[2, t, jmW] + r[2, t, jpW]);
          sz -= (1.0 - w) * (r[3, t, jmW] + r[3, t, jpW]);
        }

        var denom = 2.0 * (k: real + w) - 1.0;
        if denom == 0 then denom = 1.0;
        r1[1, t, j] = sx / denom;
        r1[2, t, j] = sy / denom;
        r1[3, t, j] = sz / denom;
      }
    }
    return r1;
  }

  /*
   * Calculate sine of the register angles.
   */
  proc sinreg(inName: string, inNbp: int, inNstep: int,
              ref r: [1..3, 1..inNstep, 1..inNbp] real,
              ref r1: [1..3, 1..inNstep, 1..inNbp] real) {
    var sinRegArr: [1..inNstep, 0..inNbp] real;
    forall j in 1..inNbp {
      for t in 1..inNstep {
        var jPrev = ((j - 2 + inNbp) % inNbp) + 1;
        var jNext = (j % inNbp) + 1;

        var v0 = mkVec(r1[1, t, j] - r1[1, t, jPrev],
                       r1[2, t, j] - r1[2, t, jPrev],
                       r1[3, t, j] - r1[3, t, jPrev]);
        var v1 = mkVec(r1[1, t, jNext] - r1[1, t, j],
                       r1[2, t, jNext] - r1[2, t, j],
                       r1[3, t, jNext] - r1[3, t, j]);

        var planeVec = vecCross(v1, v0);
        var minorGroove = mkVec(r[1, t, jPrev] - r1[1, t, jPrev],
                                r[2, t, jPrev] - r1[2, t, jPrev],
                                r[3, t, jPrev] - r1[3, t, jPrev]);

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
      var writer = f.writer();
      for t in 1..inNstep {
        for j in 0..inNbp {
          writer.write(sinRegArr[t, j], " ");
        }
        writer.writeln();
      }
      writer.close();
    }
  }

  /*
   * Output files
   */
  proc makeFiles(inName: string, inNbp: int, inNstep: int,
                 ref r: [1..3, 1..inNstep, 1..inNbp] real,
                 ref r1: [1..3, 1..inNstep, 1..inNbp] real) {
    coforall task in 1..4 {
      try! {
        select task { 
          when 1 {
            var f = open(inName + "/C.xyz", ioMode.cw);
            var w = f.writer();
            for i in 1..inNstep {
              w.writeln(inNbp); 
              w.writeln();
              for j in 1..inNbp do
                w.writef("H %8.3dr %8.3dr %8.3dr \n",
                        r[1,i,j], r[2,i,j], r[3,i,j]);
          }
        } when 2 {
          var f = open(inName + "/C1.xyz", ioMode.cw);
          var w = f.writer();
            for i in 1..inNstep {
              w.writeln(inNbp);
              w.writeln();
              for j in 1..inNbp do
                w.writef("H %8.3dr %8.3dr %8.3dr \n",
                        r1[1,i,j], r1[2,i,j], r1[3,i,j]);
          }
        } when 3 {
            var f = open(inName + "/C.3col", ioMode.cw);
            var w = f.writer();
            for i in 1..inNstep do for j in 1..inNbp do
              w.writef("%8.3dr %8.3dr %8.3dr \n",
                      r[1,i,j], r[2,i,j], r[3,i,j]);
        } when 4 {
            var f = open(inName + "/C1.3col", ioMode.cw);
            var w = f.writer();
            for i in 1..inNstep do for j in 1..inNbp do
              w.writef("%8.3dr %8.3dr %8.3dr \n",
                     r1[1,i,j], r1[2,i,j], r1[3,i,j]);
          }
        }
      }
    }
  }

  /*
   * Writhe
   */
  proc wr(ref X: [] real, t: int, l: int): real {
    var ySize = if isCircular then l + 1 else l;
    var y: [1..ySize, 1..3] real;
    for j in 1..l do for c in 1..3 do y[j, c] = X[t, j, c];
    if isCircular then for c in 1..3 do y[l+1, c] = X[t, 1, c];

    var Wr: real = 0.0;
    var loopLimit = if isCircular then l else l - 1;
    forall j in 1..loopLimit with (+ reduce Wr) {
      for k in 1..j - 1 {
        var tj = mkVec(y[j+1,1]-y[j,1], y[j+1,2]-y[j,2], y[j+1,3]-y[j,3]);
        var tk = mkVec(y[k+1,1]-y[k,1], y[k+1,2]-y[k,2], y[k+1,3]-y[k,3]);
        var rjk = mkVec(y[j,1]-y[k,1], y[j,2]-y[k,2], y[j,3]-y[k,3]);
        var dotVal = vecDot(rjk, vecCross(tj, tk));
        var distSq = vecDot(rjk, rjk);
        Wr += dotVal / (distSq * sqrt(distSq) * 2.0 * pi);
      }
    }
    return Wr;
  }

  proc calcWrithe(inName: string, inNbp: int, inNstep: int) {
    var X = readThreeCol(inName + "/C1.3col", inNbp, inNstep);
    var writheArr: [1..inNstep, 1..2] real;
    forall t in 1..inNstep {
      writheArr[t, 1] = t: real;
      writheArr[t, 2] = wr(X, t, inNbp);
    }
    try! {
      var f = open(inName + "/writhe.ser", ioMode.cw);
      var writer = f.writer();
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
      var (rA, rB, r) = readMdcrd(nbp, nstep, traj);
      var rC = haxis(nbp, nstep, r);
      var tw = calcTwist(name, nbp, nstep, rA, rB, rC);
      var r1 = caxis(nbp, nstep, r, tw);
      sinreg(name, nbp, nstep, r, r1);
      makeFiles(name, nbp, nstep, r, r1);
      calcWrithe(name, nbp, nstep);
      writeln("job ", name, " done!!! XD!! :) ");
    }
  }
}
