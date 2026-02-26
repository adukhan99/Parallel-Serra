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

  // -----------------------------------------------------------------------
  // TUPLE-BASED VECTOR MATH (stack-allocated, no heap overhead)
  // -----------------------------------------------------------------------

  type vec3 = 3*real;      // (real, real, real) — stack-allocated
  type mat3 = 3*(3*real);  // 3 row-vectors — stack-allocated

  inline proc mkVec(x: real, y: real, z: real): vec3 {
    return (x, y, z);
  }

  inline proc vecAdd(a: vec3, b: vec3): vec3 {
    return (a(0) + b(0), a(1) + b(1), a(2) + b(2));
  }

  inline proc vecSub(a: vec3, b: vec3): vec3 {
    return (a(0) - b(0), a(1) - b(1), a(2) - b(2));
  }

  inline proc vecScale(a: vec3, s: real): vec3 {
    return (a(0) * s, a(1) * s, a(2) * s);
  }

  inline proc vecDot(a: vec3, b: vec3): real {
    return a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
  }

  inline proc vecNorm(a: vec3): real {
    return sqrt(a(0)**2 + a(1)**2 + a(2)**2);
  }

  inline proc vecNormalize(a: vec3): vec3 {
    var n = vecNorm(a);
    return (a(0) / n, a(1) / n, a(2) / n);
  }

  inline proc vecCross(a: vec3, b: vec3): vec3 {
    return (a(1) * b(2) - a(2) * b(1),
            a(2) * b(0) - a(0) * b(2),
            a(0) * b(1) - a(1) * b(0));
  }

  // Matrix-vector multiply: M * v
  inline proc matVec(m: mat3, v: vec3): vec3 {
    return (vecDot(m(0), v),
            vecDot(m(1), v),
            vecDot(m(2), v));
  }

  // Matrix-matrix multiply: A * B
  inline proc matMul(a: mat3, b: mat3): mat3 {
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
  inline proc arctan360(x: real, y: real): real {
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
   * Return the rotation matrix that rotates vector V to the z axis.
   */
  inline proc setZ(vin: vec3): mat3 {
    var v = vecNormalize(vin);
    var x = v(0), y = v(1), z = v(2);

    var A = arctan360(z, y);
    var B = atan(-x / sqrt(y**2 + z**2));

    var cosA = cos(A), sinA = sin(A);
    var cosB = cos(B), sinB = sin(B);

    var Rx: mat3 = ((1.0, 0.0, 0.0),
                    (0.0, cosA, -sinA),
                    (0.0, sinA, cosA));
    var Ry: mat3 = ((cosB, 0.0, sinB),
                    (0.0, 1.0, 0.0),
                    (-sinB, 0.0, cosB));
    return matMul(Ry, Rx);
  }

  /*
   * Return the rotation matrix about z that puts V on the x axis.
   */
  inline proc setX(vin: vec3): mat3 {
    var v = vecNormalize(vin);
    var x = v(0), y = v(1);
    var C = -arctan360(x, y);
    var cosC = cos(C), sinC = sin(C);

    return ((cosC, -sinC, 0.0),
            (sinC, cosC, 0.0),
            (0.0, 0.0, 1.0));
  }

  /*
   * Calculate twist angle between two base pairs.
   */
  inline proc twistAngle(P11: vec3, P12: vec3,
                         P21: vec3, P22: vec3,
                         Zin: vec3): real {
    var r1 = vecNormalize(vecSub(P12, P11));
    var r2 = vecNormalize(vecSub(P22, P21));
    var Z = vecNormalize(Zin);

    var Rxy = setZ(Z);
    var r1r = matVec(Rxy, r1);
    var r2r = matVec(Rxy, r2);

    var Rz = setX(r1r);
    var r2rr = matVec(Rz, r2r);

    return arctan360(r2rr(0), r2rr(1)) * 180.0 / pi;
  }

  // -----------------------------------------------------------------------
  // I/O OPERATIONS
  // -----------------------------------------------------------------------

  /*
   * Read the AMBER .mdcrd trajectory file and create 3D arrays of
   * atomic coordinates.
   *
   * Returns (rA, rB, r) where:
   *   rA: [1..3, 1..nstep, 1..nbp] -- strand A C1' atoms
   *   rB: [1..3, 1..nstep, 1..nbp] -- strand B C1' atoms (reversed)
   *   r:  [1..3, 1..nstep, 1..nbp] -- midpoints of neighbouring bp steps
   */
  proc readMdcrd(inName: string, inNbp: int, inNstep: int) {
    try! {
      var fPath = inName + "/C.mdcrd";
      var f = open(fPath, ioMode.r);
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
      // atom order per timestep: atom 0..2*nbp-1
      var atomsPerStep = inNbp * 2;
      for t in 1..inNstep {
        for j in 1..inNbp {
          // Strand A: atoms 0..nbp-1 (0-based) → j-1 (0-based)
          var aIdx = ((t - 1) * atomsPerStep + (j - 1)) * 3;
          rA[1, t, j] = coords[aIdx];
          rA[2, t, j] = coords[aIdx + 1];
          rA[3, t, j] = coords[aIdx + 2];

          // Strand B: atoms reversed, Python: 2*nbp-1-k for k=0..nbp-1
          var bAtom = 2 * inNbp - j; // 0-based atom index
          var bIdx = ((t - 1) * atomsPerStep + bAtom) * 3;
          rB[1, t, j] = coords[bIdx];
          rB[2, t, j] = coords[bIdx + 1];
          rB[3, t, j] = coords[bIdx + 2];
        }
      }

      // r = midpoints of neighboring bp steps (circular)
      var r: [1..3, 1..inNstep, 1..inNbp] real;
      for j in 1..inNbp {
        var jNext = (j % inNbp) + 1;
        for t in 1..inNstep {
          for c in 1..3 {
            r[c, t, j] = 0.25 * (rA[c, t, j] + rA[c, t, jNext]
                                 + rB[c, t, j] + rB[c, t, jNext]);
          }
        }
      }

      return (rA, rB, r);
    }
  }

  // -----------------------------------------------------------------------
  // MAIN DATA PROCESSING — parallelized over base pairs
  // -----------------------------------------------------------------------

  /*
   * 1st-order helical axis: average position of 2×5 neighbours of
   * C1'-midpoints. Parallelized over base pairs.
   */
  proc haxis(inNbp: int, inNstep: int,
             ref r: [1..3, 1..inNstep, 1..inNbp] real) {
    var rC: [1..3, 1..inNstep, 1..inNbp] real;
    forall j in 1..inNbp {
      for t in 1..inNstep {
        var sx = r[1, t, j], sy = r[2, t, j], sz = r[3, t, j];
        for k in 1..5 {
          var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
          var jp = ((j - 1 + k) % inNbp) + 1;
          sx += r[1, t, jm] + r[1, t, jp];
          sy += r[2, t, jm] + r[2, t, jp];
          sz += r[3, t, jm] + r[3, t, jp];
        }
        rC[1, t, j] = sx / 11.0;
        rC[2, t, j] = sy / 11.0;
        rC[3, t, j] = sz / 11.0;
      }
    }
    return rC;
  }

  /*
   * Calculate twist from the C1' coordinates rA, rB and the
   * 1st-order helical axis rC. Parallelized over base pairs.
   * Saves tw.ser file.
   */
  proc calcTwist(inName: string, inNbp: int, inNstep: int,
                 ref rA: [1..3, 1..inNstep, 1..inNbp] real,
                 ref rB: [1..3, 1..inNstep, 1..inNbp] real,
                 ref rC: [1..3, 1..inNstep, 1..inNbp] real) {
    var tw: [1..inNstep, 1..inNbp] real;
    forall j in 1..inNbp {
      for t in 1..inNstep {
        var jNext = (j % inNbp) + 1;
        var jPrev = ((j - 2 + inNbp) % inNbp) + 1;
        var Z = mkVec(rC[1, t, jNext] - rC[1, t, jPrev],
                      rC[2, t, jNext] - rC[2, t, jPrev],
                      rC[3, t, jNext] - rC[3, t, jPrev]);
        var P11 = mkVec(rA[1, t, j], rA[2, t, j], rA[3, t, j]);
        var P12 = mkVec(rB[1, t, j], rB[2, t, j], rB[3, t, j]);
        var P21 = mkVec(rA[1, t, jNext], rA[2, t, jNext],
                        rA[3, t, jNext]);
        var P22 = mkVec(rB[1, t, jNext], rB[2, t, jNext],
                        rB[3, t, jNext]);
        tw[t, j] = twistAngle(P11, P12, P21, P22, Z);
      }
    }

    // Write tw.ser (serial I/O)
    try! {
      var f = open(inName + "/tw.ser", ioMode.cw);
      var writer = f.writer();
      for t in 1..inNstep {
        for j in 1..inNbp {
          writer.writef("%8.3dr", tw[t, j]);
        }
        writer.writeln();
      }
      writer.close();
    }

    return tw;
  }

  /*
   * Calculate central helical axis (CAXIS).
   * Running average of each bp with its 2*k neighbours,
   * weighted by the excess base pair beyond 360 degrees of twist.
   * Parallelized over base pairs.
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
          var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
          var jp = ((j - 1 + k) % inNbp) + 1;
          Tw += tw[t, jm] + tw[t, jp];
          sx += r[1, t, jm] + r[1, t, jp];
          sy += r[2, t, jm] + r[2, t, jp];
          sz += r[3, t, jm] + r[3, t, jp];
        }
        var w = (360.0 - prev) / (Tw - prev);
        var jmW = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
        var jpW = ((j - 1 + k) % inNbp) + 1;
        sx -= (1.0 - w) * (r[1, t, jmW] + r[1, t, jpW]);
        sy -= (1.0 - w) * (r[2, t, jmW] + r[2, t, jpW]);
        sz -= (1.0 - w) * (r[3, t, jmW] + r[3, t, jpW]);
        var denom = 2.0 * (k: real + w) - 1.0;
        r1[1, t, j] = sx / denom;
        r1[2, t, j] = sy / denom;
        r1[3, t, j] = sz / denom;
      }
    }
    return r1;
  }

  /*
   * Calculate sine of the register angles.
   * Parallelized over base pairs.
   * Saves sinreg.ser file.
   */
  proc sinreg(inName: string, inNbp: int, inNstep: int,
              ref r: [1..3, 1..inNstep, 1..inNbp] real,
              ref r1: [1..3, 1..inNstep, 1..inNbp] real) {
    var sinRegArr: [1..inNstep, 0..inNbp] real;

    forall j in 1..inNbp {
      var m0 = ((j - 2 + inNbp) % inNbp) + 1;
      var m1 = j;
      var m2 = (j % inNbp) + 1;

      for t in 1..inNstep {
        // v0 = r1[:,t,m1] - r1[:,t,m0]
        var v0 = mkVec(r1[1, t, m1] - r1[1, t, m0],
                       r1[2, t, m1] - r1[2, t, m0],
                       r1[3, t, m1] - r1[3, t, m0]);
        var v1 = mkVec(r1[1, t, m2] - r1[1, t, m1],
                       r1[2, t, m2] - r1[2, t, m1],
                       r1[3, t, m2] - r1[3, t, m1]);
        var M = mkVec(r[1, t, m0] - r1[1, t, m0],
                      r[2, t, m0] - r1[2, t, m0],
                      r[3, t, m0] - r1[3, t, m0]);

        var C = vecCross(v1, v0);
        var MC = vecCross(M, C);
        var sizeM = vecNorm(M);
        var sizeC = vecNorm(C);
        var sizeMC = vecNorm(MC);

        var val = sizeMC / (sizeM * sizeC);

        // Sign: (+) for minor groove into the circle
        var diff = vecSub(v1, v0);
        if vecDot(diff, M) < 0.0 then val = -val;

        sinRegArr[t, j] = val;
      }
    }

    // First column: time scale (serial)
    for t in 1..inNstep {
      sinRegArr[t, 0] = 0.01 * t: real;
    }

    // Write sinreg.ser (serial I/O)
    try! {
      var f = open(inName + "/sinreg.ser", ioMode.cw);
      var writer = f.writer();
      for t in 1..inNstep {
        for j in 0..inNbp {
          writer.writef("%8.3dr", sinRegArr[t, j]);
        }
        writer.writeln();
      }
      writer.close();
    }
  }

  /*
   * Create output xyz and 3col files (serial I/O — no benefit from
   * parallelism for sequential file writes).
   */
  proc makeFiles(inName: string, inNbp: int, inNstep: int,
                 ref r: [1..3, 1..inNstep, 1..inNbp] real,
                 ref r1: [1..3, 1..inNstep, 1..inNbp] real) {
    try! {
      var fxyz = open(inName + "/C.xyz", ioMode.cw);
      var wxyz = fxyz.writer();
      var fxyz1 = open(inName + "/C1.xyz", ioMode.cw);
      var wxyz1 = fxyz1.writer();
      for i in 1..inNstep {
        wxyz.writeln(inNbp);
        wxyz.writeln();
        wxyz1.writeln(inNbp);
        wxyz1.writeln();
        for j in 1..inNbp {
          wxyz.writef("H %8.3dr %8.3dr %8.3dr \n",
                      r[1, i, j], r[2, i, j], r[3, i, j]);
          wxyz1.writef("H %8.3dr %8.3dr %8.3dr \n",
                       r1[1, i, j], r1[2, i, j], r1[3, i, j]);
        }
      }
      wxyz.close();
      wxyz1.close();

      var fc = open(inName + "/C.3col", ioMode.cw);
      var wc = fc.writer();
      var fc1 = open(inName + "/C1.3col", ioMode.cw);
      var wc1 = fc1.writer();
      for i in 1..inNstep {
        for j in 1..inNbp {
          wc.writef("%8.3dr %8.3dr %8.3dr \n",
                    r[1, i, j], r[2, i, j], r[3, i, j]);
          wc1.writef("%8.3dr %8.3dr %8.3dr \n",
                     r1[1, i, j], r1[2, i, j], r1[3, i, j]);
        }
      }
      wc.close();
      wc1.close();
    }
  }

  // -----------------------------------------------------------------------
  // WRITHE CALCULATION — heavily parallelized
  // -----------------------------------------------------------------------

  /*
   * Writhe for a single timestep. The O(n²) double loop is
   * parallelized over j with a reduction on Wr.
   */
  proc wr(ref X: [?D] real, t: int, l: int): real {
    // Build y: copy frame t and append first point (closure)
    var y: [1..l + 1, 1..3] real;
    for j in 1..l {
      for c in 1..3 do y[j, c] = X[t, j, c];
    }
    for c in 1..3 do y[l + 1, c] = X[t, 1, c];

    // Parallel reduction over outer loop
    var Wr: real = 0.0;
    forall j in 1..l with (+ reduce Wr) {
      for k in 1..j - 1 {
        var tjx = y[j + 1, 1] - y[j, 1];
        var tjy = y[j + 1, 2] - y[j, 2];
        var tjz = y[j + 1, 3] - y[j, 3];
        var tkx = y[k + 1, 1] - y[k, 1];
        var tky = y[k + 1, 2] - y[k, 2];
        var tkz = y[k + 1, 3] - y[k, 3];
        var rjkx = y[j, 1] - y[k, 1];
        var rjky = y[j, 2] - y[k, 2];
        var rjkz = y[j, 3] - y[k, 3];

        // cross(tj, tk)
        var cx = tjy * tkz - tjz * tky;
        var cy = tjz * tkx - tjx * tkz;
        var cz = tjx * tky - tjy * tkx;

        var dotVal = rjkx * cx + rjky * cy + rjkz * cz;
        var normRjk = sqrt(rjkx**2 + rjky**2 + rjkz**2);
        Wr += dotVal / (normRjk ** 3.0) / (2.0 * pi);
      }
    }
    return Wr;
  }

  /*
   * Read a 3-column coordinate file and split it by timestep.
   */
  proc readThreeCol(filename: string, inNbp: int,
                    inNstep: int): [1..inNstep, 1..inNbp, 1..3] real {
    var X: [1..inNstep, 1..inNbp, 1..3] real;
    try! {
      var f = open(filename, ioMode.r);
      var reader = f.reader();
      for i in 1..inNstep {
        for j in 1..inNbp {
          reader.read(X[i, j, 1]);
          reader.read(X[i, j, 2]);
          reader.read(X[i, j, 3]);
        }
      }
      f.close();
    }
    return X;
  }

  /*
   * Main writhe calculation: reads C1.3col, computes writhe for each
   * timestep in parallel, and saves to writhe.ser.
   */
  proc calcWrithe(inName: string, inNbp: int, inNstep: int) {
    writeln("Calculating Writhe...");
    var inp = inName + "/C1.3col";
    var X = readThreeCol(inp, inNbp, inNstep);
    var l = inNbp;

    // Calculate writhe for each timestep — embarrassingly parallel
    var writheArr: [1..inNstep, 1..2] real;
    forall t in 1..inNstep {
      writheArr[t, 1] = t: real;
      writheArr[t, 2] = wr(X, t, l);
    }

    // Write writhe.ser (serial I/O)
    try! {
      var f = open(inName + "/writhe.ser", ioMode.cw);
      var writer = f.writer();
      for t in 1..inNstep {
        writer.writef("%5i %9.4dr\n", writheArr[t, 1]: int,
                      writheArr[t, 2]);
      }
      writer.close();
    }
  }

  // -----------------------------------------------------------------------
  // MAIN
  // -----------------------------------------------------------------------

  proc main() {
    try! {
      if name == "" || nbp == 0 || nstep == 0 {
        writeln("Usage: WrLINE --name=<name> --traj=<traj> ",
                "--nbp=<nbp> --nstep=<nstep> [--top=<top>] ",
                "[--doStrip=<bool>]");
        writeln("  name:    project/directory name");
        writeln("  top:     topology file (.prmtop)");
        writeln("  traj:    trajectory file (.mdcrd)");
        writeln("  nbp:     number of base pairs");
        writeln("  nstep:   number of timesteps");
        writeln("  doStrip: run stripC.sh (default false)");
        return;
      }

      // Optionally run the strip script
      if doStrip {
        writeln("Stripping trajectory to get C1' coordinates...");
        use Subprocess;
        var cmd = spawn(["bash", "stripC.sh", name, top, traj]);
        cmd.wait();
      }

      // Create output directory if needed
      if !isDir(name) then mkdir(name, parents=true);

      writeln("start processing ", name);
      writeln("reading files and initialising coordinate arrays...");
      var (rA, rB, r) = readMdcrd(name, nbp, nstep);

      writeln("calculating first order helical axis for ", name, "...");
      var rC = haxis(nbp, nstep, r);

      writeln("calculating twist for ", name, "...");
      var tw = calcTwist(name, nbp, nstep, rA, rB, rC);

      writeln("calculating helical axis for ", name, "...");
      var r1 = caxis(nbp, nstep, r, tw);

      writeln("calculating register angles for ", name, "...");
      sinreg(name, nbp, nstep, r, r1);

      writeln("now making output .xyz and .3col files for ", name,
              "...");
      makeFiles(name, nbp, nstep, r, r1);

      // Writhe calculation
      calcWrithe(name, nbp, nstep);

      writeln("job ", name, " done!!! XD!! :) ");
    }
  }
}
