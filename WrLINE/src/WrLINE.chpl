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
  use LinearAlgebra;
  use IO;
  use FileSystem;

  config const name: string = "";
  config const top: string = "";
  config const traj: string = "";
  config const nbp: int = 0;
  config const nstep: int = 0;
  config const doStrip: bool = false;

  // -----------------------------------------------------------------------
  // MATHEMATICAL AND VECTOR OPERATIONS (from caxislib.py)
  // -----------------------------------------------------------------------

  /*
   * Find arctan of y/x and return values between -pi and pi
   * (equivalent to atan2, but matching the Python source explicitly).
   */
  proc arctan360(x: real, y: real): real {
    var A: real;
    if x >= 0.0 && y >= 0.0 {
      A = atan(y / x);        // Q1
    } else if x < 0.0 && y >= 0.0 {
      A = atan(y / x) + pi;   // Q2
    } else if x < 0.0 && y < 0.0 {
      A = atan(y / x) - pi;   // Q3
    } else {
      A = atan(y / x);        // Q4
    }
    return A;
  }

  /*
   * Return the rotation matrix that rotates vector V to the z axis.
   */
  proc setZ(Vin: [1..3] real): [1..3, 1..3] real {
    var V = Vin / norm(Vin);
    var x = V[1];
    var y = V[2];
    var z = V[3];

    var A = arctan360(z, y);
    var B = atan(-x / sqrt(y**2 + z**2));

    var Rx: [1..3, 1..3] real;
    Rx[1, 1] = 1.0;
    Rx[1, 2] = 0.0;
    Rx[1, 3] = 0.0;
    Rx[2, 1] = 0.0;
    Rx[2, 2] = cos(A);
    Rx[2, 3] = -sin(A);
    Rx[3, 1] = 0.0;
    Rx[3, 2] = sin(A);
    Rx[3, 3] = cos(A);

    var Ry: [1..3, 1..3] real;
    Ry[1, 1] = cos(B);
    Ry[1, 2] = 0.0;
    Ry[1, 3] = sin(B);
    Ry[2, 1] = 0.0;
    Ry[2, 2] = 1.0;
    Ry[2, 3] = 0.0;
    Ry[3, 1] = -sin(B);
    Ry[3, 2] = 0.0;
    Ry[3, 3] = cos(B);

    var Rxy: [1..3, 1..3] real = dot(Ry, Rx);
    return Rxy;
  }

  /*
   * Return the rotation matrix about z that puts V on the x axis.
   */
  proc setX(Vin: [1..3] real): [1..3, 1..3] real {
    var V = Vin / norm(Vin);
    var x = V[1];
    var y = V[2];

    var C = -arctan360(x, y);

    var Rz: [1..3, 1..3] real;
    Rz[1, 1] = cos(C);
    Rz[1, 2] = -sin(C);
    Rz[1, 3] = 0.0;
    Rz[2, 1] = sin(C);
    Rz[2, 2] = cos(C);
    Rz[2, 3] = 0.0;
    Rz[3, 1] = 0.0;
    Rz[3, 2] = 0.0;
    Rz[3, 3] = 1.0;

    return Rz;
  }

  /*
   * Calculate twist angle between two base pairs.
   *  Z is the vector defining the local Z axis.
   *  P11-P12 are a base pair on strand A/B at step j,
   *  P21-P22 are a base pair on strand A/B at step j+1.
   */
  proc twistAngle(P11: [1..3] real, P12: [1..3] real,
                  P21: [1..3] real, P22: [1..3] real,
                  Zin: [1..3] real): real {
    var r1 = P12 - P11;
    var r2 = P22 - P21;
    r1 /= norm(r1);
    r2 /= norm(r2);
    var Z = Zin / norm(Zin);

    // Rotate Z to align with z axis
    var Rxy = setZ(Z);
    var r1r: [1..3] real = dot(Rxy, r1);
    var r2r: [1..3] real = dot(Rxy, r2);

    // Rotate r1 about z so its y-component becomes 0, rotate r2 along
    var Rz = setX(r1r);
    var r2rr: [1..3] real = dot(Rz, r2r);

    return arctan360(r2rr[1], r2rr[2]) * 180.0 / pi;
  }

  /*
   * Element-wise cross product for 3-component column arrays
   * indexed [component] with values for a single vector.
   */
  proc crossVec(A: [1..3] real, B: [1..3] real): [1..3] real {
    var C: [1..3] real;
    C[1] = A[2] * B[3] - A[3] * B[2];
    C[2] = A[3] * B[1] - A[1] * B[3];
    C[3] = A[1] * B[2] - A[2] * B[1];
    return C;
  }

  /*
   * Element-wise cross product for arrays indexed [component, timestep].
   * Returns an array of the same shape.
   */
  proc crossTS(A: [?D] real, B: [D] real): [D] real
    where D.rank == 2 {
    var nsteps = D.dim(1).size;
    var C: [D] real;
    for t in 1..nsteps {
      C[1, t] = A[2, t] * B[3, t] - A[3, t] * B[2, t];
      C[2, t] = A[3, t] * B[1, t] - A[1, t] * B[3, t];
      C[3, t] = A[1, t] * B[2, t] - A[2, t] * B[1, t];
    }
    return C;
  }

  /*
   * Element-wise dot product for arrays indexed [component, timestep].
   * Returns a 1D array indexed [timestep].
   */
  proc dotTS(A: [?D] real, B: [D] real): [D.dim(1)] real
    where D.rank == 2 {
    var nsteps = D.dim(1).size;
    var result: [D.dim(1)] real;
    for t in 1..nsteps {
      result[t] = A[1, t] * B[1, t] + A[2, t] * B[2, t]
                  + A[3, t] * B[3, t];
    }
    return result;
  }

  /*
   * Element-wise magnitude for arrays indexed [component, timestep].
   * Returns a 1D array indexed [timestep].
   */
  proc sizeTS(A: [?D] real): [D.dim(1)] real
    where D.rank == 2 {
    var nsteps = D.dim(1).size;
    var s: [D.dim(1)] real;
    for t in 1..nsteps {
      s[t] = sqrt(A[1, t]**2 + A[2, t]**2 + A[3, t]**2);
    }
    return s;
  }

  // -----------------------------------------------------------------------
  // I/O OPERATIONS (from caxislib.py)
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

      // Read all coordinate values (8 characters wide, fixed format)
      use List;
      var coordList: list(real);
      var line: string;
      while reader.readLine(line) {
        var l = line.size;
        // Strip trailing newline
        if l > 0 && line[l - 1] == "\n" then l -= 1;
        var n = l / 8;
        for i in 0..#n {
          var token = line[i * 8..#8];
          coordList.pushBack(token: real);
        }
      }
      f.close();

      var la = coordList.size;
      var nCoords = la / 3; // total (x,y,z) triplets

      // Reshape into [nCoords, 3] then split by timestep
      // Python: a = reshape(a, (la/3, 3))
      //         x = reshape(a[:,0], (nstep, nbp*2))
      //         y = reshape(a[:,1], (nstep, nbp*2))
      //         z = reshape(a[:,2], (nstep, nbp*2))
      var x: [1..inNstep, 1..inNbp * 2] real;
      var y: [1..inNstep, 1..inNbp * 2] real;
      var z: [1..inNstep, 1..inNbp * 2] real;
      var idx = 0;
      for coord in 0..#nCoords {
        var cx = coordList[idx]; idx += 1;
        var cy = coordList[idx]; idx += 1;
        var cz = coordList[idx]; idx += 1;
        // coord maps to (timestep, atom) in row-major
        var t = coord / (inNbp * 2) + 1;
        var a = coord % (inNbp * 2) + 1;
        x[t, a] = cx;
        y[t, a] = cy;
        z[t, a] = cz;
      }

      // Build strand arrays
      // Python (0-indexed): rA = [x[:,0:nbp:1], y[:,0:nbp:1], z[:,0:nbp:1]]
      //   => columns 0..nbp-1 => Chapel 1..nbp
      // Python (0-indexed): rB = [x[:,2*nbp-1:nbp-1:-1], ...]
      //   => columns (2*nbp-1) down to nbp => Chapel (2*nbp) down to (nbp+1)
      var rA: [1..3, 1..inNstep, 1..inNbp] real;
      var rB: [1..3, 1..inNstep, 1..inNbp] real;
      for t in 1..inNstep {
        for j in 1..inNbp {
          rA[1, t, j] = x[t, j];
          rA[2, t, j] = y[t, j];
          rA[3, t, j] = z[t, j];
          // reversed second strand:
          // Python idx (0-based): 2*nbp-1-k for k=0..nbp-1
          // Chapel idx (1-based): 2*nbp - (j-1) = 2*nbp - j + 1
          var revIdx = 2 * inNbp - j + 1;
          rB[1, t, j] = x[t, revIdx];
          rB[2, t, j] = y[t, revIdx];
          rB[3, t, j] = z[t, revIdx];
        }
      }

      // r = midpoints of neighboring bp steps (circular)
      // Python: r[:,:,j] = 0.25*(rA[:,:,j]+rA[:,:,(j+1)%nbp]
      //                          +rB[:,:,j]+rB[:,:,(j+1)%nbp])
      var r: [1..3, 1..inNstep, 1..inNbp] real;
      for j in 1..inNbp {
        var jNext = (j % inNbp) + 1; // (j+1)%nbp => 1-based circular
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
  // MAIN DATA PROCESSING (from caxislib.py)
  // -----------------------------------------------------------------------

  /*
   * 1st-order helical axis: average position of 2x5 neighbours of
   * C1'-midpoints. Used for twist calculation.
   */
  proc haxis(inNbp: int, inNstep: int,
             r: [1..3, 1..inNstep, 1..inNbp] real) {
    var rC: [1..3, 1..inNstep, 1..inNbp] real;
    for j in 1..inNbp {
      if j % 100 == 0 then
        writeln("now working on basepairs ", j, "s...");
      for t in 1..inNstep {
        var sumVal: [1..3] real;
        for c in 1..3 do sumVal[c] = r[c, t, j];
        var k = 0;
        while k < 5 {
          k += 1;
          // Circular indexing: (j-k)%nbp and (j+k)%nbp
          var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
          var jp = ((j - 1 + k) % inNbp) + 1;
          for c in 1..3 {
            sumVal[c] += r[c, t, jm] + r[c, t, jp];
          }
        }
        for c in 1..3 {
          rC[c, t, j] = sumVal[c] / (2.0 * k + 1.0);
        }
      }
    }
    return rC;
  }

  /*
   * Calculate twist from the C1' coordinates rA, rB and the
   * 1st-order helical axis rC. Saves tw.ser file.
   */
  proc calcTwist(inName: string, inNbp: int, inNstep: int,
                 rA: [1..3, 1..inNstep, 1..inNbp] real,
                 rB: [1..3, 1..inNstep, 1..inNbp] real,
                 rC: [1..3, 1..inNstep, 1..inNbp] real) {
    var tw: [1..inNstep, 1..inNbp] real;
    for j in 1..inNbp {
      if j % 100 == 0 then
        writeln("now working on basepairs ", j, "s...");
      for t in 1..inNstep {
        var jNext = (j % inNbp) + 1;
        var jPrev = ((j - 2 + inNbp) % inNbp) + 1;
        var Z: [1..3] real;
        for c in 1..3 do Z[c] = rC[c, t, jNext] - rC[c, t, jPrev];
        var P11: [1..3] real; for c in 1..3 do P11[c] = rA[c, t, j];
        var P12: [1..3] real; for c in 1..3 do P12[c] = rB[c, t, j];
        var P21: [1..3] real; for c in 1..3 do P21[c] = rA[c, t, jNext];
        var P22: [1..3] real; for c in 1..3 do P22[c] = rB[c, t, jNext];
        tw[t, j] = twistAngle(P11, P12, P21, P22, Z);
      }
    }

    // Write tw.ser
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
   */
  proc caxis(inNbp: int, inNstep: int,
             r: [1..3, 1..inNstep, 1..inNbp] real,
             tw: [1..inNstep, 1..inNbp] real) {
    var r1: [1..3, 1..inNstep, 1..inNbp] real;
    for j in 1..inNbp {
      if j % 100 == 0 then
        writeln("now working on basepairs ", j, "s...");
      for t in 1..inNstep {
        var Tw = tw[t, j];
        var sumVal: [1..3] real;
        for c in 1..3 do sumVal[c] = r[c, t, j];
        var k = 0;
        var prev = Tw; // track twist total before last addition
        while Tw < 360.0 {
          k += 1;
          prev = Tw;
          var jm = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
          var jp = ((j - 1 + k) % inNbp) + 1;
          Tw += tw[t, jm] + tw[t, jp];
          for c in 1..3 {
            sumVal[c] += r[c, t, jm] + r[c, t, jp];
          }
        }
        // Weight: how much of the last added pair to keep
        var w = (360.0 - prev) / (Tw - prev);
        var jmW = ((j - 1 - k) % inNbp + inNbp) % inNbp + 1;
        var jpW = ((j - 1 + k) % inNbp) + 1;
        for c in 1..3 {
          sumVal[c] -= (1.0 - w) * (r[c, t, jmW] + r[c, t, jpW]);
        }
        for c in 1..3 {
          r1[c, t, j] = sumVal[c] / (2.0 * (k:real + w) - 1.0);
        }
      }
    }
    return r1;
  }

  /*
   * Calculate sine of the register angles.
   * Saves sinreg.ser file.
   */
  proc sinreg(inName: string, inNbp: int, inNstep: int,
              r: [1..3, 1..inNstep, 1..inNbp] real,
              r1: [1..3, 1..inNstep, 1..inNbp] real) {
    // SinReg: [1..nstep, 1..nbp+1], first column is time index
    var sinRegArr: [1..inNstep, 0..inNbp] real;
    for j in 1..inNbp {
      if j % 100 == 0 then
        writeln("now working on basepairs ", j, "s...");
      // m = [j-1, j, j+1] % nbp (1-based circular)
      var m0 = ((j - 2 + inNbp) % inNbp) + 1;
      var m1 = j;
      var m2 = (j % inNbp) + 1;

      // v0 = r1[:,:,m1] - r1[:,:,m0]
      // v1 = r1[:,:,m2] - r1[:,:,m1]
      var v0: [1..3, 1..inNstep] real;
      var v1: [1..3, 1..inNstep] real;
      var M: [1..3, 1..inNstep] real; // minor groove vector
      for t in 1..inNstep {
        for c in 1..3 {
          v0[c, t] = r1[c, t, m1] - r1[c, t, m0];
          v1[c, t] = r1[c, t, m2] - r1[c, t, m1];
          M[c, t] = r[c, t, m0] - r1[c, t, m0];
        }
      }

      var crossV1V0 = crossTS(v1, v0);       // plane vector C
      var crossMC = crossTS(M, crossV1V0);
      var sizeM = sizeTS(M);
      var sizeC = sizeTS(crossV1V0);
      var sizeMC = sizeTS(crossMC);

      for t in 1..inNstep {
        sinRegArr[t, j] = sizeMC[t] / (sizeM[t] * sizeC[t]);
      }

      // Determine sign: (+) for minor groove into the circle
      var diff: [1..3, 1..inNstep] real;
      for t in 1..inNstep {
        for c in 1..3 {
          diff[c, t] = v1[c, t] - v0[c, t];
        }
      }
      var dotVal = dotTS(diff, M);
      for t in 1..inNstep {
        if dotVal[t] < 0.0 then
          sinRegArr[t, j] = -sinRegArr[t, j];
      }
    }

    // First column: time scale (0.01 * step)
    for t in 1..inNstep {
      sinRegArr[t, 0] = 0.01 * t: real;
    }

    // Write sinreg.ser
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
   * Create output xyz and 3col files from r (midpoints) and r1 (CAXIS).
   */
  proc makeFiles(inName: string, inNbp: int, inNstep: int,
                 r: [1..3, 1..inNstep, 1..inNbp] real,
                 r1: [1..3, 1..inNstep, 1..inNbp] real) {
    try! {
      // Write xyz files
      var fxyz = open(inName + "/C.xyz", ioMode.cw);
      var wxyz = fxyz.writer();
      var fxyz1 = open(inName + "/C1.xyz", ioMode.cw);
      var wxyz1 = fxyz1.writer();
      for i in 1..inNstep {
        if i % 1000 == 0 then
          writeln(i, "steps has been written for .xyz output");
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

      // Write 3col files
      var fc = open(inName + "/C.3col", ioMode.cw);
      var wc = fc.writer();
      var fc1 = open(inName + "/C1.3col", ioMode.cw);
      var wc1 = fc1.writer();
      for i in 1..inNstep {
        if i % 1000 == 0 then
          writeln(i, "steps has been written for .3col output");
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
  // WRITHE CALCULATION (from writhe.py)
  // -----------------------------------------------------------------------

  /*
   * Read a 3-column coordinate file and split it by timestep.
   * Returns X: [1..nstep, 1..nbp, 1..3].
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
   * Writhe for a single timestep.
   * x: [1..nstep, 1..nbp, 1..3], t: timestep index, l: number of bp.
   * Returns the writhe using the discretized Gauss linking integral.
   */
  proc wr(x: [] real, t: int, l: int): real {
    // Build y: copy frame t and append first point at end (closure)
    var y: [1..l + 1, 1..3] real;
    for j in 1..l {
      for c in 1..3 do y[j, c] = x[t, j, c];
    }
    for c in 1..3 do y[l + 1, c] = x[t, 1, c];

    var Wr = 0.0;
    for j in 1..l {
      for k in 1..j - 1 {
        var tj: [1..3] real; // tangent at j
        var tk: [1..3] real; // tangent at k
        var rjk: [1..3] real; // vector from k to j
        for c in 1..3 {
          tj[c] = y[j + 1, c] - y[j, c];
          tk[c] = y[k + 1, c] - y[k, c];
          rjk[c] = y[j, c] - y[k, c];
        }
        var crossTjTk = crossVec(tj, tk);
        var dotVal = + reduce(rjk * crossTjTk);
        var normRjk = norm(rjk);
        var W = dotVal / (normRjk ** 3.0) / (2.0 * pi);
        Wr += W;
      }
    }
    return Wr;
  }

  /*
   * Main writhe calculation: reads C1.3col, computes writhe for each
   * timestep, and saves to writhe.ser.
   */
  proc calcWrithe(inName: string, inNbp: int, inNstep: int) {
    writeln("Calculating Writhe...");
    var inp = inName + "/C1.3col";
    var X = readThreeCol(inp, inNbp, inNstep);
    var l = inNbp;

    // Calculate writhe for each timestep
    var writheArr: [1..inNstep, 1..2] real;
    for t in 1..inNstep {
      if t % 100 == 0 then
        writeln("now working on steps ", t, "s...");
      writheArr[t, 1] = t: real;
      writheArr[t, 2] = wr(X, t, l);
    }

    // Write writhe.ser
    try! {
      var f = open(inName + "/writhe.ser", ioMode.cw);
      var writer = f.writer();
      for t in 1..inNstep {
        writer.writef("%5i %9.4dr\n", writheArr[t, 1]: int, writheArr[t, 2]);
      }
      writer.close();
    }
  }

  // -----------------------------------------------------------------------
  // MAIN (from WrLINE.py)
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

      writeln("now making output .xyz and .3col files for ", name, "...");
      makeFiles(name, nbp, nstep, r, r1);

      // Writhe calculation
      calcWrithe(name, nbp, nstep);

      writeln("job ", name, " done!!! XD!! :) ");
    }
  }
}
