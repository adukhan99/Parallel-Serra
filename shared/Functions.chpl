/*
 * Copyright (C) 2026 Alexander Dukhan
 * Credit owed to Victor Velasco and Agnes Noy whose fortran software
 * this is a direct re-write of.
 */
module Functions {
  use Parms;
  use Math;
  use LinearAlgebra;
  use IO;
  use CustomIO;

  // -----------------------------------------------------------------------
  // TUPLE-BASED VECTOR MATH (stack-allocated, no heap overhead)
  // -----------------------------------------------------------------------

  proc mkVec3(x: real, y: real, z: real): vec3 {
    return (x, y, z);
  }

  proc vec3Add(a: vec3, b: vec3): vec3 {
    return (a(0) + b(0), a(1) + b(1), a(2) + b(2));
  }

  proc vec3Sub(a: vec3, b: vec3): vec3 {
    return (a(0) - b(0), a(1) - b(1), a(2) - b(2));
  }

  proc vec3Scale(a: vec3, s: real): vec3 {
    return (a(0) * s, a(1) * s, a(2) * s);
  }

  proc vec3Dot(a: vec3, b: vec3): real {
    return a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
  }

  proc vec3Norm(a: vec3): real {
    return sqrt(a(0)**2 + a(1)**2 + a(2)**2);
  }

  proc vec3Norm2(a: vec3): real {
    return a(0)**2 + a(1)**2 + a(2)**2;
  }

  proc vec3Normalize(a: vec3): vec3 {
    var n = vec3Norm(a);
    if n < eps then return a;
    return (a(0) / n, a(1) / n, a(2) / n);
  }

  proc vec3Cross(a: vec3, b: vec3): vec3 {
    return (a(1) * b(2) - a(2) * b(1),
            a(2) * b(0) - a(0) * b(2),
            a(0) * b(1) - a(1) * b(0));
  }

  // Matrix-vector multiply: M * v
  proc mat3Vec(m: mat3, v: vec3): vec3 {
    return (vec3Dot(m(0), v),
            vec3Dot(m(1), v),
            vec3Dot(m(2), v));
  }

  // Matrix-matrix multiply: A * B
  proc mat3Mul(a: mat3, b: mat3): mat3 {
    // Columns of B
    var b0 = mkVec3(b(0)(0), b(1)(0), b(2)(0)),
        b1 = mkVec3(b(0)(1), b(1)(1), b(2)(1)),
        b2 = mkVec3(b(0)(2), b(1)(2), b(2)(2));
    return ((vec3Dot(a(0), b0), vec3Dot(a(0), b1), vec3Dot(a(0), b2)),
            (vec3Dot(a(1), b0), vec3Dot(a(1), b1), vec3Dot(a(1), b2)),
            (vec3Dot(a(2), b0), vec3Dot(a(2), b1), vec3Dot(a(2), b2)));
  }

  proc mat3Transpose(m: mat3): mat3 {
    return ((m(0)(0), m(1)(0), m(2)(0)),
            (m(0)(1), m(1)(1), m(2)(1)),
            (m(0)(2), m(1)(2), m(2)(2)));
  }

  // -----------------------------------------------------------------------
  // GENERAL FUNCTIONS AND SUBROUTINES
  // -----------------------------------------------------------------------

  /* Average and standard deviation of an array */
  proc averageStd(X: [] real) {
    var n = X.size;
    var res: [1..2] real;

    // Get Average
    res[1] = + reduce(X) / n:real;
    var a = res[1];

    // Now Standard deviation
    var s = + reduce ((X - a) ** 2);
    res[2] = sqrt(s / n:real);

    return res;
  }

  /* Absolute value (norm) of a vector or promoted expression */
  proc absv(a) {
    return sqrt(+ reduce(a**2));
  }

  /* Square of squared norm of a vector or promoted expression */
  proc absv2(a) {
    return + reduce(a**2);
  }

  /* Normalize a vector or promoted expression */
  proc normalizeVector(v) {
    var vsum = sqrt(+ reduce(v**2));
    if vsum < eps then return v;
    return v / vsum;
  }

  /* Cross product of two 3D vectors */
  proc crossProduct3(U: [1..3] real, V: [1..3] real) {
    var res: [1..3] real;
    res[1] = U[2]*V[3] - U[3]*V[2];
    res[2] = U[3]*V[1] - U[1]*V[3];
    res[3] = U[1]*V[2] - U[2]*V[1];
    return res;
  }

  /* General rotation matrix R of magnitude 'a' about unit vector 'u' */
  proc generalRotationMatrix(u: [1..3] real, a: real) {
    var R: [1..3, 1..3] real;
    var ca = cos(a),
        sa = sin(a),
        omca = 1.0 - ca;

    R[1,1] = ca + omca * u[1]**2;
    R[1,2] = omca * u[1] * u[2] - u[3] * sa;
    R[1,3] = omca * u[1] * u[3] + u[2] * sa;
    R[2,1] = omca * u[1] * u[2] + u[3] * sa;
    R[2,2] = ca + omca * u[2]**2;
    R[2,3] = omca * u[2] * u[3] - u[1] * sa;
    R[3,1] = omca * u[1] * u[3] - u[2] * sa;
    R[3,2] = omca * u[2] * u[3] + u[1] * sa;
    R[3,3] = ca + omca * u[3]**2;

    return R;
  }

  /* Jacobi diagonalization method for real symmetric matrix M */
  proc diagonalizationJacobi(M: [?D] real) {
    var ndim = D.dim(0).size;
    var S = M;
    var V: [1..ndim, 1..ndim] real;
    for i in 1..ndim do V[i,i] = 1.0;

    var upS = S;
    var cond = 4.0;
    var r = 0;
    const cont = 50;

    while cond > eps && r < cont {
      r += 1;
      var larg = -1.0;
      var i, j: int;
      for k in 1..ndim-1 {
        for l in k+1..ndim {
          if larg < abs(S[k,l]) {
            larg = abs(S[k,l]);
            i = k;
            j = l;
          }
        }
      }

      var th: real;
      if S[i,i] == S[j,j] {
        th = pi / 4.0;
      } else {
        th = 0.5 * atan(2.0 * S[i,j] / (S[j,j] - S[i,i]));
      }

      var co = cos(th),
          si = sin(th);

      var G: [1..ndim, 1..ndim] real;
      for k in 1..ndim do G[k,k] = 1.0;
      G[i,i] = co;
      G[j,j] = co;
      G[i,j] = si;
      G[j,i] = -si;

      V = dot(V, G);
      // S = G^T * S * G (optimized manually in 4x4, here using dot)
      upS = dot(transpose(G), S);
      upS = dot(upS, G);
      S = upS;

      if larg < eps then break;
    }

    return (S, V);
  }

  /* Projection of vector u on unit vector v */
  proc projvU(u: [1..3] real, v: [1..3] real) {
    var a = dot(u, v);
    return a * v;
  }

  /* Non-unit projection helper (same as projvU if v is unit) */
  proc nprojvU(u: [1..3] real, v: [1..3] real) {
    return dot(u, v) * v;
  }

  /* Matrix inversion for 3x3 matrix (tuple-based) */
  proc mat3Inversion(A: mat3): mat3 {
    var det_A = 
      A(0)(0)*A(1)(1)*A(2)(2) +
      A(0)(1)*A(1)(2)*A(2)(0) +
      A(0)(2)*A(1)(0)*A(2)(1) -
      A(2)(0)*A(1)(1)*A(0)(2) -
      A(2)(1)*A(1)(2)*A(0)(0) -
      A(2)(2)*A(1)(0)*A(0)(1);

    if abs(det_A) < eps then halt("0 determinant of matrix A");

    var inv: mat3;
    inv(0) = (A(1)(1)*A(2)(2) - A(1)(2)*A(2)(1),
              A(0)(2)*A(2)(1) - A(0)(1)*A(2)(2),
              A(0)(1)*A(1)(2) - A(0)(2)*A(1)(1));
    inv(1) = (A(1)(2)*A(2)(0) - A(1)(0)*A(2)(2),
              A(0)(0)*A(2)(2) - A(0)(2)*A(2)(0),
              A(0)(2)*A(1)(0) - A(0)(0)*A(1)(2));
    inv(2) = (A(1)(0)*A(2)(1) - A(1)(1)*A(2)(0),
              A(0)(1)*A(2)(0) - A(0)(0)*A(2)(1),
              A(0)(0)*A(1)(1) - A(0)(1)*A(1)(0));

    var detInv = 1.0 / det_A;
    return ((inv(0)(0) * detInv, inv(0)(1) * detInv, inv(0)(2) * detInv),
            (inv(1)(0) * detInv, inv(1)(1) * detInv, inv(1)(2) * detInv),
            (inv(2)(0) * detInv, inv(2)(1) * detInv, inv(2)(2) * detInv));
  }

  /* Matrix inversion for 3x3 matrix */
  proc matrixInversion3X3(A: [1..3, 1..3] real) {
    var tA: mat3 = ((A[1,1], A[1,2], A[1,3]),
                    (A[2,1], A[2,2], A[2,3]),
                    (A[3,1], A[3,2], A[3,3]));
    var res = mat3Inversion(tA);
    var R: [1..3, 1..3] real;
    for i in 1..3 do
      for j in 1..3 do
        R[i,j] = res(i-1)(j-1);
    return R;
  }

  /* Analytic inversion of 4x4 matrix */
  proc inverseMatrixAnalytic4X4(V: [1..4, 1..4] real) {
    var det = 
    V[1,4]*V[2,3]*V[3,2]*V[4,1] - V[1,3]*V[2,4]*V[3,2]*V[4,1] - 
    V[1,4]*V[2,2]*V[3,3]*V[4,1] + V[1,2]*V[2,4]*V[3,3]*V[4,1] + 
    V[1,3]*V[2,2]*V[3,4]*V[4,1] - V[1,2]*V[2,3]*V[3,4]*V[4,1] - 
    V[1,4]*V[2,3]*V[3,1]*V[4,2] + V[1,3]*V[2,4]*V[3,1]*V[4,2] + 
    V[1,4]*V[2,1]*V[3,3]*V[4,2] - V[1,1]*V[2,4]*V[3,3]*V[4,2] - 
    V[1,3]*V[2,1]*V[3,4]*V[4,2] + V[1,1]*V[2,3]*V[3,4]*V[4,2] + 
    V[1,4]*V[2,2]*V[3,1]*V[4,3] - V[1,2]*V[2,4]*V[3,1]*V[4,3] - 
    V[1,4]*V[2,1]*V[3,2]*V[4,3] + V[1,1]*V[2,4]*V[3,2]*V[4,3] + 
    V[1,2]*V[2,1]*V[3,4]*V[4,3] - V[1,1]*V[2,2]*V[3,4]*V[4,3] - 
    V[1,3]*V[2,2]*V[3,1]*V[4,4] + V[1,2]*V[2,3]*V[3,1]*V[4,4] + 
    V[1,3]*V[2,1]*V[3,2]*V[4,4] - V[1,1]*V[2,3]*V[3,2]*V[4,4] - 
    V[1,2]*V[2,1]*V[3,3]*V[4,4] + V[1,1]*V[2,2]*V[3,3]*V[4,4];

    if abs(det) < eps {
      writeln("Error in Inverse of Covariance, DET=0., showing V matrix");
      writeln(V);
      halt();
    }

    var b: [1..4, 1..4] real;
    b[1,1] = V[2,3]*V[3,4]*V[4,2]-V[2,4]*V[3,3]*V[4,2]+V[2,4]*V[3,2]*V[4,3]-
             V[2,2]*V[3,4]*V[4,3]-V[2,3]*V[3,2]*V[4,4]+V[2,2]*V[3,3]*V[4,4];
    b[1,2] = V[1,4]*V[3,3]*V[4,2]-V[1,3]*V[3,4]*V[4,2]-V[1,4]*V[3,2]*V[4,3]+
             V[1,2]*V[3,4]*V[4,3]+V[1,3]*V[3,2]*V[4,4]-V[1,2]*V[3,3]*V[4,4];
    b[1,3] = V[1,3]*V[2,4]*V[4,2]-V[1,4]*V[2,3]*V[4,2]+V[1,4]*V[2,2]*V[4,3]-
             V[1,2]*V[2,4]*V[4,3]-V[1,3]*V[2,2]*V[4,4]+V[1,2]*V[2,3]*V[4,4];
    b[1,4] = V[1,4]*V[2,3]*V[3,2]-V[1,3]*V[2,4]*V[3,2]-V[1,4]*V[2,2]*V[3,3]+
             V[1,2]*V[2,4]*V[3,3]+V[1,3]*V[2,2]*V[3,4]-V[1,2]*V[2,3]*V[3,4];
    b[2,1] = V[2,4]*V[3,3]*V[4,1]-V[2,3]*V[3,4]*V[4,1]-V[2,4]*V[3,1]*V[4,3]+
             V[2,1]*V[3,4]*V[4,3]+V[2,3]*V[3,1]*V[4,4]-V[2,1]*V[3,3]*V[4,4];
    b[2,2] = V[1,3]*V[3,4]*V[4,1]-V[1,4]*V[3,3]*V[4,1]+V[1,4]*V[3,1]*V[4,3]-
             V[1,1]*V[3,4]*V[4,3]-V[1,3]*V[3,1]*V[4,4]+V[1,1]*V[3,3]*V[4,4];
    b[2,3] = V[1,4]*V[2,3]*V[4,1]-V[1,3]*V[2,4]*V[4,1]-V[1,4]*V[2,1]*V[4,3]+
             V[1,1]*V[2,4]*V[4,3]+V[1,3]*V[2,1]*V[4,4]-V[1,1]*V[2,3]*V[4,4];
    b[2,4] = V[1,3]*V[2,4]*V[3,1]-V[1,4]*V[2,3]*V[3,1]+V[1,4]*V[2,1]*V[3,3]-
             V[1,1]*V[2,4]*V[3,3]-V[1,3]*V[2,1]*V[3,4]+V[1,1]*V[2,3]*V[3,4];
    b[3,1] = V[2,2]*V[3,4]*V[4,1]-V[2,4]*V[3,2]*V[4,1]+V[2,4]*V[3,1]*V[4,2]-
             V[2,1]*V[3,4]*V[4,2]-V[2,2]*V[3,1]*V[4,4]+V[2,1]*V[3,2]*V[4,4];
    b[3,2] = V[1,4]*V[3,2]*V[4,1]-V[1,2]*V[3,4]*V[4,1]-V[1,4]*V[3,1]*V[4,2]+
             V[1,1]*V[3,4]*V[4,2]+V[1,2]*V[3,1]*V[4,4]-V[1,1]*V[3,2]*V[4,4];
    b[3,3] = V[1,2]*V[2,4]*V[4,1]-V[1,4]*V[2,2]*V[4,1]+V[1,4]*V[2,1]*V[4,2]-
             V[1,1]*V[2,4]*V[4,2]-V[1,2]*V[2,1]*V[4,4]+V[1,1]*V[2,2]*V[4,4];
    b[3,4] = V[1,4]*V[2,2]*V[3,1]-V[1,2]*V[2,4]*V[3,1]-V[1,4]*V[2,1]*V[3,2]+
             V[1,1]*V[2,4]*V[3,2]+V[1,2]*V[2,1]*V[3,4]-V[1,1]*V[2,2]*V[3,4];
    b[4,1] = V[2,3]*V[3,2]*V[4,1]-V[2,2]*V[3,3]*V[4,1]-V[2,3]*V[3,1]*V[4,2]+
             V[2,1]*V[3,3]*V[4,2]+V[2,2]*V[3,1]*V[4,3]-V[2,1]*V[3,2]*V[4,3];
    b[4,2] = V[1,2]*V[3,3]*V[4,1]-V[1,3]*V[3,2]*V[4,1]+V[1,3]*V[3,1]*V[4,2]-
             V[1,1]*V[3,3]*V[4,2]-V[1,2]*V[3,1]*V[4,3]+V[1,1]*V[3,2]*V[4,3];
    b[4,3] = V[1,3]*V[2,2]*V[4,1]-V[1,2]*V[2,3]*V[4,1]-V[1,3]*V[2,1]*V[4,2]+
             V[1,1]*V[2,3]*V[4,2]+V[1,2]*V[2,1]*V[4,3]-V[1,1]*V[2,2]*V[4,3];
    b[4,4] = V[1,2]*V[2,3]*V[3,1]-V[1,3]*V[2,2]*V[3,1]+V[1,3]*V[2,1]*V[3,2]-
             V[1,1]*V[2,3]*V[3,2]-V[1,2]*V[2,1]*V[3,3]+V[1,1]*V[2,2]*V[3,3];

    return b / det;
  }

  /* Simple linear regression */
  proc simpleLinearRegression(Y: [] real, X: [] real) {
    var n = X.size;
    if Y.size != n then halt("Different sizes of X & Y in linear regression");

    var avx = + reduce(X) / n:real;
    var avy = + reduce(Y) / n:real;

    var vx = 0.0, cxy = 0.0;
    for i in X.domain {
      vx += (X[i] - avx)**2;
      cxy += (X[i] - avx) * (Y[i] - avy);
    }

    var b = cxy / vx;
    var a = avy - b * avx;

    var e2 = 0.0;
    for i in X.domain {
      e2 += (Y[i] - a - b * X[i])**2;
    }

    if n - 2 <= 0 then halt("n not large enough in linear regression");
    var sb = sqrt((e2 / (n - 2)) / vx);
    var tn_2 = tstudentThQuantil(n - 2);
    var c = tn_2 * sb;

    return (a, b, c);
  }

  /* Simple linear regression with intercept a=0 */
  proc simpleLinearRegressionA0(Y: [] real, X: [] real) {
    var n = X.size;
    if Y.size != n then halt("Different sizes of X & Y in linear regression");

    var avx2 = + reduce(X**2) / n:real;
    var avxy = + reduce(X * Y) / n:real;

    var b = avxy / avx2;
    var a = 0.0;

    var avx = + reduce(X) / n:real;
    var vx = + reduce((X - avx)**2);
    var e2 = + reduce((Y - a - b * X)**2);

    if n - 2 <= 0 then halt("n not large enough in linear regression");
    var sb = sqrt((e2 / (n - 2)) / vx);
    var tn_2 = tstudentThQuantil(n - 2);
    var c = tn_2 * sb;

    return (a, b, c);
  }

  /* Simple linear regression with intercept a=1 */
  proc simpleLinearRegressionA1(Y: [] real, X: [] real) {
    var n = X.size;
    if Y.size != n then halt("Different sizes of X & Y in linear regression");

    var avx = + reduce(X) / n:real;
    var avx2 = + reduce(X**2) / n:real;
    var avxy = + reduce(X * Y) / n:real;

    var b = (avxy - avx) / avx2;
    var a = 1.0;

    var vx = + reduce((X - avx)**2);
    var e2 = + reduce((Y - a - b * X)**2);

    if n - 2 <= 0 then halt("n not large enough in linear regression");
    var sb = sqrt((e2 / (n - 2)) / vx);
    var tn_2 = tstudentThQuantil(n - 2);
    var c = tn_2 * sb;

    return (a, b, c);
  }

  /* t-student distribution quantil */
  proc tstudentThQuantil(df: int) {
    if df < 1 then halt("Error in tstudentThQuantil: df < 1");

    var n = if df > tQuantil.dim(1).size then tQuantil.dim(1).size else df;
    
    return tQuantil[1, n];
  }

  // -----------------------------------------------------------------------
  // SerraNA FUNCTIONS
  // -----------------------------------------------------------------------

  /* Jacobi diagonalization method for real symmetric mat4 (tuple-based) */
  proc diagonalizationJacobi4(M: mat4) {
    var S = M;
    var V: mat4 = ((1.0, 0.0, 0.0, 0.0),
                   (0.0, 1.0, 0.0, 0.0),
                   (0.0, 0.0, 1.0, 0.0),
                   (0.0, 0.0, 0.0, 1.0));

    var cond = 4.0;
    var r = 0;
    const cont = 50;

    while cond > eps && r < cont {
      r += 1;
      var larg = -1.0;
      var i, j: int;
      for k in 0..2 {
        for l in k+1..3 {
          if larg < abs(S(k)(l)) {
            larg = abs(S(k)(l));
            i = k;
            j = l;
          }
        }
      }

      var th: real;
      if S(i)(i) == S(j)(j) {
        th = pi / 4.0;
      } else {
        th = 0.5 * atan(2.0 * S(i)(j) / (S(j)(j) - S(i)(i)));
      }

      var co = cos(th),
          si = sin(th);

      // V = V * G (optimized)
      for rIdx in 0..3 {
        var row = V(rIdx);
        var v_i = row(i);
        var v_j = row(j);
        var resRow = row;
        resRow(i) = v_i * co - v_j * si;
        resRow(j) = v_i * si + v_j * co;
        V(rIdx) = resRow;
      }

      // S = G^T * S * G (optimized)
      var newS = S;
      for k in 0..3 {
        if k == i || k == j then continue;
        var s_ki = S(k)(i);
        var s_kj = S(k)(j);
        newS(i)(k) = s_ki * co - s_kj * si;
        newS(k)(i) = newS(i)(k);
        newS(j)(k) = s_ki * si + s_kj * co;
        newS(k)(j) = newS(j)(k);
      }
      var s_ii = S(i)(i);
      var s_jj = S(j)(j);
      var s_ij = S(i)(j);
      newS(i)(i) = s_ii * co**2 + s_jj * si**2 - 2.0 * s_ij * si * co;
      newS(j)(j) = s_ii * si**2 + s_jj * co**2 + 2.0 * s_ij * si * co;
      newS(i)(j) = 0.0;
      newS(j)(i) = 0.0;
      S = newS;

      if larg < eps then break;
    }

    return (S, V);
  }

  /* Orientation matrix R and origin vector O of a base (optimized) */
  proc getRotationROriginO(Ex: [] real, St: [] real) {
    var av_Ex, av_St: vec3 = (0.0, 0.0, 0.0);

    var N_val = Ex.dim(1).size;
    var N_inv = 1.0 / N_val:real;

    for k in 1..N_val {
      av_Ex(0) += Ex[1, k];
      av_Ex(1) += Ex[2, k];
      av_Ex(2) += Ex[3, k];
      av_St(0) += St[1, k];
      av_St(1) += St[2, k];
      av_St(2) += St[3, k];
    }
    av_Ex = vec3Scale(av_Ex, N_inv);
    av_St = vec3Scale(av_St, N_inv);

    var C: mat3 = ((0.0,0.0,0.0),(0.0,0.0,0.0),(0.0,0.0,0.0));
    var N_1_inv = 1.0 / (N_val - 1):real;

    for k in 1..N_val {
      var dEx = mkVec3(Ex[1, k] - av_Ex(0), Ex[2, k] - av_Ex(1),
                       Ex[3, k] - av_Ex(2));
      var dSt = mkVec3(St[1, k] - av_St(0), St[2, k] - av_St(1),
                       St[3, k] - av_St(2));
      // C += dSt * dEx^T
      for i in 0..2 {
        var row = C(i);
        for j in 0..2 {
          row(j) += dSt(i) * dEx(j);
        }
        C(i) = row;
      }
    }
    // Scale C
    for i in 0..2 do C(i) = vec3Scale(C(i), N_1_inv);

    var M = mat3ToMat4M(C);
    var (S, V) = diagonalizationJacobi4(M);
    var q = largestEigenvec4(S, V);
    var tR = getRotationRTuple(q);

    var R: [1..3, 1..3] real;
    for i in 0..2 do
      for j in 0..2 do
        R[i+1, j+1] = tR(i)(j);

    var tO = vec3Sub(av_Ex, mat3Vec(mat3Transpose(tR), av_St));
    var O: [1..3] real = [tO(0), tO(1), tO(2)];

    return (R, O);
  }

  /* Real symmetric matrix M from 3x3 matrix C (tuple-based) */
  proc mat3ToMat4M(C: mat3): mat4 {
    var M: mat4;
    M(0) = (C(0)(0) + C(1)(1) + C(2)(2),
            C(1)(2) - C(2)(1),
            C(2)(0) - C(0)(2),
            C(0)(1) - C(1)(0));
    M(1) = (M(0)(1),
            C(0)(0) - C(1)(1) - C(2)(2),
            C(0)(1) + C(1)(0),
            C(2)(0) + C(0)(2));
    M(2) = (M(0)(2),
            M(1)(2),
            -C(0)(0) + C(1)(1) - C(2)(2),
            C(1)(2) + C(2)(1));
    M(3) = (M(0)(3),
            M(1)(3),
            M(2)(3),
            -C(0)(0) - C(1)(1) + C(2)(2));
    return M;
  }

  /* Real symmetric matrix M from 3x3 matrix C (array-based) */
  proc rsMatrixM(C: [1..3, 1..3] real) {
    var tC: mat3 = ((C[1,1], C[1,2], C[1,3]),
                    (C[2,1], C[2,2], C[2,3]),
                    (C[3,1], C[3,2], C[3,3]));
    var res = mat3ToMat4M(tC);
    var M: [1..4, 1..4] real;
    for i in 1..4 do
      for j in 1..4 do
        M[i,j] = res(i-1)(j-1);
    return M;
  }

  /* Largest eigenvector of matrix M from diagonalization results S and V */
  proc largestEigenvec4(S: mat4, V: mat4): vec4 {
    var larg = eps;
    var j = 0;
    var k = 0;
    for i in 0..3 {
      if S(i)(i) < 0.0 then k += 1;
      if k == 4 then halt("Error: all eigenvalues negative");
      if larg < S(i)(i) {
        larg = S(i)(i);
        j = i;
      }
    }
    return (V(0)(j), V(1)(j), V(2)(j), V(3)(j));
  }

  /* Rotation matrix R from eigenvector q (tuple-based) */
  proc getRotationRTuple(q: vec4): mat3 {
    var R: mat3;
    R(0) = (q(0)**2 + q(1)**2 - q(2)**2 - q(3)**2,
            2.0 * (q(1)*q(2) - q(0)*q(3)),
            2.0 * (q(1)*q(3) + q(0)*q(2)));
    R(1) = (2.0 * (q(2)*q(1) + q(0)*q(3)),
            q(0)**2 - q(1)**2 + q(2)**2 - q(3)**2,
            2.0 * (q(2)*q(3) - q(0)*q(1)));
    R(2) = (2.0 * (q(3)*q(1) - q(0)*q(2)),
            2.0 * (q(3)*q(2) + q(0)*q(1)),
            q(0)**2 - q(1)**2 - q(2)**2 + q(3)**2);
    return R;
  }

  /* Rotation matrix R from eigenvector q (array-based) */
  proc getRotationR(q: [1..4] real) {
    var tq: vec4 = (q[1], q[2], q[3], q[4]);
    var res = getRotationRTuple(tq);
    var R: [1..3, 1..3] real;
    for i in 1..3 do
      for j in 1..3 do
        R[i,j] = res(i-1)(j-1);
    return R;
  }

  /* Rotation matrix R from axis u and angle a (tuple-based) */
  proc mat3Rotation(u: vec3, a: real): mat3 {
    var ca = cos(a);
    var sa = sin(a);
    var omca = 1.0 - ca;
    var R: mat3;
    R(0) = (ca + omca * u(0)**2,
            omca * u(0) * u(1) - u(2) * sa,
            omca * u(0) * u(2) + u(1) * sa);
    R(1) = (omca * u(0) * u(1) + u(2) * sa,
            ca + omca * u(1)**2,
            omca * u(1) * u(2) - u(0) * sa);
    R(2) = (omca * u(0) * u(2) - u(1) * sa,
            omca * u(1) * u(2) + u(0) * sa,
            ca + omca * u(2)**2);
    return R;
  }

  /* Base-pair parameters (optimized) */
  proc basepairParameters(O1in: [1..3] real, O2in: [1..3] real,
                          R1in: [1..3, 1..3] real, R2in: [1..3, 1..3] real) {
    var O1: vec3 = (O1in[1], O1in[2], O1in[3]);
    var O2: vec3 = (O2in[1], O2in[2], O2in[3]);
    var R1: mat3 = ((R1in[1,1], R1in[1,2], R1in[1,3]),
                    (R1in[2,1], R1in[2,2], R1in[2,3]),
                    (R1in[3,1], R1in[3,2], R1in[3,3]));
    var R2: mat3 = ((R2in[1,1], R2in[1,2], R2in[1,3]),
                    (R2in[2,1], R2in[2,2], R2in[2,3]),
                    (R2in[3,1], R2in[3,2], R2in[3,3]));

    var BPP: [1..6] real;
    var Tmbt: [1..3, 1..3] real;
    var Ombt: [1..3] real;

    // Fortran uses column 3 (R(:,3)) for delta/bp
    var R1_3 = mkVec3(R1(0)(2), R1(1)(2), R1(2)(2));
    var R2_3 = mkVec3(R2(0)(2), R2(1)(2), R2(2)(2));
    var dotR2R1 = vec3Dot(R2_3, R1_3);
    var delta = acos(clamp(dotR2R1, -1.0, 1.0));
    
    var bp = vec3Normalize(vec3Cross(R2_3, R1_3));

    var T1 = mat3Mul(mat3Rotation(bp, -0.5 * delta), R1);
    var T2 = mat3Mul(mat3Rotation(bp, 0.5 * delta), R2);

    var tTmbt: mat3;
    for i in 0..2 {
      var row = vec3Scale(vec3Add(T1(i), T2(i)), 0.5);
      tTmbt(i) = row;
    }

    var col0 = vec3Normalize(mkVec3(tTmbt(0)(0), tTmbt(1)(0), tTmbt(2)(0)));
    var col1 = vec3Normalize(mkVec3(tTmbt(0)(1), tTmbt(1)(1), tTmbt(2)(1)));
    var col2 = vec3Normalize(mkVec3(tTmbt(0)(2), tTmbt(1)(2), tTmbt(2)(2)));
    tTmbt(0) = (col0(0), col1(0), col2(0));
    tTmbt(1) = (col0(1), col1(1), col2(1));
    tTmbt(2) = (col0(2), col1(2), col2(2));

    for i in 0..2 do
      for j in 0..2 do
        Tmbt[i+1, j+1] = tTmbt(i)(j);

    var tOmbt = vec3Scale(vec3Add(O1, O2), 0.5);
    Ombt = [tOmbt(0), tOmbt(1), tOmbt(2)];

    // Fortran uses column 2 (index 1) for Opening
    var T1_2 = mkVec3(T1(0)(1), T1(1)(1), T1(2)(1));
    var T2_2 = mkVec3(T2(0)(1), T2(1)(1), T2(2)(1));
    var dotT2T1 = vec3Dot(T2_2, T1_2);
    BPP[6] = acos(clamp(dotT2T1, -1.0, 1.0));
    
    if vec3Dot(vec3Cross(T2_2, T1_2), col2) < 0.0
             then BPP[6] = -BPP[6];

    // psi uses bp and col1 (Tmbt(:,2))
    var dotBpCol1 = vec3Dot(bp, col1);
    var psi = acos(clamp(dotBpCol1, -1.0, 1.0));
    if vec3Dot(vec3Cross(bp, col1), col2) < 0.0
              then psi = -psi;

    BPP[5] = delta * cos(psi);
    BPP[4] = delta * sin(psi);

    var dO = vec3Sub(O1, O2);
    // BPP(1:3)=MATMUL(O1-O2,Tmbt) -> BPP[i] = dot(dO, col_{i-1})
    BPP[1] = vec3Dot(dO, col0);
    BPP[2] = vec3Dot(dO, col1);
    BPP[3] = vec3Dot(dO, col2);

    for i in 4..6 do BPP[i] *= rad_to_deg;

    return (BPP, Tmbt, Ombt);
  }

  /* Base-step parameters (optimized) */
  proc basestepParameters(O1in: [1..3] real, O2in: [1..3] real,
                          R1in: [1..3, 1..3] real, R2in: [1..3, 1..3] real,
                          totaltwist: real) {
    var O1: vec3 = (O1in[1], O1in[2], O1in[3]);
    var O2: vec3 = (O2in[1], O2in[2], O2in[3]);
    var R1: mat3 = ((R1in[1,1], R1in[1,2], R1in[1,3]),
                    (R1in[2,1], R1in[2,2], R1in[2,3]),
                    (R1in[3,1], R1in[3,2], R1in[3,3]));
    var R2: mat3 = ((R2in[1,1], R2in[1,2], R2in[1,3]),
                    (R2in[2,1], R2in[2,2], R2in[2,3]),
                    (R2in[3,1], R2in[3,2], R2in[3,3]));

    var BSP: [1..9] real;

    // column 3 is index 2
    var R1_3 = mkVec3(R1(0)(2), R1(1)(2), R1(2)(2));
    var R2_3 = mkVec3(R2(0)(2), R2(1)(2), R2(2)(2));

    BSP[8] = clamp(vec3Dot(R1_3, R2_3), -1.0, 1.0);
    BSP[7] = acos(BSP[8]);

    var rt = vec3Normalize(vec3Cross(R1_3, R2_3));
    
    // Fortran: T1 = MATMUL(Rrt, R1) with Rrt(rt, 0.5*BSP(7))
    //          T2 = MATMUL(Rrt, R2) with Rrt(rt, -0.5*BSP(7))
    var T1 = mat3Mul(mat3Rotation(rt, 0.5 * BSP[7]), R1);
    var T2 = mat3Mul(mat3Rotation(rt, -0.5 * BSP[7]), R2);

    // Tmst = 0.5 * (T1 + T2) then normalize columns
    var col0 = vec3Normalize(vec3Scale(
      vec3Add(mkVec3(T1(0)(0), T1(1)(0), T1(2)(0)),
              mkVec3(T2(0)(0), T2(1)(0), T2(2)(0))), 0.5));
    var col1 = vec3Normalize(vec3Scale(
      vec3Add(mkVec3(T1(0)(1), T1(1)(1), T1(2)(1)),
              mkVec3(T2(0)(1), T2(1)(1), T2(2)(1))), 0.5));
    var col2 = vec3Normalize(mkVec3(0.5*(T1(0)(2)+T2(0)(2)), 
                                    0.5*(T1(1)(2)+T2(1)(2)), 
                                    0.5*(T1(2)(2)+T2(2)(2))));
    col2 = vec3Normalize(col2);

    var T1_2 = mkVec3(T1(0)(1), T1(1)(1), T1(2)(1));
    var T2_2 = mkVec3(T2(0)(1), T2(1)(1), T2(2)(1));
    BSP[6] = acos(clamp(vec3Dot(T1_2, T2_2), -1.0, 1.0));
    if vec3Dot(vec3Cross(T1_2, T2_2), col2) < 0.0 then BSP[6] = -BSP[6];

    var nbpInt = (totaltwist / 180.0):int;
    if nbpInt % 2 != 0 then BSP[6] += 2.0 * pi;

    nbpInt = (totaltwist / 360.0):int;
    BSP[6] += nbpInt * 2.0 * pi;

    if totaltwist - rad_to_deg * BSP[6] > 180.0 then BSP[6] += 2.0 * pi;
    if totaltwist - rad_to_deg * BSP[6] < -180.0 then BSP[6] -= 2.0 * pi;

    nbpInt = ((BSP[6] + pi) / (2.0 * pi)):int;
    var sign = if nbpInt % 2 != 0 then -1.0 else 1.0;
    col0 = vec3Scale(col0, sign);
    col1 = vec3Scale(col1, sign);

    var phi = acos(clamp(vec3Dot(rt, col1), -1.0, 1.0));
    if vec3Dot(vec3Cross(rt, col1), col2) < 0.0 then phi = -phi;

    BSP[4] = BSP[7] * sin(phi);
    BSP[5] = BSP[7] * cos(phi);

    var dO = vec3Sub(O2, O1);
    BSP[1] = vec3Dot(dO, col0);
    BSP[2] = vec3Dot(dO, col1);
    BSP[3] = vec3Dot(dO, col2);

    for i in 4..7 do BSP[i] *= rad_to_deg;
    BSP[9] = BSP[7]**2;

    return BSP;
  }

  /* Covariance matrix of deformation variables */
  proc deformationCovariance(Roll: [] real,
          Tilt: [] real,
          Twist: [] real,
          Streatch: [] real) {
    var n = Roll.size;
    var X: [1..n, 1..4] real;
    X[1..n, 1] = Roll;
    X[1..n, 2] = Tilt;
    X[1..n, 3] = Twist;
    X[1..n, 4] = Streatch;

    var av: [1..4] real;
    for j in 1..4 do av[j] = + reduce(X[1..n, j]) / n:real;

    var Y: [1..n, 1..4] real;
    for j in 1..4 do Y[1..n, j] = X[1..n, j] - av[j];

    var V: [1..4, 1..4] real;
    for i in 1..4 {
      for j in 1..4 {
        V[i,j] = + reduce(Y[1..n, i] * Y[1..n, j]) / n:real;
      }
    }
    return V;
  }

  /* Dynamic persistence length */
  proc dynamicPersistenceLength2(tilt: real, roll: real) {
    return 2.0 * (tilt * roll / (tilt + roll));
  }

  /* Rebuilding algorithm for base triads */
  proc reverseAlgorithmBaseTriads(BSP: [] real) {
    var nbp = BSP.dim(1).size;
    var Tg_i: [1..3, 1..3, 1..nbp] real;
    var rgI: [1..3, 1..nbp] real;

    var x: [1..3] real = [1.0, 0.0, 0.0];
    var y: [1..3] real = [0.0, 1.0, 0.0];
    var z: [1..3] real = [0.0, 0.0, 1.0];

    rgI[1..3, 1] = 0.0;
    Tg_i[1..3, 1, 1] = x;
    Tg_i[1..3, 2, 1] = y;
    Tg_i[1..3, 3, 1] = z;

    for i in 2..nbp {
      var an12 = sqrt(BSP[4, i-1]**2 + BSP[5, i-1]**2) * deg_to_rad;
      var ax12: [1..3] real = [BSP[4, i-1]*deg_to_rad/an12, 
               BSP[5, i-1]*deg_to_rad/an12, 0.0];
      var phi = acos(ax12[2]);
      var an3 = BSP[6, i-1] * deg_to_rad;
      var h = crossProduct3(ax12, y);
      if dot(h, z) > 0.0 then phi = abs(phi); else phi = -abs(phi);

      var R1 = generalRotationMatrix(z, an3/2.0 - phi);
      var R2 = generalRotationMatrix(y, an12);
      var R3 = generalRotationMatrix(z, an3/2.0 + phi);
      Tg_i[1..3, 1..3, i] = dot(dot(R1, R2), R3);

      var R2_mid = generalRotationMatrix(y, an12/2.0);
      var R3_mid = generalRotationMatrix(z, phi);
      var Ti_mst = dot(dot(R1, R2_mid), R3_mid);

      var hPos: [1..3] real = [BSP[1, i-1], BSP[2, i-1], BSP[3, i-1]];
      rgI[1..3, i] = dot(hPos, transpose(Ti_mst));

      Tg_i[1..3, 1..3, i] = dot(Tg_i[1..3, 1..3, i-1], Tg_i[1..3, 1..3, i]);
      rgI[1..3, i] = 
        rgI[1..3, i-1] + dot(
          rgI[1..3, i],
          transpose(reshape(Tg_i[1..3, 1..3, i-1], {1..3, 1..3}))
        );
    }
    return (Tg_i, rgI);
  }

  /* Central fragment calculation */
  proc centralFragment(dat: [] real, n1: int, n2: int, nbp: int, str: int) {
    if (n2 < n1) && (str != 2) then 
      halt("Second range less than first makes no sense");
    if n2 > nbp || n1 > nbp then 
      halt("Ranges cannot be greater than number of bases");

    var overallSize = if str == 2 then 
      if n1 < n2 then n2-n1 else nbp-n1+n2 else n2-n1;
    var overall: [1..overallSize] real;
    var l = 1, p = 1;

    if str != 2 {
      for j in 1..nbp-1 {
        for i in 1..nbp-j {
          var s = i + j;
          if i >= n1 && i <= n2 && s <= n2 {
            overall[l] += dat[p];
            if s == n2 then l += 1;
          }
          p += 1;
        }
      }
      for i in 1..n2-n1 do overall[i] /= (n2-n1-i+1):real;
    } else {
      if n1 < n2 {
        for j in 1..nbp-1 {
          for i in 1..nbp {
            var s = i + j;
            if i >= n1 && i <= n2 && s <= n2 {
              overall[l] += dat[p];
              if s == n2 then l += 1;
            }
            p += 1;
          }
        }
        for i in 1..n2-n1 do overall[i] /= (n2-n1-i+1):real;
      } else {
        for j in 1..nbp-1 {
          for i in 1..nbp {
            var s = i + j;
            if s > nbp then s -= nbp;
            if i >= n1 && i >= n2 && s <= n2 {
              overall[l] += dat[p];
              if s == n2 then l += 1;
            }
            p += 1;
          }
        }
        var sRange = nbp - n1 + n2;
        for i in 1..sRange do overall[i] /= (sRange-i+1):real;
      }
    }
    return overall;
  }

  /* Extract sublength from 3D data */
  proc extractSublength3D(datain: [?D] real,
                          nbp: int,
                          sublength: int,
                          str: int) {
    var subdataSize = if str == 2 then nbp else nbp - sublength;
    var subdata: [1..D.dim(0).size, 1..D.dim(1).size, 1..subdataSize] real;

    if str == 1 {
      var l = 0;
      for j in 1..nbp-1 {
        if j == sublength {
          for k in 1..nbp-sublength {
            l += 1;
            subdata[1..D.dim(0).size,
              1..D.dim(1).size, k] = datain[1..D.dim(0).size,
              1..D.dim(1).size, l];
          }
          break;
        }
        l += nbp - j;
      }
    } else if str == 2 {
      var l = nbp * (sublength - 1);
      for k in 1..nbp {
        l += 1;
        subdata[1..D.dim(0).size,
          1..D.dim(1).size, k] = datain[1..D.dim(0).size,
          1..D.dim(1).size, l];
      }
    } else {
      halt("Error in extractSublength3D: Unidentified type of structure");
    }

    return subdata;
  }

  // -----------------------------------------------------------------------
  // SerraLINE FUNCTIONS
  // -----------------------------------------------------------------------

  /* Project coordinates onto a global plane */
  proc projectCoordinatesGPlane(ref coords: [] real,
          nFitting: int,
          bpFitting: [] int) {
    var nbp = coords.dim(1).size;
    var numFrames = coords.dim(2).size;
    var G_n: [1..numFrames, 1..3, 1..3] real;
    var best: [1..numFrames] int;
    var avgDistFrame: [1..numFrames] real;
    var maxDistFrame: [1..numFrames] real;

    forall k in 1..numFrames {
      var c: [1..3] real;
      var A: [1..3, 1..nbp] real;
      
      if bpFitting.size > 0 {
        for i in 1..nFitting do c += coords[1..3, bpFitting[i], k];
        c /= nFitting:real;
        for i in 1..nbp do A[1..3, i] = coords[k, i, 1..3] - c;
        
        var A_fit: [1..3, 1..nFitting] real;
        for i in 1..nFitting do A_fit[1..3, i] = A[1..3, bpFitting[i]];
        
        var (S, U) = diagonalizationJacobi(dot(A_fit, transpose(A_fit)));
        var d: [1..3] real = [S[1,1], S[2,2], S[3,3]];
        var minVal = min reduce d;
        for i in 1..3 do if d[i] == minVal then best[k] = i;

        var distToPlane: [1..nbp] real;
        for i in 1..nbp {
          var A_p = projvU(A[1..3, i], U[1..3, best[k]]);
          coords[k, i, 1..3] = A[1..3, i] - A_p;
          distToPlane[i] = absv(A_p);
          avgDistFrame[k] += distToPlane[i];
        }
        avgDistFrame[k] /= nbp:real;
        maxDistFrame[k] = max reduce distToPlane;
        G_n[k, 1..3, 1..3] = U;
      } else {
        for i in 1..nbp do c += coords[k, i, 1..3];
        c /= nbp:real;
        for i in 1..nbp do A[1..3, i] = coords[k, i, 1..3] - c;
        
        var (S, U) = diagonalizationJacobi(dot(A, transpose(A)));
        var d: [1..3] real = [S[1,1], S[2,2], S[3,3]];
        var minVal = min reduce d;
        for i in 1..3 do if d[i] == minVal then best[k] = i;

        var distToPlane: [1..nbp] real;
        for i in 1..nbp {
          var A_p = projvU(A[1..3, i], U[1..3, best[k]]);
          coords[k, i, 1..3] = A[1..3, i] - A_p;
          distToPlane[i] = absv(A_p);
          avgDistFrame[k] += distToPlane[i];
        }
        avgDistFrame[k] /= nbp:real;
        maxDistFrame[k] = max reduce distToPlane;
        G_n[k, 1..3, 1..3] = U;
      }
    }
    return (G_n, best, avgDistFrame, maxDistFrame);
  }

  /* Calculate tangent vectors */
  proc getTangentVectors(ref coords: [] real,
          circleStr: bool,
          tLength: int) {
    var nbp = coords.dim(1).size;
    var numFrames = coords.dim(2).size;
    var tangents_n = if circleStr then nbp else nbp - tLength;
    var tangents: [1..numFrames, 1..tangents_n, 1..3] real;
    forall k in 1..numFrames {
      for i in 1..nbp-tLength {
        tangents[k, i, 1..3] = 
          normalizeVector(coords[k, i+tLength, 1..3] - coords[k, i, 1..3]);
      }
      if circleStr {
        var j = 0;
        for i in nbp-tLength+1..nbp {
          j += 1;
          tangents[k, i, 1..3] = 
            normalizeVector(coords[k, j, 1..3] - coords[k, i, 1..3]);
        }
      }
    }
    return tangents;
  }

  /* Central fragment calculation with mean and std */
  proc centralFragmentMeanStd(dat: [] real,
          n1: int,
          n2: int,
          nbp: int,
          str: int) {
    var overall = centralFragment(dat, n1, n2, nbp, str);
    var std: [1..overall.size] real;
    var l = 1, p = 1;

    if str != 2 {
      for j in 1..nbp-1 {
        for i in 1..nbp-j {
          var s = i + j;
          if i >= n1 && i <= n2 && s <= n2 {
            std[l] += (dat[p] - overall[l])**2;
            if s == n2 then l += 1;
          }
          p += 1;
        }
      }
      for i in 1..n2-n1 do std[i] = sqrt(std[i] / (n2-n1-i+1):real);
    } else {
      if n1 < n2 {
        for j in 1..nbp-1 {
          for i in 1..nbp {
            var s = i + j;
            if i >= n1 && i <= n2 && s <= n2 {
              std[l] += (dat[p] - overall[l])**2;
              if s == n2 then l += 1;
            }
            p += 1;
          }
        }
        for i in 1..n2-n1 do std[i] = sqrt(std[i] / (n2-n1-i+1):real);
      } else {
        for j in 1..nbp-1 {
          for i in 1..nbp {
            var s = i + j;
            if s > nbp then s -= nbp;
            if i >= n1 && i >= n2 && s <= n2 {
              std[l] += (dat[p] - overall[l])**2;
              if s == n2 then l += 1;
            }
            p += 1;
          }
        }
        var sRange = nbp - n1 + n2;
        for i in 1..sRange do std[i] = sqrt(std[i] / (sRange-i+1):real);
      }
    }
    return (overall, std);
  }

  /* Extract sublength from 2D data */
  proc extractSublength2D(datain: [?D] real,
          nbp: int,
          sublength: int,
          str: int) {
    var subdataSize = if str == 2 then nbp else nbp - sublength;
    var subdata: [1..D.dim(0).size, 1..subdataSize] real;

    if str == 1 {
      var l = 0;
      for j in 1..nbp-1 {
        if j == sublength {
          for k in 1..nbp-sublength {
            l += 1;
            subdata[1..D.dim(0).size, k] = datain[1..D.dim(0).size, l];
          }
          break;
        }
        l += nbp - j;
      }
    } else if str == 2 {
      var l = nbp * (sublength - 1);
      for k in 1..nbp {
        l += 1;
        subdata[1..D.dim(0).size, k] = datain[1..D.dim(0).size, l];
      }
    } else {
      halt("Error in extractSublength2D: Unidentified type of structure");
    }

    return subdata;
  }

  /* Calculate widths, heights and aspect ratios */
  proc getWidthHeightAratio(ndim: int,
          mdim: int,
          numFrames: int,
          ref G_n: [1..3,
          1..3,
          1..numFrames] real,
          ref best: [1..numFrames] int,
          ref coords: [] real) {
    var width: [1..numFrames] real;
    var height: [1..numFrames] real;
    var aratio: [1..numFrames] real;

    forall k in 1..numFrames {
      var dOp = -1.0;
      var dR: [1..3] real;

      for l in 1..ndim-1 {
        for i in 1..ndim-l {
          var j = i + l;
          var r = coords[k, j, 1..3] - coords[k, i, 1..3];
          var d = absv2(r);
          if d > dOp {
            dOp = d;
            dR = r;
          }
        }
      }

      var rUnit = normalizeVector(dR);
      var G_k: [1..3, 1..2] real;
      G_k[1..3, 1] = rUnit;
      G_k[1..3, 2] = crossProduct3(rUnit, G_n[1..3, best[k], k]);

      var awi: [1..mdim] real;
      var ahe: [1..mdim] real;
      var s = 0;
      for l in 1..ndim-1 {
        for i in 1..ndim-l {
          s += 1;
          var j = i + l;
          var rVec = coords[k, j, 1..3] - coords[k, i, 1..3];
          var xy1 = projvU(rVec, G_k[1..3, 1]);
          var xy2 = projvU(rVec, G_k[1..3, 2]);
          ahe[s] = absv(xy1);
          awi[s] = absv(xy2);
        }
      }
      width[k] = max reduce awi;
      height[k] = max reduce ahe;
      aratio[k] = width[k] / height[k];
    }
    return (width, height, aratio);
  }

  /* Calculate widths, heights and aspect ratios for closed structures */
  proc getCWidthHeightAratio(ndim: int,
          numFrames: int,
          ref G_n: [1..3,
          1..3,
          1..numFrames] real,
          ref best: [1..numFrames] int,
          ref coords: [] real) {
    var width: [1..numFrames] real;
    var height: [1..numFrames] real;
    var aratio: [1..numFrames] real;

    forall k in 1..numFrames {
      var dOp = -1.0;
      var dR: [1..3] real;

      for l in 1..ndim-1 {
        for i in 1..ndim {
          var s = i + l;
          var j = if s > ndim then s - ndim else s;
          var r = coords[k, j, 1..3] - coords[k, i, 1..3];
          var d = absv2(r);
          if d > dOp {
            dOp = d;
            dR = r;
          }
        }
      }

      var rUnit = normalizeVector(dR);
      var G_k: [1..3, 1..2] real;
      G_k[1..3, 1] = rUnit;
      G_k[1..3, 2] = crossProduct3(rUnit, G_n[1..3, best[k], k]);

      var awi: [1..ndim, 1..ndim-1] real;
      var ahe: [1..ndim, 1..ndim-1] real;
      for l in 1..ndim-1 {
        for i in 1..ndim {
          var sIdx = i + l;
          var j = if sIdx > ndim then sIdx - ndim else sIdx;
          var rVec = coords[k, j, 1..3] - coords[k, i, 1..3];
          var xy1 = projvU(rVec, G_k[1..3, 1]);
          var xy2 = projvU(rVec, G_k[1..3, 2]);
          ahe[i,l] = absv(xy1);
          awi[i,l] = absv(xy2);
        }
      }
      width[k] = max reduce awi;
      height[k] = max reduce ahe;
      aratio[k] = width[k] / height[k];
    }
    return (width, height, aratio);
  }

  /* Calculate bending angles between tangent vectors */
  proc getBendings(ndim: int,
          numFrames: int,
          ref G_n: [1..3,
          1..3,
          1..numFrames] real,
          ref best: [1..numFrames] int,
          ref tangents: [] real,
          bends_size: int) {
    var bends: [1..numFrames, 1..bends_size] real;
    forall k in 1..numFrames {
      var l = 0;
      var jVal = 1;
      for i in 1..ndim-jVal {
        l += 1;
        var c = 
          crossProduct3(tangents[k, i, 1..3], tangents[k, i+jVal, 1..3]);
        var d = dot(tangents[k, i, 1..3], tangents[k, i+jVal, 1..3]);
        var b: real;
        if d > 1.0 then b = 0.0;
        else if d < -1.0 then b = pi;
          else b = acos(d);

        if dot(c, G_n[1..3, best[k], k]) < 0.0 then b = -b;
        bends[k,l] = rad_to_deg * b;
      }

      for j in 2..ndim-1 {
        for i in 1..ndim-j {
          l += 1;
          bends[k,l] = + reduce(bends[k, i..i+j]);
        }
      }
    }
    for k in 1..numFrames do 
      for i in 1..bends_size do bends[k,i] = abs(bends[k,i]);
    return bends;
  }

  /* Calculate bending angles without fitting a plane */
  proc getBendingsNofit(ndim: int,
          numFrames: int,
          ref tangents: [] real,
          bends_size: int) {
    var bends: [1..numFrames, 1..bends_size] real;
    forall k in 1..numFrames {
      var l = 0;
      for j in 1..ndim-1 {
        for i in 1..ndim-j {
          l += 1;
          var a = dot(tangents[k, i, 1..3], tangents[k, i+j, 1..3]);
          var b = acos(clamp(a, -1.0, 1.0));
          bends[k,l] = rad_to_deg * b;
        }
      }
    }
    return bends;
  }

  /* Calculate bending angles for closed structures */
  proc getCBendings(ndim: int,
          numFrames: int,
          ref G_n: [1..3,
          1..3,
          1..numFrames] real,
          ref best: [1..numFrames] int,
          ref tangents: [] real) {
    var bends: [1..numFrames, 1..ndim, 1..ndim-1] real;
    forall k in 1..numFrames {
      var l = 1;
      for i in 1..ndim-l {
        var j = i + l;
        var c = crossProduct3(tangents[k, i, 1..3], tangents[k, j, 1..3]);
        var d = dot(tangents[k, i, 1..3], tangents[k, j, 1..3]);
        var b: real;
        if d > 1.0 then b = 0.0;
        else if d < -1.0 then b = pi;
          else b = acos(d);

        if dot(c, G_n[1..3, best[k], k]) < 0.0 then b = -b;
        bends[k,i,l] = rad_to_deg * b;
      }
      for i in ndim-l+1..ndim {
        var j = i + l - ndim;
        var c = crossProduct3(tangents[k, i, 1..3], tangents[k, j, 1..3]);
        var d = dot(tangents[k, i, 1..3], tangents[k, j, 1..3]);
        var b: real;
        if d > 1.0 then b = 0.0;
        else if d < -1.0 then b = pi;
          else b = acos(d);

        if dot(c, G_n[1..3, best[k], k]) < 0.0 then b = -b;
        bends[k,i,l] = rad_to_deg * b;
      }

      for l_val in 2..ndim-1 {
        for i in 1..ndim-l_val {
          bends[k, i, l_val] = + reduce(bends[k, i..i+l_val-1, 1]);
        }
        for i in ndim-l_val+1..ndim {
          var j = i + l_val - ndim;
          bends[k, i, l_val] = 
            + reduce(bends[k, i..ndim, 1]) + + reduce(bends[k, 1..j-1, 1]);
        }
      }
    }
    for k in 1..numFrames do 
      for i in 1..ndim do 
        for j in 1..ndim-1 do 
          bends[k,i,j] = abs(bends[k,i,j]);
    return bends;
  }

  /* Calculate bending angles for closed structures without fitting a plane */
  proc getCBendingsNofit(ndim: int,
          numFrames: int,
          ref tangents: [] real) {
    var bends: [1..numFrames, 1..ndim, 1..ndim-1] real;
    forall k in 1..numFrames {
      for l in 1..ndim-1 {
        for i in 1..ndim-l {
          var a = dot(tangents[k, i, 1..3], tangents[k, i+l, 1..3]);
          bends[k,i,l] = rad_to_deg * acos(clamp(a, -1.0, 1.0));
        }
        for i in ndim-l+1..ndim {
          var a = dot(tangents[k, i, 1..3], tangents[k, i+l-ndim, 1..3]);
          bends[k,i,l] = rad_to_deg * acos(clamp(a, -1.0, 1.0));
        }
      }
    }
    return bends;
  }

  /* Write projected coordinates */
  proc writeProjectedCoords(coords: [] real,
          G_n: [] real,
          best: [] int,
          seq: [] string,
          xyz: bool) {
    try! {
      var nbp = coords.dim(1).size;
      var numFrames = coords.dim(2).size;
      if xyz {
        var f = open("projected.xyz", ioMode.cw);
        var writer = f.writer();
        for k in 1..numFrames {
          writer.writeln(nbp);
          writer.writeln("   ");
          var b = G_n[1..3, best[k], k];
          var e3_arr: [1..3] real = [0.0, 0.0, 1.0];
          var angle = acos(clamp(dot(b, e3_arr), -1.0, 1.0));
          var d = normalizeVector(crossProduct3(b, e3_arr));
          var Rot = generalRotationMatrix(d, angle);
          for i in 1..nbp {
            var pos = dot(Rot, coords[k, i, 1..3]);
            writer.writeln(seq[i], " ", pos[1], " ", pos[2], " ", pos[3]);
          }        }
        writer.close();
        f.close();
      } else {
        var f = open("projected.crd", ioMode.cw);
        var writer = f.writer();
        writer.writeln("Generated by SerraLINE");
        for k in 1..numFrames {
          var b = G_n[1..3, best[k], k];
          var e3_arr: [1..3] real = [0.0, 0.0, 1.0];
          var angle = acos(clamp(dot(b, e3_arr), -1.0, 1.0));
          var d = normalizeVector(crossProduct3(b, e3_arr));
          var Rot = generalRotationMatrix(d, angle);
          var rotatedCoords = dot(Rot, coords[k, 1..nbp, 1..3]);
          
          var count = 0;
          for i in 1..nbp {
             var pos = dot(Rot, coords[k, i, 1..3]);
             for j in 1..3 {
               writer.write(pos[j], " ");
               count += 1;
               if count % 10 == 0 then writer.writeln();
             }
           }
          if count % 10 != 0 then writer.writeln();
        }
        writer.close();
        f.close();
      }
    }
  }

  /* Fit a plane to a set of points */
  proc getGPlane(points: [] real) {
    var s = points.dim(1).size;
    var At: [1..3, 1..s] real;
    At[1, 1..s] = points[1, 1..s];
    At[2, 1..s] = points[2, 1..s];
    At[3, 1..s] = 1.0;
    var B = points[3, 1..s];
    
    var C = dot(At, transpose(At));
    var I_C = matrixInversion3X3(C);
    var D = dot(At, B);
    var n = dot(I_C, D);
    return normalizeVector(n);
  }

  /* Helper to clamp values */
  proc clamp(x: real, lo: real, hi: real) {
    if x < lo then return lo;
    if x > hi then return hi;
    return x;
  }
}
