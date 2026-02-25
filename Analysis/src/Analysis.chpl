/*
 * Copyright (C) 2026 Alexander Dukhan
 * Credit owed to Victor Velasco and Agnes Noy whose fortran software
 * this is a direct re-write of.
 *
 * SerraNA: Analysis
 * Calculates overall elastic constants and different definitions of
 * persistence lengths. Direct Chapel port of sources/SerraNA/Analysis.f90.
 */
module Analysis {
  use Parms;
  use Functions;
  use CustomIO;
  use Math;
  use IO;

  // Command-line configuration -- maps to the ov_NA.in inputs in the Fortran.
  config const elasFile: string = "elastic.out";
  config const strucFile: string = "structural.out";

  // Persistence-length subfragment [rA_a..rA_b] and sublength [rA_l1..rA_l2]
  config const rA_a:  int = 0;
  config const rA_b:  int = 0;
  config const rA_l1: int = 0;
  config const rA_l2: int = 0;

  // Twist subfragment and sublength
  config const rT_a:  int = 0;
  config const rT_b:  int = 0;
  config const rT_l1: int = 0;
  config const rT_l2: int = 0;

  // Stretch subfragment and sublength
  config const rS_a:  int = 0;
  config const rS_b:  int = 0;
  config const rS_l1: int = 0;
  config const rS_l2: int = 0;

  proc main() {
    // -----------------------------------------------------------------------
    // READING SECTION
    // -----------------------------------------------------------------------
    // Read elastic parameters from file.
    var (_, enbp, eframes, _, estr, _, _, _,
         _, _, elasp, ovElasp, _, _, _, _, _, _) =
      readElasticParms(elasFile);

    // Read structural parameters from file.
    var (_, _, _, _, _, _, _, _,
         _, _, _, _, strucp, avstrp, ovStrucp, ovAvstrp, _, _) =
      readStructuralParms(strucFile);

    var nbp    = enbp;
    var frames = eframes;
    var str    = estr;

    if elasp.dim(1).size != strucp.dim(2).size then
      halt("Are you sure you are giving me parameters for same system?");

    // -----------------------------------------------------------------------
    // BUILD RANGE ARRAYS
    // -----------------------------------------------------------------------
    var rA: [1..2, 1..2] int;
    rA[1, 1] = rA_a;
    rA[2, 1] = rA_b;
    rA[1, 2] = rA_l1;
    rA[2, 2] = rA_l2;

    var rT: [1..2, 1..2] int;
    rT[1, 1] = rT_a;
    rT[2, 1] = rT_b;
    rT[1, 2] = rT_l1;
    rT[2, 2] = rT_l2;

    var rS: [1..2, 1..2] int;
    rS[1, 1] = rS_a;
    rS[2, 1] = rS_b;
    rS[1, 2] = rS_l1;
    rS[2, 2] = rS_l2;

    // -----------------------------------------------------------------------
    // VALIDATE RANGES
    // -----------------------------------------------------------------------
    var whole_Af = (rA[1, 1] == 0 && rA[2, 1] == 0);
    var whole_Al = (rA[1, 2] == 0 && rA[2, 2] == 0);
    validateSubfrag(rA[1, 1], rA[2, 1], nbp, str, "Persistence length");
    validateSublen(rA[1, 2], rA[2, 2], rA[1, 1], rA[2, 1],
                   nbp, str, "Persistence length");

    var whole_Tf = (rT[1, 1] == 0 && rT[2, 1] == 0);
    var whole_Tl = (rT[1, 2] == 0 && rT[2, 2] == 0);
    validateSubfrag(rT[1, 1], rT[2, 1], nbp, str, "Twist");
    validateSublen(rT[1, 2], rT[2, 2], rT[1, 1], rT[2, 1],
                   nbp, str, "Twist");

    // Stretch subfragment: apply default midpoint selection when not given.
    var whole_Sf: bool;
    if rS[1, 1] == 0 && rS[2, 1] == 0 {
      if str == 2 {
        whole_Sf = true;
      } else if nbp >= 18 {
        var s = nbp / 2;
        rS[1, 1] = s - 8;
        rS[2, 1] = s + 9;
        whole_Sf = false;
      } else {
        whole_Sf = true;
      }
    } else {
      whole_Sf = false;
    }
    var whole_Sl = (rS[1, 2] == 0 && rS[2, 2] == 0);
    if !whole_Sf then
      validateSubfrag(rS[1, 1], rS[2, 1], nbp, str, "Stretch");
    validateSublen(rS[1, 2], rS[2, 2], rS[1, 1], rS[2, 1],
                   nbp, str, "Stretch");

    // -----------------------------------------------------------------------
    // COMPUTE AVERAGING VECTOR SIZES
    // -----------------------------------------------------------------------
    var A_n = rangeN(rA[1, 1], rA[2, 1], nbp);
    var T_n = rangeN(rT[1, 1], rT[2, 1], nbp);
    var S_n = rangeN(rS[1, 1], rS[2, 1], nbp);

    // -----------------------------------------------------------------------
    // REARRANGE DATA INTO 1-D AVERAGING VECTORS
    // -----------------------------------------------------------------------

    // Tilt: elastic param index 4
    var av_tilt: [1..A_n] real;
    if whole_Af {
      av_tilt = ovElasp[1, 4, 1..A_n];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[4, i];
      av_tilt = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Roll: elastic param index 3
    var av_roll: [1..A_n] real;
    if whole_Af {
      av_roll = ovElasp[1, 3, 1..A_n];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[3, i];
      av_roll = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Dynamic persistence length from tilt & roll: elastic param index 11
    var av_Ad_c: [1..A_n] real;
    if whole_Af {
      av_Ad_c = ovElasp[1, 11, 1..A_n];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[11, i];
      av_Ad_c = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Twist: elastic param index 2
    var av_twist: [1..T_n] real;
    if whole_Tf {
      av_twist = ovElasp[1, 2, 1..T_n];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[2, i];
      av_twist = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Partial variance of end-to-end distance: elastic param index 13
    var av_E_E: [1..S_n] real;
    if whole_Sf {
      av_E_E = ovElasp[1, 13, 1..S_n];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[13, i];
      av_E_E = centralFragment(col, rS[1, 1], rS[2, 1], nbp, str);
    }

    // Squared bending angle (strucp index 10):
    //   transform to 1 - angle²/2  [in rad²]
    var av_b2: [1..A_n] real;
    if whole_Af {
      for i in 1..A_n do
        av_b2[i] = 1.0 - 0.5 * ovStrucp[1, 10, i]
                             * deg_to_rad * deg_to_rad;
    } else {
      var nBsp = strucp.dim(2).size;
      var col: [1..nBsp] real;
      for i in 1..nBsp do
        col[i] = 1.0 - 0.5 * strucp[1, 10, i] * deg_to_rad * deg_to_rad;
      av_b2 = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Average-structure bend² (avstrp index 2):
    //   transform similarly
    var b2_avstr: [1..A_n] real;
    if whole_Af {
      for i in 1..A_n do
        b2_avstr[i] = 1.0 - 0.5 * ovAvstrp[1, 2, i]
                                 * deg_to_rad * deg_to_rad;
    } else {
      var nBsp = avstrp.dim(1).size;
      var col: [1..nBsp] real;
      for i in 1..nBsp do
        col[i] = 1.0 - 0.5 * avstrp[2, i] * deg_to_rad * deg_to_rad;
      b2_avstr = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // -----------------------------------------------------------------------
    // SET SUBLENGTH RANGES WHEN USING WHOLE-FRAGMENT DEFAULTS
    // -----------------------------------------------------------------------
    if whole_Al {
      if A_n >= 21 {
        rA[1, 2] = 11;
        rA[2, 2] = A_n - 10;
      } else {
        rA[1, 2] = 1;
        rA[2, 2] = A_n;
      }
    }

    if whole_Tl {
      if T_n >= 21 {
        rT[1, 2] = 11;
        rT[2, 2] = T_n - 10;
      } else {
        rT[1, 2] = 1;
        rT[2, 2] = T_n;
      }
    }

    if whole_Sl {
      if S_n >= 17 {
        rS[1, 2] = 8;
        rS[2, 2] = 17;
      } else if S_n >= 10 {
        rS[1, 2] = 8;
        rS[2, 2] = S_n;
      } else {
        rS[1, 2] = 1;
        rS[2, 2] = S_n;
      }
    }

    // -----------------------------------------------------------------------
    // CALCULATE ELASTIC CONSTANTS
    // -----------------------------------------------------------------------

    // Tilt
    var lA = rA[2, 2] - rA[1, 2] + 1;
    var Tilt = averageStd(av_tilt[rA[1, 2]..rA[2, 2]]);
    Tilt[2] /= sqrt(lA:real);

    // Roll
    var Roll = averageStd(av_roll[rA[1, 2]..rA[2, 2]]);
    Roll[2] /= sqrt(lA:real);

    // Dynamic persistence length (c) from tilt & roll averages
    var Ad_c = averageStd(av_Ad_c[rA[1, 2]..rA[2, 2]]);
    Ad_c[2] /= sqrt(lA:real);

    // Twist
    var lT = rT[2, 2] - rT[1, 2] + 1;
    var Twist = averageStd(av_twist[rT[1, 2]..rT[2, 2]]);
    Twist[2] /= sqrt(lT:real);

    // Stretch via simple linear regression
    var lS = rS[2, 2] - rS[1, 2] + 1;
    var tempS: [1..lS] real;
    for i in 0..lS-1 do tempS[i+1] = (rS[1, 2] + i):real;

    var Stretch: [1..4] real;
    var (_, bS, cS) = simpleLinearRegression(
      av_E_E[rS[1, 2]..rS[2, 2]], tempS);
    Stretch[3] = bS;
    Stretch[4] = cS;
    Stretch[1] = bnm * 1.0e23 * bKTpN / Stretch[3];
    var r_aux = bnm * 1.0e23 * bKTpN * Stretch[4];
    Stretch[4] = abs(r_aux / Stretch[3]**2);

    // Persistence lengths — linear fit forced through intercept = 1
    var aux_r = rA[1, 2];
    rA[1, 2] = 1;
    var lA2 = rA[2, 2] - rA[1, 2] + 1;
    var tempA: [1..lA2] real;
    for i in 0..lA2-1 do tempA[i+1] = (rA[1, 2] + i):real;

    // (a) total persistence length
    var A_a: [1..4] real;
    var (_, bAa, cAa) = simpleLinearRegressionA1(
      av_b2[rA[1, 2]..rA[2, 2]], tempA);
    A_a[3] = bAa;
    A_a[4] = cAa;
    A_a[1] = -bnm / A_a[3];
    r_aux = bnm * A_a[4];
    A_a[4] = abs(r_aux / A_a[3]**2);

    // (a) static persistence length
    var As_a: [1..4] real;
    var (_, bAsa, cAsa) = simpleLinearRegressionA1(
      b2_avstr[rA[1, 2]..rA[2, 2]], tempA);
    As_a[3] = bAsa;
    As_a[4] = cAsa;
    As_a[1] = -bnm / As_a[3];
    r_aux = bnm * As_a[4];
    As_a[4] = abs(r_aux / As_a[3]**2);

    // Dynamic contribution b2_d = 1 + av_b2 - b2_avstr
    var b2_d: [1..A_n] real;
    for i in 1..A_n do b2_d[i] = 1.0 + av_b2[i] - b2_avstr[i];

    // (a) dynamic persistence length
    var Ad_a: [1..4] real;
    var (_, bAda, cAda) = simpleLinearRegressionA1(
      b2_d[rA[1, 2]..rA[2, 2]], tempA);
    Ad_a[3] = bAda;
    Ad_a[4] = cAda;
    Ad_a[1] = -bnm / Ad_a[3];
    r_aux = bnm * Ad_a[4];
    Ad_a[4] = abs(r_aux / Ad_a[3]**2);

    // (b) harmonic mean of static and dynamic (a)
    var A_bb: [1..2] real;
    A_bb[1] = Ad_a[1] * As_a[1] / (Ad_a[1] + As_a[1]);

    // (d) harmonic mean using static (a) and dynamic (c)
    var A_d: [1..2] real;
    A_d[1] = Ad_c[1] * As_a[1] / (Ad_c[1] + As_a[1]);

    rA[1, 2] = aux_r; // restore aux_r

    // -----------------------------------------------------------------------
    // PRINT RESULTS
    // -----------------------------------------------------------------------
    printElasticConstants(Tilt, Roll, Twist, Stretch,
                          A_a, As_a, Ad_a, A_bb, Ad_c, A_d,
                          rA, rS, rT, str, frames, nbp);
  }

  // -------------------------------------------------------------------------
  // HELPER PROCS
  // -------------------------------------------------------------------------

  /*
   * Size of the averaging vector for a subfragment [a, b] over nbp base-pairs.
   * Mirrors Analysis.f90's A_n / T_n / S_n calculation.
   */
  private proc rangeN(a: int, b: int, nbp: int): int {
    if a < b then return b - a;
    else if b > a then return nbp + b - a;
    else return nbp - 1; // a == b == 0: whole fragment
  }

  /*
   * Validate a subfragment selection [a, b] against nbp.
   * Halts with a descriptive message if invalid.
   */
  private proc validateSubfrag(a: int, b: int, nbp: int,
                                str: int, name: string) {
    if a == b && b != 0 && a != 0 then
      halt("Invalid subfragment selection for " + name);
    if a < 1 && b != 0 then
      halt("Invalid subfragment selection for " + name);
    if b < 1 && a != 0 then
      halt("Invalid subfragment selection for " + name);
    if a > nbp then halt("Invalid subfragment selection for " + name);
    if b > nbp then halt("Invalid subfragment selection for " + name);
    if str == 1 && b < a then
      halt("Invalid subfragment selection for " + name);
  }

  /*
   * Validate a sublength selection [l1, l2] for a given quantity.
   * Also validates that [l1,l2] fits inside the subfragment [a, b].
   */
  private proc validateSublen(l1: int, l2: int, a: int, b: int,
                               nbp: int, str: int, name: string) {
    if l1 == l2 && l1 != 0 && l2 != 0 then
      halt("Invalid sublength selection for " + name);
    if l1 < 1 && l2 != 0 then
      halt("Invalid sublength selection for " + name);
    if l2 < 1 && l1 != 0 then
      halt("Invalid sublength selection for " + name);
    if l2 < l1 then halt("Invalid sublength selection for " + name);

    if a != 0 && b != 0 && l1 != 0 && l2 != 0 {
      var s = if b < a then b + nbp else b;
      var span = abs(s - a);
      if l1 < 1 || l1 > span then
        halt("Invalid sublength selection for " + name);
      if l2 < 1 || l2 > span then
        halt("Invalid sublength selection for " + name);
    }
  }

  /*
   * Print a formatted summary of all computed elastic constants to stdout.
   * Mirrors io_mod.f90::print_elastic_constants.
   */
  private proc printElasticConstants(Tilt: [] real, Roll: [] real,
                                     Twist: [] real, Stretch: [] real,
                                     A_a: [] real, As_a: [] real,
                                     Ad_a: [] real, A_bb: [] real,
                                     Ad_c: [] real, A_d: [] real,
                                     rA: [] int, rS: [] int, rT: [] int,
                                     str: int, frames: int, nbp: int) {
    if str == 2 then writeln("CLOSED STRUCTURE ANALYSED");
    else writeln("LINEAR STRUCTURE ANALYSED");
    writeln("BASE-PAIRS: ", nbp);
    writeln("FRAMES:     ", frames);
    writeln("");

    writef("%-15s %-15s %-15s %-15s %-15s %-15s\n",
           "", "Elastic cte", "Intercept", "Slope",
           "Confidence-I", "Strd Error");
    writeln("-------------------------------------------------------------------------------------------");

    writef("%-15s %15.3f %15s %15s %15s %15.3f\n",
           "Tilt (nm):", Tilt[1], "#", "#", "#", Tilt[2]);
    writef("%-15s %15.3f %15s %15s %15s %15.3f\n",
           "Roll (nm):", Roll[1], "#", "#", "#", Roll[2]);
    writef("%-15s %15.3f %15s %15s %15s %15.3f\n",
           "Twist (nm):", Twist[1], "#", "#", "#", Twist[2]);
    writef("%-15s %15.3f %15.5f %15.5f %15.3f %15s\n",
           "Stretch (pN):", Stretch[1], Stretch[2], Stretch[3],
           Stretch[4], "#");
    writef("%-15s %15.3f %15.5f %15.5f %15.3f %15s\n",
           "A [a] (nm):", A_a[1], A_a[2], A_a[3], A_a[4], "#");
    writef("%-15s %15.3f %15.5f %15.5f %15.3f %15s\n",
           "As [a] (nm):", As_a[1], As_a[2], As_a[3], As_a[4], "#");
    writef("%-15s %15.3f %15.5f %15.5f %15.3f %15s\n",
           "Ad [a] (nm):", Ad_a[1], Ad_a[2], Ad_a[3], Ad_a[4], "#");
    writef("%-15s %15.3f %15s %15s %15s %15s\n",
           "A [b] (nm):", A_bb[1], "#", "#", "#", "#");
    writef("%-15s %15.3f %15s %15s %15s %15.3f\n",
           "Ad [c] (nm):", Ad_c[1], "#", "#", "#", Ad_c[2]);
    writef("%-15s %15.3f %15s %15s %15s %15s\n",
           "A [d] (nm):", A_d[1], "#", "#", "#", "#");
    writeln("");

    if rA[1, 1] != rA[2, 1] {
      writef("Persistence lengths, Tilt and Roll calculated from base: " +
             "%i to base %i and from lengths: %i to %i\n",
             rA[1, 1], rA[2, 1], rA[1, 2], rA[2, 2]);
    } else {
      writef("Persistence lengths, Tilt and Roll calculated from whole " +
             "fragment and from lengths: %i to %i\n",
             rA[1, 2], rA[2, 2]);
    }

    if rT[1, 1] != rT[2, 1] {
      writef("Twist calculated from base: %i to base %i " +
             "and from lengths: %i to %i\n",
             rT[1, 1], rT[2, 1], rT[1, 2], rT[2, 2]);
    } else {
      writef("Twist calculated from whole fragment and from lengths:" +
             " %i to %i\n", rT[1, 2], rT[2, 2]);
    }

    if rS[1, 1] != rS[2, 1] {
      writef("Stretch calculated from base: %i to base %i " +
             "and from lengths: %i to %i\n",
             rS[1, 1], rS[2, 1], rS[1, 2], rS[2, 2]);
    } else {
      writef("Stretch calculated from whole fragment and from lengths:" +
             " %i to %i\n", rS[1, 2], rS[2, 2]);
    }

    writeln("");
    writeln("A = Persistence length, As = Static persistence length," +
            " Ad = Dynamic persistence length");
    writeln("(a) => overall fitting");
    writeln("(b) => A = As*Ad/(As+Ad)");
    writeln("(c) => Ad = Tilt*Roll/(Tilt+Roll)");
    writeln("(d) => A = As*Ad/(As+Ad) with Ad obtained by (c)");
  }
}
