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
  use Map;

  // Command-line configuration -- maps to the ov_NA.in inputs in the Fortran.
  config var configFile: string = "";
  config var elasFile: string = "elastic.out";
  config var strucFile: string = "structural.out";

  // Persistence-length subfragment [rAa..rAb] and sublength [rAl1..rAl2]
  config var rAa: int = 0;
  config var rAb: int = 0;
  config var rAl1: int = 0;
  config var rAl2: int = 0;

  // Twist subfragment and sublength
  config var rTa: int = 0;
  config var rTb: int = 0;
  config var rTl1: int = 0;
  config var rTl2: int = 0;

  // Stretch subfragment and sublength
  config var rSa: int = 0;
  config var rSb: int = 0;
  config var rSl1: int = 0;
  config var rSl2: int = 0;

  /*
   * Parse a simple key = value config file.
   * Lines beginning with 'c' or '#' are treated as comments.
   * Returns a map of key -> value strings.
   */
  proc parseConfigFile(path: string) {
    var cfg: map(string, string);
    try! {
      var f = open(path, ioMode.r);
      var reader = f.reader();
      var line: string;
      while reader.readLine(line) {
        line = line.strip();
        if line == "" then continue;
        if line[0] == "c" || line[0] == "#" then continue;
        var eqIdx = line.find("="):int;
        if eqIdx == -1 then continue;
        var k = line[0..eqIdx-1].strip();
        var v = line[eqIdx+1..].strip();
        cfg.add(k, v);
      }
      f.close();
    }
    return cfg;
  }

  proc main() {
    // Load config file values as defaults (CLI flags override them).
    var cfgElasFile = elasFile,
        cfgStrucFile = strucFile;
    var cfgRaA = rAa, cfgRaB = rAb,
        cfgRaL1 = rAl1, cfgRaL2 = rAl2;
    var cfgRtA = rTa, cfgRtB = rTb,
        cfgRtL1 = rTl1, cfgRtL2 = rTl2;
    var cfgRsA = rSa, cfgRsB = rSb,
        cfgRsL1 = rSl1, cfgRsL2 = rSl2;

    if configFile != "" {
      var cfg = parseConfigFile(configFile);
      if elasFile == "elastic.out" && cfg.contains("elasFile") then
        cfgElasFile = try! cfg["elasFile"];
      if strucFile == "structural.out" && cfg.contains("strucFile") then
        cfgStrucFile = try! cfg["strucFile"];
      if rAa == 0 && cfg.contains("rAa") then
        cfgRaA = try! cfg["rAa"]:int;
      if rAb == 0 && cfg.contains("rAb") then
        cfgRaB = try! cfg["rAb"]:int;
      if rAl1 == 0 && cfg.contains("rAl1") then
        cfgRaL1 = try! cfg["rAl1"]:int;
      if rAl2 == 0 && cfg.contains("rAl2") then
        cfgRaL2 = try! cfg["rAl2"]:int;
      if rTa == 0 && cfg.contains("rTa") then
        cfgRtA = try! cfg["rTa"]:int;
      if rTb == 0 && cfg.contains("rTb") then
        cfgRtB = try! cfg["rTb"]:int;
      if rTl1 == 0 && cfg.contains("rTl1") then
        cfgRtL1 = try! cfg["rTl1"]:int;
      if rTl2 == 0 && cfg.contains("rTl2") then
        cfgRtL2 = try! cfg["rTl2"]:int;
      if rSa == 0 && cfg.contains("rSa") then
        cfgRsA = try! cfg["rSa"]:int;
      if rSb == 0 && cfg.contains("rSb") then
        cfgRsB = try! cfg["rSb"]:int;
      if rSl1 == 0 && cfg.contains("rSl1") then
        cfgRsL1 = try! cfg["rSl1"]:int;
      if rSl2 == 0 && cfg.contains("rSl2") then
        cfgRsL2 = try! cfg["rSl2"]:int;
    }

    // -----------------------------------------------------------------------
    // READING SECTION
    // -----------------------------------------------------------------------
    var enbp, eframes, estr: int;
    var elaspDom = {1..13, 1..0};
    var elasp: [elaspDom] real;
    var ovElaspDom = {1..2, 1..13, 1..0};
    var ovElasp: [ovElaspDom] real;

    var strucpDom = {1..2, 1..11, 1..0};
    var strucp: [strucpDom] real;
    var avstrpDom = {1..3, 1..0};
    var avstrp: [avstrpDom] real;
    var ovStrucpDom = {1..2, 1..11, 1..0};
    var ovStrucp: [ovStrucpDom] real;
    var ovAvstrpDom = {1..2, 1..3, 1..0};
    var ovAvstrp: [ovAvstrpDom] real;

    coforall task in 1..2
        with (ref enbp, ref eframes, ref estr,
              ref elaspDom, ref elasp, ref ovElaspDom, ref ovElasp,
              ref strucpDom, ref strucp, ref avstrpDom, ref avstrp,
              ref ovStrucpDom, ref ovStrucp, ref ovAvstrpDom, ref ovAvstrp) {
      if task == 1 {
        // Read elastic parameters from file.
        var (_, t_enbp, t_eframes, _, t_estr, _, _, _,
             _, _, t_elasp, t_ovElasp, _, _, _, _, _, _) =
          readElasticParms(cfgElasFile);
        enbp = t_enbp;
        eframes = t_eframes;
        estr = t_estr;
        elaspDom = t_elasp.domain;
        elasp = t_elasp;
        ovElaspDom = t_ovElasp.domain;
        ovElasp = t_ovElasp;
      } else {
        // Read structural parameters from file.
        var (_, _, _, _, _, _, _, _,
             _, _, _, _, t_strucp, t_avstrp, t_ovStrucp, t_ovAvstrp, _, _) =
          readStructuralParms(cfgStrucFile);
        strucpDom = t_strucp.domain;
        strucp = t_strucp;
        avstrpDom = t_avstrp.domain;
        avstrp = t_avstrp;
        ovStrucpDom = t_ovStrucp.domain;
        ovStrucp = t_ovStrucp;
        ovAvstrpDom = t_ovAvstrp.domain;
        ovAvstrp = t_ovAvstrp;
      }
    }

    var nbp    = enbp;
    var frames = eframes;
    var str    = estr;

    if elasp.dim(1).size != strucp.dim(2).size then
      halt("Are you sure you are giving me parameters for same system?");

    // -----------------------------------------------------------------------
    // BUILD RANGE ARRAYS
    // -----------------------------------------------------------------------
    var rA: [1..2, 1..2] int;
    rA[1, 1] = cfgRaA;
    rA[2, 1] = cfgRaB;
    rA[1, 2] = cfgRaL1;
    rA[2, 2] = cfgRaL2;

    var rT: [1..2, 1..2] int;
    rT[1, 1] = cfgRtA;
    rT[2, 1] = cfgRtB;
    rT[1, 2] = cfgRtL1;
    rT[2, 2] = cfgRtL2;

    var rS: [1..2, 1..2] int;
    rS[1, 1] = cfgRsA;
    rS[2, 1] = cfgRsB;
    rS[1, 2] = cfgRsL1;
    rS[2, 2] = cfgRsL2;

    // -----------------------------------------------------------------------
    // VALIDATE RANGES
    // -----------------------------------------------------------------------
    var wholeAf = (rA[1, 1] == 0 && rA[2, 1] == 0);
    var wholeAl = (rA[1, 2] == 0 && rA[2, 2] == 0);
    validateSubfrag(rA[1, 1], rA[2, 1], nbp, str, "Persistence length");
    validateSublen(rA[1, 2], rA[2, 2], rA[1, 1], rA[2, 1],
                   nbp, "Persistence length");

    var wholeTf = (rT[1, 1] == 0 && rT[2, 1] == 0);
    var wholeTl = (rT[1, 2] == 0 && rT[2, 2] == 0);
    validateSubfrag(rT[1, 1], rT[2, 1], nbp, str, "Twist");
    validateSublen(rT[1, 2], rT[2, 2], rT[1, 1], rT[2, 1],
                   nbp, "Twist");

    // Stretch subfragment: apply default midpoint selection when not given.
    var wholeSf: bool;
    if rS[1, 1] == 0 && rS[2, 1] == 0 {
      if str == 2 {
        wholeSf = true;
      } else if nbp >= 18 {
        var s = nbp / 2;
        rS[1, 1] = s - 8;
        rS[2, 1] = s + 9;
        wholeSf = false;
      } else {
        wholeSf = true;
      }
    } else {
      wholeSf = false;
    }
    var wholeSl = (rS[1, 2] == 0 && rS[2, 2] == 0);
    if !wholeSf then
      validateSubfrag(rS[1, 1], rS[2, 1], nbp, str, "Stretch");
    validateSublen(rS[1, 2], rS[2, 2], rS[1, 1], rS[2, 1],
                   nbp, "Stretch");

    // -----------------------------------------------------------------------
    // COMPUTE AVERAGING VECTOR SIZES
    // -----------------------------------------------------------------------
    var aN = rangeN(rA[1, 1], rA[2, 1], nbp);
    var tN = rangeN(rT[1, 1], rT[2, 1], nbp);
    var sN = rangeN(rS[1, 1], rS[2, 1], nbp);

    // -----------------------------------------------------------------------
    // REARRANGE DATA INTO 1-D AVERAGING VECTORS
    // -----------------------------------------------------------------------

    // Tilt: elastic param index 4
    var avTilt: [1..aN] real;
    if wholeAf {
      avTilt = ovElasp[1, 4, 1..aN];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[4, i];
      avTilt = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Roll: elastic param index 3
    var avRoll: [1..aN] real;
    if wholeAf {
      avRoll = ovElasp[1, 3, 1..aN];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[3, i];
      avRoll = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Dynamic persistence length from tilt & roll: elastic param index 11
    var avAdC: [1..aN] real;
    if wholeAf {
      avAdC = ovElasp[1, 11, 1..aN];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[11, i];
      avAdC = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Twist: elastic param index 2
    var avTwist: [1..tN] real;
    if wholeTf {
      avTwist = ovElasp[1, 2, 1..tN];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[2, i];
      avTwist = centralFragment(col, rT[1, 1], rT[2, 1], nbp, str);
    }

    // Partial variance of end-to-end distance: elastic param index 13
    var avEE: [1..sN] real;
    if wholeSf {
      avEE = ovElasp[1, 13, 1..sN];
    } else {
      var col: [1..elasp.dim(1).size] real;
      for i in 1..elasp.dim(1).size do col[i] = elasp[13, i];
      avEE = centralFragment(col, rS[1, 1], rS[2, 1], nbp, str);
    }

    // Squared bending angle (strucp index 10):
    //   transform to 1 - angle²/2  [in rad²]
    var avB2: [1..aN] real;
    if wholeAf {
      for i in 1..aN do
        avB2[i] = 1.0 - 0.5 * ovStrucp[1, 10, i]
                             * deg_to_rad * deg_to_rad;
    } else {
      var nBsp = strucp.dim(2).size;
      var col: [1..nBsp] real;
      for i in 1..nBsp do
        col[i] = 1.0 - 0.5 * strucp[1, 10, i] * deg_to_rad * deg_to_rad;
      avB2 = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // Average-structure bend² (avstrp index 2):
    //   transform similarly
    var b2Avstr: [1..aN] real;
    if wholeAf {
      for i in 1..aN do
        b2Avstr[i] = 1.0 - 0.5 * ovAvstrp[1, 2, i]
                                 * deg_to_rad * deg_to_rad;
    } else {
      var nBsp = avstrp.dim(1).size;
      var col: [1..nBsp] real;
      for i in 1..nBsp do
        col[i] = 1.0 - 0.5 * avstrp[2, i] * deg_to_rad * deg_to_rad;
      b2Avstr = centralFragment(col, rA[1, 1], rA[2, 1], nbp, str);
    }

    // -----------------------------------------------------------------------
    // SET SUBLENGTH RANGES WHEN USING WHOLE-FRAGMENT DEFAULTS
    // -----------------------------------------------------------------------
    if wholeAl {
      if aN >= 21 {
        rA[1, 2] = 11;
        rA[2, 2] = aN - 10;
      } else {
        rA[1, 2] = 1;
        rA[2, 2] = aN;
      }
    }

    if wholeTl {
      if tN >= 21 {
        rT[1, 2] = 11;
        rT[2, 2] = tN - 10;
      } else {
        rT[1, 2] = 1;
        rT[2, 2] = tN;
      }
    }

    if wholeSl {
      if sN >= 17 {
        rS[1, 2] = 8;
        rS[2, 2] = 17;
      } else if sN >= 10 {
        rS[1, 2] = 8;
        rS[2, 2] = sN;
      } else {
        rS[1, 2] = 1;
        rS[2, 2] = sN;
      }
    }

    // -----------------------------------------------------------------------
    // CALCULATE ELASTIC CONSTANTS
    // -----------------------------------------------------------------------

    // Tilt
    var lA = rA[2, 2] - rA[1, 2] + 1;
    var Tilt = averageStd(avTilt[rA[1, 2]..rA[2, 2]]);
    Tilt[2] /= sqrt(lA:real);

    // Roll
    var Roll = averageStd(avRoll[rA[1, 2]..rA[2, 2]]);
    Roll[2] /= sqrt(lA:real);

    // Dynamic persistence length (c) from tilt & roll averages
    var adC = averageStd(avAdC[rA[1, 2]..rA[2, 2]]);
    adC[2] /= sqrt(lA:real);

    // Twist
    var lT = rT[2, 2] - rT[1, 2] + 1;
    var Twist = averageStd(avTwist[rT[1, 2]..rT[2, 2]]);
    Twist[2] /= sqrt(lT:real);

    // Stretch via simple linear regression
    var lS = rS[2, 2] - rS[1, 2] + 1;
    var tempS: [1..lS] real;
    for i in 0..lS-1 do tempS[i+1] = (rS[1, 2] + i):real;

    var Stretch: [1..4] real;
    var (_, bS, cS) = simpleLinearRegression(
      avEE[rS[1, 2]..rS[2, 2]], tempS);
    Stretch[3] = bS;
    Stretch[4] = cS;
    Stretch[1] = bnm * 1.0e23 * bKTpN / Stretch[3];
    var rAux = bnm * 1.0e23 * bKTpN * Stretch[4];
    Stretch[4] = abs(rAux / Stretch[3]**2);

    // Persistence lengths — linear fit forced through intercept = 1
    var auxR = rA[1, 2];
    rA[1, 2] = 1;
    var lA2 = rA[2, 2] - rA[1, 2] + 1;
    var tempA: [1..lA2] real;
    for i in 0..lA2-1 do tempA[i+1] = (rA[1, 2] + i):real;

    // (a) total persistence length
    var aA: [1..4] real;
    var (_, bAa, cAa) = simpleLinearRegressionA1(
      avB2[rA[1, 2]..rA[2, 2]], tempA);
    aA[3] = bAa;
    aA[4] = cAa;
    aA[1] = -bnm / aA[3];
    rAux = bnm * aA[4];
    aA[4] = abs(rAux / aA[3]**2);

    // (a) static persistence length
    var asA: [1..4] real;
    var (_, bAsa, cAsa) = simpleLinearRegressionA1(
      b2Avstr[rA[1, 2]..rA[2, 2]], tempA);
    asA[3] = bAsa;
    asA[4] = cAsa;
    asA[1] = -bnm / asA[3];
    rAux = bnm * asA[4];
    asA[4] = abs(rAux / asA[3]**2);

    // Dynamic contribution b2_d = 1 + av_b2 - b2_avstr
    var b2D: [1..aN] real;
    for i in 1..aN do b2D[i] = 1.0 + avB2[i] - b2Avstr[i];

    // (a) dynamic persistence length
    var adA: [1..4] real;
    var (_, bAda, cAda) = simpleLinearRegressionA1(
      b2D[rA[1, 2]..rA[2, 2]], tempA);
    adA[3] = bAda;
    adA[4] = cAda;
    adA[1] = -bnm / adA[3];
    rAux = bnm * adA[4];
    adA[4] = abs(rAux / adA[3]**2);

    // (b) harmonic mean of static and dynamic (a)
    var aBb: [1..2] real;
    aBb[1] = adA[1] * asA[1] / (adA[1] + asA[1]);

    // (d) harmonic mean using static (a) and dynamic (c)
    var aD: [1..2] real;
    aD[1] = adC[1] * asA[1] / (adC[1] + asA[1]);

    rA[1, 2] = auxR; // restore auxR

    // -----------------------------------------------------------------------
    // PRINT RESULTS
    // -----------------------------------------------------------------------
    printElasticConstants(Tilt, Roll, Twist, Stretch,
                          aA, asA, adA, aBb, adC, aD,
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
    else if b < a then return nbp + b - a;
    else
      return nbp - 1; // a == b == 0: whole fragment
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
                               nbp: int, name: string) {
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
                                     aA: [] real, asA: [] real,
                                     adA: [] real, aBb: [] real,
                                     adC: [] real, aD: [] real,
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
    writeln("--------------------------------------------------------");

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
           "A [a] (nm):", aA[1], aA[2], aA[3], aA[4], "#");
    writef("%-15s %15.3f %15.5f %15.5f %15.3f %15s\n",
           "As [a] (nm):", asA[1], asA[2], asA[3], asA[4], "#");
    writef("%-15s %15.3f %15.5f %15.5f %15.3f %15s\n",
           "Ad [a] (nm):", adA[1], adA[2], adA[3], adA[4], "#");
    writef("%-15s %15.3f %15s %15s %15s %15s\n",
           "A [b] (nm):", aBb[1], "#", "#", "#", "#");
    writef("%-15s %15.3f %15s %15s %15s %15.3f\n",
           "Ad [c] (nm):", adC[1], "#", "#", "#", adC[2]);
    writef("%-15s %15.3f %15s %15s %15s %15s\n",
           "A [d] (nm):", aD[1], "#", "#", "#", "#");
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
