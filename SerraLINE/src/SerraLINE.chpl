/* This is a re-write of the SerraLINE package developed
   by Victor Manuel Velasco Berrelleza and Agnes Noy. All Chapel
   code is written by Alexander Dukhan. */

module SerraLINE {
  use Parms;
  use CustomIO;
  use Functions;
  use Math;
  use IO;
  use Map;
  use List;

  config var configFile: string = "";
  config var top: string = "";
  config var traj: string = "";
  config var isCircular: bool = false; // circle_str
  // str: 1 for single, 2 for double, 0 for none
  config var strandsType: int = 1;
  config var nbpArg: int = 0; // if no topology, specify nbp
  config var fitplane: bool = true;
  // comma-separated indices or "0" for all
  config var bpFitting: string = "0";
  config var tLength: int = 1;
  config var printProj: bool = false;
  config var xyzFormat: bool = false;

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

  proc writeAvParms(nbp: int, nldim: int, frames: int, ref seq: [] string,
                    ref avstd_bends: [] real, ref avstd_width: [] real,
                    ref avstd_height: [] real, ref avstd_aratio: [] real,
                    t_length: int, ref avg_dist: [] real, ref max_dist: [] real,
                    ref avg_dist_rel: [] real, ref max_dist_rel: [] real) {
    try! {
      var f = open("SerraLINE.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln("PARAMETERS\n");
      writer.writeln("OPENED STRUCTURE");
      writer.writeln("METHOD: PROJECTION");
      writer.writeln("BASE PAIRS ", nbp);
      writer.write("SEQUENCE\n");
      for s in seq do writer.write(s);
      writer.writeln();
      writer.writeln("SNAPSHOTS ANALYSED    ", frames);
      writer.writeln("TANGENT LENGTH:       ", t_length);
      writer.writeln();
      writer.writeln(
        "First column averages, second column standard deviations");
      writer.writeln("WIDTH (Angstroms): ",
                     avstd_width[1], " ", avstd_width[2]);
      writer.writeln("HEIGHT (Angstroms): ",
                     avstd_height[1], " ", avstd_height[2]);
      writer.writeln("ASPECT RATIO: ",
                     avstd_aratio[1], " ", avstd_aratio[2]);
      writer.writeln("AVERAGE OF DISTANCES TO PLANE (Angstroms): ",
                     avg_dist[1], " ", avg_dist[2]);
      writer.writeln("AVERAGE OF MAXIMUM DISTANCES TO PLANE (Angstroms): ",
                     max_dist[1], " ", max_dist[2]);
      writer.writeln("RELATIVE AVERAGE OF DISTANCES TO PLANE (percent): ",
                     avg_dist_rel[1], " ", avg_dist_rel[2]);
      writer.writeln(
        "RELATIVE AVERAGE OF MAXIMUM DISTANCES TO PLANE (percent): ",
        max_dist_rel[1], " ", max_dist_rel[2]);
      writer.writeln("\n");
      var l = 0;
      for j in 1..nldim-1 {
        writer.writeln(j+1, "mer");
        writer.writeln("base-step", " ", "Bending angle");
        writer.writeln("-" * 38);
        var endLimit = nldim - j;
        if endLimit >= 1 {
          for i in 1..endLimit {
            l += 1;
            writer.writeln(i, "-", i+j, " ",
                           seq[i], " ", seq[i+j], " ",
                           avstd_bends[1, l], " ", avstd_bends[2, l]);
          }
        }
        writer.writeln();
      }
      writer.close();
    }
  }

  proc writeCAvParms(nbp: int, frames: int, ref seq: [] string,
                     ref avstd_bends: [] real, ref avstd_width: [] real,
                     ref avstd_height: [] real, ref avstd_aratio: [] real,
                     t_length: int, ref avg_dist: [] real,
                     ref max_dist: [] real,
                     ref avg_dist_rel: [] real, ref max_dist_rel: [] real) {
    try! {
      var c_seq: [1..2*nbp] string;
      c_seq[1..nbp] = seq;
      c_seq[nbp+1..2*nbp] = seq;
      var f = open("SerraLINE.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln("PARAMETERS\n");
      writer.writeln("CLOSED STRUCTURE");
      writer.writeln("METHOD: PROJECTION");
      writer.writeln("BASE PAIRS ", nbp);
      writer.write("SEQUENCE\n");
      for s in seq do writer.write(s);
      writer.writeln();
      writer.writeln("SNAPSHOTS ANALYSED    ", frames);
      writer.writeln("TANGENT LENGTH:       ", t_length);
      writer.writeln();
      writer.writeln(
        "First column averages, second column standard deviations");
      writer.writeln("WIDTH (Angstroms): ",
                     avstd_width[1], " ", avstd_width[2]);
      writer.writeln("HEIGHT (Angstroms): ",
                     avstd_height[1], " ", avstd_height[2]);
      writer.writeln("ASPECT RATIO: ",
                     avstd_aratio[1], " ", avstd_aratio[2]);
      writer.writeln("AVERAGE OF DISTANCES TO PLANE (Angstroms): ",
                     avg_dist[1], " ", avg_dist[2]);
      writer.writeln("AVERAGE OF MAXIMUM DISTANCES TO PLANE (Angstroms): ",
                     max_dist[1], " ", max_dist[2]);
      writer.writeln("RELATIVE AVERAGE OF DISTANCES TO PLANE (percent): ",
                     avg_dist_rel[1], " ", avg_dist_rel[2]);
      writer.writeln(
        "RELATIVE AVERAGE OF MAXIMUM DISTANCES TO PLANE (percent): ",
        max_dist_rel[1], " ", max_dist_rel[2]);
      writer.writeln("\n");
      for l in 1..nbp-1 {
        writer.writeln(l+1, "mer");
        writer.writeln("base-step", " ", "Bending angle");
        writer.writeln("-" * 38);
        for i in 1..nbp {
          var s = i + l;
          if s > nbp then s = s - nbp;
          writer.writeln(i, "-", s, " ",
                         c_seq[i], " ", c_seq[s], " ",
                         avstd_bends[1, i, l], " ", avstd_bends[2, i, l]);
        }
        writer.writeln();
      }
      writer.close();
    }
  }

  proc writeAvBends(nbp: int, nldim: int, frames: int, ref seq: [] string,
                    ref avstd_bends: [] real, t_length: int) {
    try! {
      var f = open("SerraLINE.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln(" PARAMETERS");
      writer.writeln(" ");
      writer.writeln(" OPENED STRUCTURE");
      writer.writeln(" METHOD: WITHOUT PROJECTION");
      writer.writef(" BASE PAIRS %12i\n", nbp);
      writer.writeln(" SEQUENCE");
      writer.write(" ");
      for s in seq do writer.write(s);
      writer.writeln();
      writer.writef("  SNAPSHOTS ANALYSED%10i\n", frames);
      writer.writef("     TANGENT LENGTH:%10i\n", t_length);
      writer.writeln(" ");
      writer.writeln(
        " First column averages, second column standard deviations");
      var l = 0;
      for j in 1..nldim-1 {
        writer.writef("%10i", j);
        writer.writeln("        bp");
        writer.writeln("           base-step       Bending angle");
        writer.writeln(" " + "-" * 38);
        var endLimit = nldim - j;
        if endLimit >= 1 {
          for i in 1..endLimit {
            l += 1;
            writer.writef("%4i-%4i  %s%s  %10.3dr%10.3dr\n",
                          i, i+j, seq[i], seq[i+j],
                          avstd_bends[1, l], avstd_bends[2, l]);
          }
        }
        writer.writeln(" ");
      }
      writer.close();
    }
  }

  proc writeCAvBends(nbp: int, frames: int, ref seq: [] string,
                     ref avstd_bends: [] real, t_length: int) {
    try! {
      var c_seq: [1..2*nbp] string;
      c_seq[1..nbp] = seq;
      c_seq[nbp+1..2*nbp] = seq;
      var f = open("SerraLINE.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln("PARAMETERS\n");
      writer.writeln("CLOSED STRUCTURE");
      writer.writeln("METHOD: WITHOUT PROJECTION");
      writer.writeln("BASE PAIRS ", nbp);
      writer.write("SEQUENCE\n");
      for s in seq do writer.write(s);
      writer.writeln();
      writer.writeln("SNAPSHOTS ANALYSED    ", frames);
      writer.writeln("TANGENT LENGTH:       ", t_length);
      writer.writeln();
      writer.writeln(
        "First column averages, second column standard deviations\n");
      for l in 1..nbp-1 {
        writer.writeln(l+1, "mer");
        writer.writeln("base-step", " ", "Bending angle");
        writer.writeln("-" * 38);
        for i in 1..nbp {
          var s = i + l;
          if s > nbp then s = s - nbp;
          writer.writeln(i, "-", s, " ",
                         c_seq[i], " ", c_seq[s], " ",
                         avstd_bends[1, i, l], " ", avstd_bends[2, i, l]);
        }
        writer.writeln();
      }
      writer.close();
    }
  }

  proc main() {
    try! {
      // Load config file values as defaults (CLI flags override them).
      var cfgTop = top,
        cfgTraj = traj,
        cfgIsCircular = isCircular,
        cfgStrandsType = strandsType,
        cfgNbpArg = nbpArg,
        cfgFitplane = fitplane,
        cfgBpFitting = bpFitting,
        cfgTLength = tLength,
        cfgPrintProj = printProj,
        cfgXyzFormat = xyzFormat;
  
      if configFile != "" {
        var cfg = parseConfigFile(configFile);
        // Config values serve as defaults; CLI flags always win.
        if top == "" then
          cfgTop = cfg.get("top", cfgTop);
        if traj == "" then
          cfgTraj = cfg.get("traj", cfgTraj);
        if !isCircular && cfg.contains("isCircular") then
          cfgIsCircular = cfg["isCircular"] == "true";
        if strandsType == 1 && cfg.contains("strandsType") then
          cfgStrandsType = cfg["strandsType"]:int;
        if nbpArg == 0 && cfg.contains("nbpArg") then
          cfgNbpArg = cfg["nbpArg"]:int;
        if fitplane && cfg.contains("fitplane") then
          cfgFitplane = cfg["fitplane"] == "true";
        if bpFitting == "0" then
          cfgBpFitting = cfg.get("bpFitting", cfgBpFitting);
        if tLength == 1 && cfg.contains("tLength") then
          cfgTLength = cfg["tLength"]:int;
        if !printProj && cfg.contains("printProj") then
          cfgPrintProj = cfg["printProj"] == "true";
        if !xyzFormat && cfg.contains("xyzFormat") then
          cfgXyzFormat = cfg["xyzFormat"] == "true";
      }
  
      if cfgTraj == "" {
        writeln(
          "Trajectory file must be specified using --traj=<file>" +
          " or via --configFile=<file>");
        return;
      }
  
      // Use named domains so that we can reassign them to empty ranges
      // to release memory, mirroring Fortran's deallocate() calls.
      var actualNbp = cfgNbpArg;
      var seqDom = {1..0};
      var finalSeq: [seqDom] string;
      var coordDom = {1..3, 1..0, 1..0};
      var coords: [coordDom] real;
      var numFrames = 0;
  
      if cfgStrandsType != 0 && cfgTop != "" {
        writeln("Reading topology file");
        var (origNbp, nAtoms, box, origSeq, ringIndices) =
          topologyAmber(cfgTop, cfgStrandsType);
  
        writeln("Reading trajectory file");
        var atomsToRead = if cfgTraj.endsWith(".xyz") 
                          then origNbp else nAtoms;
        var (rawCoords, frames) = coordinatesOther(cfgTraj, atomsToRead);
        actualNbp = origNbp;
        numFrames = frames;
        coordDom = rawCoords.domain;
        coords = rawCoords;
        seqDom = {1..actualNbp};
        finalSeq = origSeq[0..actualNbp-1];
      } else {
        writeln("Reading trajectory file (no topology)");
        if actualNbp <= 0 then
          halt("Must specify --nbpArg since no topology is provided");
        var (rawCoords, fn) = coordinatesOther(cfgTraj, actualNbp);
        numFrames = fn;
        coordDom = rawCoords.domain;
        coords = rawCoords;
        var arr: [1..actualNbp] string;
        for i in 1..actualNbp do arr[i] = "X";
        seqDom = {1..actualNbp};
        finalSeq = arr;
      }
  
      if !cfgIsCircular && cfgTLength > actualNbp - 2 then
        halt("Tangent length must be less than N-2 for opened structures");
      if cfgIsCircular && cfgTLength > actualNbp - 1 then
        halt("Tangent length must be less than N-1 for closed structures");
  
      var ndim = actualNbp - 1;
      var mdim = ndim * (ndim - 1) / 2;
      var nldim = actualNbp - cfgTLength;
      var mldim = nldim * (nldim - 1) / 2;
      if cfgIsCircular { mdim = actualNbp; }
  
      // Parse bpFitting: "0" or empty means use all base-pairs;
      // otherwise it is a comma-separated list of 1-based indices.
      var bpFittingDom = {1..0};
      var bpFittingArr: [bpFittingDom] int;
      var nFitting = 0;
      if cfgBpFitting != "0" && cfgBpFitting != "" {
        var idxList: list(int);
        for tok in cfgBpFitting.split(",") {
          var idx = tok.strip():int;
          if idx >= 1 && idx <= actualNbp then idxList.pushBack(idx);
        }
        if idxList.size > 0 {
          bpFittingDom = {1..idxList.size};
          for i in 1..idxList.size do bpFittingArr[i] = idxList[i-1];
          nFitting = bpFittingDom.size;
        }
      }
  
      // ---------------------------------------------------------------
      // PLANE SECTION
      // Fortran: allocate G_n, best, avg/max_dist_frame, project,
      //          then deallocate(bp_fitting).
      // Each array gets its own named domain so it can be freed
      // independently, mirroring Fortran's per-variable deallocate().
      // ---------------------------------------------------------------
      var gnDom = {1..3, 1..3, 1..numFrames};
      var G_n: [gnDom] real;
      var bestDom = {1..numFrames};
      var best: [bestDom] int;
      var distDom = {1..numFrames};
      var avgDistFrame, maxDistFrame: [distDom] real;
  
      if cfgFitplane {
        writeln("Fitting planes");
        var res = projectCoordinatesGPlane(
          coords, nFitting, bpFittingArr);
        G_n = res(0);
        best = res(1);
        avgDistFrame = res(2);
        maxDistFrame = res(3);
        if cfgPrintProj then
          writeProjectedCoords(
            coords, G_n, best, finalSeq, cfgXyzFormat);
      }
      // Fortran: deallocate(bp_fitting)
      bpFittingDom = {1..0};
  
      // ---------------------------------------------------------------
      // WIDTH, HEIGHT AND ASPECT RATIO SECTION
      // Fortran: allocate width/height/aratio, compute.
      // ---------------------------------------------------------------
      var widthDom = {1..numFrames};
      var width: [widthDom] real;
      var heightDom = {1..numFrames};
      var height: [heightDom] real;
      var aratioDom = {1..numFrames};
      var aspectRatio: [aratioDom] real;
  
      if cfgFitplane {
        writeln("Calculating widths, heights and aspect ratios");
        if cfgIsCircular {
          (width, height, aspectRatio) =
            getCWidthHeightAratio(mdim, numFrames, G_n, best, coords);
        } else {
          (width, height, aspectRatio) =
            getWidthHeightAratio(ndim, mdim, numFrames,
                                G_n, best, coords);
        }
      }
  
      // ---------------------------------------------------------------
      // TANGENT VECTORS SECTION
      // Fortran: compute tangents, then deallocate(coords).
      // ---------------------------------------------------------------
      writeln("Calculating tangent vectors");
      var tangentsN = if cfgIsCircular then actualNbp else nldim;
      var tangentsDom = {1..3, 1..tangentsN, 1..numFrames};
      var tangents: [tangentsDom] real;
      tangents = getTangentVectors(coords, cfgIsCircular, cfgTLength);
  
      // Fortran line 311-313: deallocate(coords) — no longer needed
      coordDom = {1..3, 1..0, 1..0};
  
      // ---------------------------------------------------------------
      // BENDINGS + RELATIVE SIZES + AVERAGES + WRITING
      //
      // The circular and non-circular paths use completely different
      // array shapes (c_bends vs bends). In the Fortran, these are
      // conditionally allocated. In Chapel, we use block scoping so
      // only the needed array is ever instantiated. This avoids the
      // catastrophic O(F*N^4) c_bends allocation in the non-circular
      // case (and vice-versa).
      //
      // The relative sizes, averaging, and writing stages are placed
      // inside the same branch so they can reference the scoped arrays
      // and so all per-frame data is freed when the branch scope exits.
      // This mirrors the Fortran's deallocate() ordering exactly.
      // ---------------------------------------------------------------
      writeln("Calculating bending angles");
  
      if cfgIsCircular {
        // c_bends: [1..F, 1..nbp, 1..nbp-1] — O(F*N^2), reasonable
        var cBends: [1..numFrames, 1..mdim, 1..mdim-1] real;
        if cfgFitplane {
          cBends = getCBendings(mdim, numFrames, G_n, best, tangents);
        } else {
          cBends = getCBendingsNofit(mdim, numFrames, tangents);
        }
        // tangents, G_n, best no longer needed — freed at scope exit
  
        // Relative sizes
        if cfgFitplane {
          var avgDistRelFrame: [1..numFrames] real;
          var maxDistRelFrame: [1..numFrames] real;
          for i in 1..numFrames {
            avgDistRelFrame[i] = 100.0 * (avgDistFrame[i] / height[i]);
            maxDistRelFrame[i] = 100.0 * (maxDistFrame[i] / height[i]);
          }
  
          // Averages
          writeln("Calculating averages and standard deviations");
          var cAvstdBends: [1..2, 1..mdim, 1..mdim-1] real;
          forall (i, l) in {1..mdim, 1..mdim-1} {
            var res = averageStd(cBends[1..numFrames, i, l]);
            cAvstdBends[1, i, l] = res[1];
            cAvstdBends[2, i, l] = res[2];
          }
  
          var avstdWidth = averageStd(width);
          var avstdHeight = averageStd(height);
          var avstdAratio = averageStd(aspectRatio);
          var avgDist = averageStd(avgDistFrame);
          var maxDist = averageStd(maxDistFrame);
          var avgDistRel = averageStd(avgDistRelFrame);
          var maxDistRel = averageStd(maxDistRelFrame);
  
          // Write
          writeln("Writing output modules");
          writeCAvParms(actualNbp, numFrames, finalSeq,
                        cAvstdBends, avstdWidth, avstdHeight, avstdAratio,
                        cfgTLength, avgDist, maxDist,
                        avgDistRel, maxDistRel);
        } else {
          // No fitplane: only bendings
          writeln("Calculating averages and standard deviations");
          var cAvstdBends: [1..2, 1..mdim, 1..mdim-1] real;
          forall (i, l) in {1..mdim, 1..mdim-1} {
            var res = averageStd(cBends[1..numFrames, i, l]);
            cAvstdBends[1, i, l] = res[1];
            cAvstdBends[2, i, l] = res[2];
          }
  
          writeln("Writing output modules");
          writeCAvBends(
            actualNbp, numFrames, finalSeq, cAvstdBends, cfgTLength);
        }
      } else {
        // bends: [1..F, 1..mldim] — O(F*N^2/2), manageable
        var bends: [1..numFrames, 1..mldim] real;
        if cfgFitplane {
          bends = getBendings(nldim, numFrames, G_n, best, tangents, mldim);
        } else {
          bends = getBendingsNofit(nldim, numFrames, tangents, mldim);
        }
        // tangents, G_n, best no longer needed — freed at scope exit
  
        // Relative sizes
        if cfgFitplane {
          var avgDistRelFrame: [1..numFrames] real;
          var maxDistRelFrame: [1..numFrames] real;
          for i in 1..numFrames {
            avgDistRelFrame[i] = 100.0 * (avgDistFrame[i] / height[i]);
            maxDistRelFrame[i] = 100.0 * (maxDistFrame[i] / height[i]);
          }
  
          // Averages
          writeln("Calculating averages and standard deviations");
          var avstdBends: [1..2, 1..mldim] real;
          forall i in 1..mldim {
            var res = averageStd(bends[1..numFrames, i]);
            avstdBends[1, i] = res[1];
            avstdBends[2, i] = res[2];
          }
  
          var avstdWidth = averageStd(width);
          var avstdHeight = averageStd(height);
          var avstdAratio = averageStd(aspectRatio);
          var avgDist = averageStd(avgDistFrame);
          var maxDist = averageStd(maxDistFrame);
          var avgDistRel = averageStd(avgDistRelFrame);
          var maxDistRel = averageStd(maxDistRelFrame);
  
          // Write
          writeln("Writing output modules");
          writeAvParms(actualNbp, nldim, numFrames, finalSeq,
                       avstdBends, avstdWidth, avstdHeight, avstdAratio,
                       cfgTLength, avgDist, maxDist,
                       avgDistRel, maxDistRel);
        } else {
          // No fitplane: only bendings
          writeln("Calculating averages and standard deviations");
          var avstdBends: [1..2, 1..mldim] real;
          forall i in 1..mldim {
            var res = averageStd(bends[1..numFrames, i]);
            avstdBends[1, i] = res[1];
            avstdBends[2, i] = res[2];
          }
  
          writeln("Writing output modules");
          writeAvBends(actualNbp, nldim, numFrames,
                       finalSeq, avstdBends, cfgTLength);
        }
      }
      // All per-frame arrays (bends/c_bends, width, height, etc.)
      // are now out of scope and freed.
    }
  }
  }
