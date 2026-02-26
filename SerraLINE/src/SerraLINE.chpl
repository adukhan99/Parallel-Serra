/* This is a re-write of the SerraLINE package developed
   by Victor Manuel Velasco Berrelleza and Agnes Noy. All Chapel
   code is written by Alexander Dukhan. */

module SerraLINE {
  use Parms;
  use CustomIO;
  use Functions;
  use Math;
  use LinearAlgebra;
  use IO;
  use Map;

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

  proc writeAvParms(nbp: int, ndim: int, frames: int, ref seq: [] string,
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
      writer.writef(F_PARM_5, "SNAPSHOTS ANALYSED", frames);
      writer.writeln();
      writer.writef(F_PARM_5, "TANGENT LENGTH:", t_length);
      writer.writeln();
      writer.writeln(
        "First column averages, second column standard deviations");
      writer.writef(
        F_PARM_6, "WIDTH (Angstroms): ", avstd_width[1], avstd_width[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6, "HEIGHT (Angstroms): ", avstd_height[1], avstd_height[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6, "ASPECT RATIO:", avstd_aratio[1], avstd_aratio[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6,
        "AVERAGE OF DISTANCES TO PLANE (Angstroms):",
        avg_dist[1], avg_dist[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6,
        "AVERAGE OF MAXIMUM DISTANCES TO PLANE (Angstroms):",
        max_dist[1], max_dist[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6,
        "RELATIVE AVERAGE OF DISTANCES TO PLANE (%):",
        avg_dist_rel[1], avg_dist_rel[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6,
        "RELATIVE AVERAGE OF MAXIMUM DISTANCES TO PLANE (%):",
        max_dist_rel[1], max_dist_rel[2]);
      writer.writeln("\n");
      var l = 0;
      for j in 1..ndim-1 {
        writer.writef(F_PARM_3, j+1, "mer");
        writer.writeln();
        writer.writef(F_PARM_1, "base-step", "Bending angle");
        writer.writeln();
        writer.writeln("-" * 38);
        var endLimit = ndim - j;
        if endLimit >= 1 {
          for i in 1..endLimit {
            l += 1;
            writer.writef(
              F_PARM_2, i, "-", i+j, " ", " ",
              seq[i], seq[i+j], " ", " ",
              avstd_bends[1, l], avstd_bends[2, l]);
            writer.writeln();
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
                     t_length: int, ref avg_dist: [] real, ref max_dist: [] real,
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
      writer.writef(F_PARM_5, "SNAPSHOTS ANALYSED", frames);
      writer.writeln();
      writer.writef(F_PARM_5, "TANGENT LENGTH:", t_length);
      writer.writeln();
      writer.writeln(
        "First column averages, second column standard deviations");
      writer.writef(
        F_PARM_6, "WIDTH (Angstroms): ", avstd_width[1], avstd_width[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6, "HEIGHT (Angstroms): ", avstd_height[1], avstd_height[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6, "ASPECT RATIO:", avstd_aratio[1], avstd_aratio[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6,
        "AVERAGE OF DISTANCES TO PLANE (Angstroms):",
        avg_dist[1], avg_dist[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6,
        "AVERAGE OF MAXIMUM DISTANCES TO PLANE (Angstroms):",
        max_dist[1], max_dist[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6,
        "RELATIVE AVERAGE OF DISTANCES TO PLANE (%):",
        avg_dist_rel[1], avg_dist_rel[2]);
      writer.writeln();
      writer.writef(
        F_PARM_6,
        "RELATIVE AVERAGE OF MAXIMUM DISTANCES TO PLANE (%):",
        max_dist_rel[1], max_dist_rel[2]);
      writer.writeln("\n");
      for l in 1..nbp-1 {
        writer.writef(F_PARM_3, l+1, "mer");
        writer.writeln();
        writer.writef(F_PARM_1, "base-step", "Bending angle");
        writer.writeln();
        writer.writeln("-" * 38);
        for i in 1..nbp {
          var s = i + l;
          if s > nbp then s = s - nbp;
          writer.writef(
            F_PARM_2, i, "-", s, " ", " ",
            c_seq[i], c_seq[s], " ", " ",
            avstd_bends[1, i, l], avstd_bends[2, i, l]);
          writer.writeln();
        }
        writer.writeln();
      }
      writer.close();
    }
  }

  proc writeAvBends(nbp: int, ndim: int, frames: int, ref seq: [] string,
                    ref avstd_bends: [] real, t_length: int) {
    try! {
      var f = open("SerraLINE.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln("PARAMETERS\n");
      writer.writeln("OPENED STRUCTURE");
      writer.writeln("METHOD: WITHOUT PROJECTION");
      writer.writeln("BASE PAIRS ", nbp);
      writer.write("SEQUENCE\n");
      for s in seq do writer.write(s);
      writer.writeln();
      writer.writef(F_PARM_5, "SNAPSHOTS ANALYSED", frames);
      writer.writeln();
      writer.writef(F_PARM_5, "TANGENT LENGTH:", t_length);
      writer.writeln();
      writer.writeln(
        "First column averages, second column standard deviations\n");
      var l = 0;
      for j in 1..ndim-1 {
        writer.writef(F_AVB_3, j+1, "mer");
        writer.writeln();
        writer.writef(F_AVB_1, "base-step", "Bending angle");
        writer.writeln();
        writer.writeln("-" * 38);
        var endLimit = ndim - j;
        if endLimit >= 1 {
          for i in 1..endLimit {
            l += 1;
            writer.writef(
              F_AVB_2, i, "-", i+j, " ", " ",
              seq[i], seq[i+j], " ", " ",
              avstd_bends[1, l], avstd_bends[2, l]);
            writer.writeln();
          }
        }
        writer.writeln();
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
      writer.writef(F_PARM_5, "SNAPSHOTS ANALYSED", frames);
      writer.writeln();
      writer.writef(F_PARM_5, "TANGENT LENGTH:", t_length);
      writer.writeln();
      writer.writeln(
        "First column averages, second column standard deviations\n");
      for l in 1..nbp-1 {
        writer.writef(F_AVB_3, l+1, "mer");
        writer.writeln();
        writer.writef(F_AVB_1, "base-step", "Bending angle");
        writer.writeln();
        writer.writeln("-" * 38);
        for i in 1..nbp {
          var s = i + l;
          if s > nbp then s = s - nbp;
          writer.writef(
            F_AVB_2, i, "-", s, " ", " ",
            c_seq[i], c_seq[s], " ", " ",
            avstd_bends[1, i, l], avstd_bends[2, i, l]);
          writer.writeln();
        }
        writer.writeln();
      }
      writer.close();
    }
  }

  proc main() {
    try! {
    // Load config file values as defaults (CLI flags override them).
    var cfgTop = top;
    var cfgTraj = traj;
    var cfgIsCircular = isCircular;
    var cfgStrandsType = strandsType;
    var cfgNbpArg = nbpArg;
    var cfgFitplane = fitplane;
    var cfgBpFitting = bpFitting;
    var cfgTLength = tLength;
    var cfgPrintProj = printProj;
    var cfgXyzFormat = xyzFormat;

    if configFile != "" {
      var cfg = parseConfigFile(configFile);
      if top == "" && cfg.contains("top") then
        cfgTop = try! cfg["top"];
      if traj == "" && cfg.contains("traj") then
        cfgTraj = try! cfg["traj"];
      if !isCircular && cfg.contains("isCircular") then
        cfgIsCircular = (try! cfg["isCircular"] == "true");
      if strandsType == 1 && cfg.contains("strandsType") then
        cfgStrandsType = try! cfg["strandsType"]:int;
      if nbpArg == 0 && cfg.contains("nbpArg") then
        cfgNbpArg = try! cfg["nbpArg"]:int;
      if fitplane && cfg.contains("fitplane") then
        cfgFitplane = (try! cfg["fitplane"] == "true");
      if bpFitting == "0" && cfg.contains("bpFitting") then
        cfgBpFitting = try! cfg["bpFitting"];
      if tLength == 1 && cfg.contains("tLength") then
        cfgTLength = try! cfg["tLength"]:int;
      if !printProj && cfg.contains("printProj") then
        cfgPrintProj = (try! cfg["printProj"] == "true");
      if !xyzFormat && cfg.contains("xyzFormat") then
        cfgXyzFormat = (try! cfg["xyzFormat"] == "true");
    }

    if cfgTraj == "" {
      writeln(
        "Trajectory file must be specified using --traj=<file>" +
        " or via --configFile=<file>");
      return;
    }

    var actualNbp = cfgNbpArg;
    var seqDom = {1..0};
    var finalSeq: [seqDom] string;
    var coordDom = {1..3, 1..0, 1..0};
    var filteredCoords: [coordDom] real;
    var numFrames = 0;

    if cfgStrandsType != 0 && cfgTop != "" {
      writeln("Reading topology file");
      var (origNbp, nAtoms, box, origSeq, ringIndices) =
        topologyAmber(cfgTop, cfgStrandsType);

      // SerraLINE uses raw atom positions (one per base-pair), not ring atoms.
      // coordinatesOther reads exactly nAtoms = origNbp atoms per frame.
      writeln("Reading trajectory file");
      var (coords, frames) = coordinatesOther(cfgTraj, nAtoms);
      actualNbp = origNbp;
      numFrames = frames;
      coordDom = coords.domain;
      filteredCoords = coords;
      seqDom = {1..actualNbp};
      finalSeq = origSeq[1..actualNbp];
    } else {
      writeln("Reading trajectory file (no topology)");
      if actualNbp <= 0 then
        halt("Must specify --nbpArg since no topology is provided");
      var (coords, fn) = coordinatesOther(cfgTraj, actualNbp);
      numFrames = fn;
      coordDom = coords.domain;
      filteredCoords = coords;
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


    var G_n: [1..3, 1..3, 1..numFrames] real;
    var best: [1..numFrames] int;
    var avg_dist_frame: [1..numFrames] real;
    var max_dist_frame: [1..numFrames] real;

    // Parse bpFitting: "0" or empty means use all base-pairs;
    // otherwise it is a comma-separated list of 1-based indices.
    var bpFittingDom = {1..0};
    var bp_fitting_arr: [bpFittingDom] int;
    var n_fitting = 0;
    if cfgBpFitting != "0" && cfgBpFitting != "" {
      use List;
      var idxList: list(int);
      for tok in cfgBpFitting.split(",") {
        var idx = tok.strip():int;
        if idx >= 1 && idx <= actualNbp then idxList.pushBack(idx);
      }
      if idxList.size > 0 {
        var arr = idxList.toArray();
        bpFittingDom = arr.domain;
        bp_fitting_arr = arr;
        n_fitting = bp_fitting_arr.size;
      }
    }

    if cfgFitplane {
      writeln("Fitting planes");
      var res = projectCoordinatesGPlane(
        filteredCoords, n_fitting, bp_fitting_arr);
      G_n = res(0);
      best = res(1);
      avg_dist_frame = res(2);
      max_dist_frame = res(3);
      if cfgPrintProj then
        writeProjectedCoords(
          filteredCoords, G_n, best, finalSeq, cfgXyzFormat);
    }

    var width: [1..numFrames] real;
    var height: [1..numFrames] real;
    var aspect_ratio: [1..numFrames] real;

    if cfgFitplane {
      writeln("Calculating widths, heights and aspect ratios");
      if cfgIsCircular {
        (width, height, aspect_ratio) =
          getCWidthHeightAratio(mdim, numFrames, G_n, best, filteredCoords);
      } else {
        (width, height, aspect_ratio) =
          getWidthHeightAratio(ndim, mdim, numFrames,
                              G_n, best, filteredCoords);
      }
    }

    writeln("Calculating tangent vectors");
    var tangents_n = if cfgIsCircular then actualNbp else nldim;
    var tangents: [1..3, 1..tangents_n, 1..numFrames] real;
    tangents = getTangentVectors(filteredCoords, cfgIsCircular, cfgTLength);

    writeln("Calculating bending angles");
    var c_bends: [1..numFrames, 1..mdim, 1..mdim-1] real;
    var bends: [1..numFrames, 1..mldim] real;

    if cfgIsCircular {
      if cfgFitplane {
        c_bends = getCBendings(mdim, numFrames, G_n, best, tangents);
      } else {
        c_bends = getCBendingsNofit(mdim, numFrames, tangents);
      }
    } else {
      if cfgFitplane {
        bends = getBendings(nldim, numFrames, G_n, best, tangents, mldim);
      } else {
        bends = getBendingsNofit(nldim, numFrames, tangents, mldim);
      }
    }

    var avg_dist_rel_frame: [1..numFrames] real;
    var max_dist_rel_frame: [1..numFrames] real;

    if cfgFitplane {
      forall i in 1..numFrames {
        avg_dist_rel_frame[i] = 100.0 * (avg_dist_frame[i] / height[i]);
        max_dist_rel_frame[i] = 100.0 * (max_dist_frame[i] / height[i]);
      }
    }

    writeln("Calculating averages and standard deviations");
    var c_avstd_bends: [1..2, 1..mdim, 1..mdim-1] real;
    var avstd_bends: [1..2, 1..mldim] real;
    var avstd_width, avstd_height, avstd_aratio,
        avg_dist, max_dist, avg_dist_rel, max_dist_rel: [1..2] real;

    if cfgIsCircular {
      forall (i, l) in {1..mdim, 1..mdim-1} {
        var res = averageStd(c_bends[1..numFrames, i, l]);
        c_avstd_bends[1, i, l] = res[1];
        c_avstd_bends[2, i, l] = res[2];
      }
    } else {
      forall i in 1..mldim {
        var res = averageStd(bends[1..numFrames, i]);
        avstd_bends[1, i] = res[1];
        avstd_bends[2, i] = res[2];
      }
    }

    if cfgFitplane {
      var wRes = averageStd(width);
      avstd_width[1] = wRes[1];
      avstd_width[2] = wRes[2];
      var hRes = averageStd(height);
      avstd_height[1] = hRes[1];
      avstd_height[2] = hRes[2];
      var arRes = averageStd(aspect_ratio);
      avstd_aratio[1] = arRes[1];
      avstd_aratio[2] = arRes[2];

      var adRes = averageStd(avg_dist_frame);
      avg_dist[1] = adRes[1];
      avg_dist[2] = adRes[2];
      var mdRes = averageStd(max_dist_frame);
      max_dist[1] = mdRes[1];
      max_dist[2] = mdRes[2];
      var ardlRes = averageStd(avg_dist_rel_frame);
      avg_dist_rel[1] = ardlRes[1];
      avg_dist_rel[2] = ardlRes[2];
      var mrdlRes = averageStd(max_dist_rel_frame);
      max_dist_rel[1] = mrdlRes[1];
      max_dist_rel[2] = mrdlRes[2];
    }

    writeln("Writing output modules");
    if cfgIsCircular {
      if cfgFitplane {
        writeCAvParms(actualNbp, numFrames, finalSeq,
                      c_avstd_bends, avstd_width, avstd_height, avstd_aratio,
                      cfgTLength, avg_dist, max_dist,
                      avg_dist_rel, max_dist_rel);
      } else {
        writeCAvBends(
          actualNbp, numFrames, finalSeq, c_avstd_bends, cfgTLength);
      }
    } else {
      if cfgFitplane {
        writeAvParms(actualNbp, ndim, numFrames, finalSeq,
                     avstd_bends, avstd_width, avstd_height, avstd_aratio,
                     cfgTLength, avg_dist, max_dist,
                     avg_dist_rel, max_dist_rel);
      } else {
        writeAvBends(actualNbp, ndim, numFrames,
                     finalSeq, avstd_bends, cfgTLength);
      }
    }
    }
  }
}
