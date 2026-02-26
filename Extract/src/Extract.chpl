/*
 * Copyright (C) 2026 Alexander Dukhan
 * Credit owed to Victor Velasco and Agnes Noy whose fortran software
 * this is a direct re-write of.
 */
module Extract {
  use Parms;
  use Functions;
  use CustomIO;
  use List;
  use IO;
  use Map;

  config var configFile: string = "";
  config var fileIn: string = "";
  config var typeExt: int = 0;
  config var sublength: int = -1;
  config var a: int = 0;
  config var b: int = 0;

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
    var cfgFileIn = fileIn;
    var cfgTypeExt = typeExt;
    var cfgSublength = sublength;
    var cfgA = a;
    var cfgB = b;

    if configFile != "" {
      var cfg = parseConfigFile(configFile);
      if fileIn == "" && cfg.contains("fileIn") then
        cfgFileIn = try! cfg["fileIn"];
      if typeExt == 0 && cfg.contains("typeExt") then
        cfgTypeExt = try! cfg["typeExt"]:int;
      if sublength == -1 && cfg.contains("sublength") then
        cfgSublength = try! cfg["sublength"]:int;
      if a == 0 && cfg.contains("a") then
        cfgA = try! cfg["a"]:int;
      if b == 0 && cfg.contains("b") then
        cfgB = try! cfg["b"]:int;
    }

    if cfgFileIn == "" {
      writeln(
        "Error: fileIn must be provided (e.g. --fileIn=...) or via " +
        "--configFile=<file>");
      return;
    }
    var s_len = if cfgSublength >= 0 then cfgSublength else 0;
    extractParameters(cfgFileIn, cfgTypeExt, s_len, cfgA, cfgB);
  }

  /*
   * Extracts data for a sublength or subfragment from a parameter file.
   */
  proc extractParameters(fileIn: string, typeExt: int, sublength: int,
                       a: int, b: int) {
    // Determine the type of file and read it
    var (typeParm, nbp, frames, strands, str, nBsp, seqI, seqII, 
         BSP, ovBsp, elasp, ovElasp, strucp, avstrp, ovStrucp, ovAvstrp, 
         BPP, ovBpp) = readExtractInputFile(fileIn);

    var N = 0;
    var mid: list(real);

    if typeExt == 0 {
      // Sublength extraction
      var sublSize = if str == 2 then nbp else nbp - sublength;
      N = sublSize;
      for i in 1..N {
        var m = (2*i + sublength):real / 2.0;
        if m > nbp then m -= nbp:real;
        mid.pushBack(m);
      }

      var midArr = mid.toArray();

      // Extract and sort by mid point
      select typeParm {
        when 2 { // Elastic
          var sublElasp = extractSublength2D(elasp, nbp, sublength, str);
          sortByMid(midArr, sublElasp);
          writeExtracted2D(midArr, sublElasp,
                           "elastic_" + (sublength+1):string + "mer.out");
        }
        when 3 { // Structural
          var sublStrucp = extractSublength3D(strucp, nbp, sublength, str);
          var sublAvstrp = extractSublength2D(avstrp, nbp, sublength, str);
          sortByMid(midArr, sublStrucp, sublAvstrp);
          writeExtracted3D2D(midArr, sublStrucp, sublAvstrp,
                             "structural_" + (sublength+1):string + "mer.out");
        }
        when 4 { // BPP
          writeExtractedBPP(BPP, "BPP_plot.out");
        }
      }
    } else {
      // Subfragment overall extraction
      select typeParm {
        when 2 { // Elastic
                var overallSize = if a == 0 && b == 0 then nBsp 
                                  else if str == 2 then (if a < b then b-a else nbp-a+b)
                                  else b-a;          var subovElaspMean: [1..13, 1..overallSize] real;
          var subovElaspStd: [1..13, 1..overallSize] real;
          forall i in 1..13 {
            var (m, s) = centralFragmentMeanStd(elasp[i, ..], a, b, nbp, str);
            subovElaspMean[i, ..] = m;
            subovElaspStd[i, ..] = s;
          }
          writeOveralls2D(subovElaspMean, subovElaspStd,
                          "elastic_plot.out");
        }
        when 3 { // Structural
                var overallSize = if a == 0 && b == 0 then nBsp 
                                  else if str == 2 then (if a < b then b-a else nbp-a+b)
                                  else b-a;          var subovStrucpMean: [1..2, 1..11, 1..overallSize] real;
          var subovStrucpStd: [1..2, 1..11, 1..overallSize] real;
          forall (i, j) in {1..2, 1..11} {
            var (m, s) = centralFragmentMeanStd(strucp[i, j, ..], a, b, nbp, str);
            subovStrucpMean[i, j, ..] = m;
            subovStrucpStd[i, j, ..] = s;
          }
          var subovAvstrpMean: [1..3, 1..overallSize] real;
          var subovAvstrpStd: [1..3, 1..overallSize] real;
          forall i in 1..3 {
            var (m, s) = centralFragmentMeanStd(avstrp[i, ..], a, b, nbp, str);
            subovAvstrpMean[i, ..] = m;
            subovAvstrpStd[i, ..] = s;
          }
          writeOveralls3D2D(subovStrucpMean, subovStrucpStd, 
                             subovAvstrpMean, subovAvstrpStd, 
                             "structural_plot.out");
        }
      }
    }
  }

  // Sorting by midpoint (bubble sort as in Fortran)
  proc sortByMid(ref mid: [] real, ref data: [?D2] real) {
    var n = mid.size;
    for j in 1..n-1 {
      for i in 1..n-j {
        if mid[i] > mid[i+1] {
          var tmpM = mid[i];
          mid[i] = mid[i+1];
          mid[i+1] = tmpM;
          
          var tmpD: [D2.dim(0)] real = data[.., i];
          data[.., i] = data[.., i+1];
          data[.., i+1] = tmpD;
        }
      }
    }
  }

  proc sortByMid(ref mid: [] real, ref data1: [?D2] real, ref data2: [?D3] real) {
    var n = mid.size;
    for j in 1..n-1 {
      for i in 1..n-j {
        if mid[i] > mid[i+1] {
          var tmpM = mid[i];
          mid[i] = mid[i+1];
          mid[i+1] = tmpM;
          
          var tmpD1: [D2.dim(0), D2.dim(1)] real = data1[.., .., i];
          data1[.., .., i] = data1[.., .., i+1];
          data1[.., .., i+1] = tmpD1;

          var tmpD2: [D3.dim(0)] real = data2[.., i];
          data2[.., i] = data2[.., i+1];
          data2[.., i+1] = tmpD2;
        }
      }
    }
  }

  // Helper functions for writing extracted data
  proc writeExtracted2D(ref mid: [] real, ref data: [] real,
                      filename: string) {
    try! {
      var f = open(filename, ioMode.cw);
      var writer = f.writer();
      for i in 1..mid.size {
        writer.writef("%10.1f", mid[i]);
        for j in 1..data.dim(0).size {
          writer.writef("%20.3f", data[j, i]);
        }
        writer.writeln();
      }
      writer.close();
      f.close();
    }
  }

  proc writeExtracted3D2D(ref mid: [] real, ref data1: [] real, ref data2: [] real,
                        filename: string) {
    try! {
      var f = open(filename, ioMode.cw);
      var writer = f.writer();
      for i in 1..mid.size {
        writer.writef("%10.1f", mid[i]);
        for j in 1..data1.dim(1).size {
          for k in 1..data1.dim(0).size {
            writer.writef("%10.3f", data1[k, j, i]);
          }
        }
        for j in 1..data2.dim(0).size {
          writer.writef("%10.3f", data2[j, i]);
        }
        writer.writeln();
      }
      writer.close();
      f.close();
    }
  }

  proc writeExtractedBPP(ref data: [] real, filename: string) {
    try! {
      var f = open(filename, ioMode.cw);
      var writer = f.writer();
      var nbp = data.dim(2).size;
      for i in 1..nbp {
        writer.writef("%5d", i);
        for j in 1..data.dim(1).size {
          for k in 1..data.dim(0).size {
            writer.writef("%10.3f", data[k, j, i]);
          }
        }
        writer.writeln();
      }
      writer.close();
      f.close();
    }
  }

  proc writeOveralls2D(ref mean: [] real, ref std: [] real,
                     filename: string) {
    try! {
      var f = open(filename, ioMode.cw);
      var writer = f.writer();
      // Assuming the last dimension is the step dimension
      var nSteps = mean.dim(mean.rank-1).size;
      for i in 1..nSteps {
        writer.writef("%5d", i+1);
        if mean.rank == 1 {
          writer.writef("%10.3f%10.3f", mean[i], std[i]);
        } else if mean.rank == 2 {
          for j in 1..mean.dim(0).size {
            writer.writef("%10.3f%10.3f", mean[j, i], std[j, i]);
          }
        }
        writer.writeln();
      }
      writer.close();
      f.close();
    }
  }

  proc writeOveralls3D2D(ref m1: [] real, ref s1: [] real, ref m2: [] real,
                       ref s2: [] real, filename: string) {
    try! {
      var f = open(filename, ioMode.cw);
      var writer = f.writer();
      var nSteps = m1.dim(m1.rank-1).size;
      for i in 1..nSteps {
        writer.writef("%5d", i+1);
        if m1.rank == 3 {
          for j in 1..m1.dim(0).size {
            for k in 1..m1.dim(1).size {
              writer.writef("%10.3f%10.3f", m1[j, k, i], s1[j, k, i]);
            }
          }
        }
        if m2.rank == 2 {
          for j in 1..m2.dim(0).size {
            writer.writef("%10.3f%10.3f", m2[j, i], s2[j, i]);
          }
        }
        writer.writeln();
      }
      writer.close();
      f.close();
    }
  }
}
