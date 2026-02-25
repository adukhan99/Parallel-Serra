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
          var (subovElaspMean, subovElaspStd) = 
            centralFragmentMeanStd(elasp, a, b, nbp, str);
          writeOveralls2d(subovElaspMean, subovElaspStd,
                          "elastic_plot.out");
        }
        when 3 { // Structural
          var (subovStrucpMean, subovStrucpStd) = 
            centralFragmentMeanStd(strucp, a, b, nbp, str);
          var (subovAvstrpMean, subovAvstrpStd) = 
            centralFragmentMeanStd(avstrp, a, b, nbp, str);
          writeOveralls3d2d(subovStrucpMean, subovStrucpStd, 
                            subovAvstrpMean, subovAvstrpStd, 
                            "structural_plot.out");
        }
      }
    }
  }

  // Sorting by midpoint (bubble sort as in Fortran)
  proc sortByMid(mid: [] real, data: [?D2] real) {
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

  proc sortByMid(mid: [] real, data1: [?D2] real, data2: [?D3] real) {
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
  proc writeExtracted2D(mid: [] real, data: [] real,
                      filename: string) {
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

  proc writeExtracted3D2D(mid: [] real, data1: [] real, data2: [] real,
                        filename: string) {
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

  proc writeExtractedBPP(data: [] real, filename: string) {
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

  proc writeOveralls2D(mean: [] real, std: [] real,
                     filename: string) {
    var f = open(filename, ioMode.cw);
    var writer = f.writer();
    for i in 1..mean.size {
      writer.writef("%5d", i);
      writer.writef("%10.3f%10.3f", mean[i], std[i]);
      writer.writeln();
    }
    writer.close();
    f.close();
  }

  proc writeOveralls3D2D(m1: [] real, s1: [] real, m2: [] real,
                       s2: [] real, filename: string) {
    var f = open(filename, ioMode.cw);
    var writer = f.writer();
    for i in 1..m1.size {
      writer.writef("%5d", i);
      writer.writef("%10.3f%10.3f", m1[i], s1[i]);
      writer.writef("%10.3f%10.3f", m2[i], s2[i]);
      writer.writeln();
    }
    writer.close();
    f.close();
  }
}
