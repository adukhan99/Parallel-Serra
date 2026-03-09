/*
 * Copyright (C) 2026 Alexander Dukhan
 * Credit owed to Victor Velasco and Agnes Noy whose fortran software
 * this is a direct re-write of.
 */
module CustomIO {
  use Parms;
  use Functions;
  use IO;
  use List;
  use FileSystem;

  /* 
   * Reads topology file from AMBER topology format.
   * Extracts number of base-pairs, atoms, residues, and the sequence.
   */
  proc topologyAmber(topPath: string, strandsType: int) {
    try! {
      var nbp, nAtoms, nRes, box: int;
      
      var f = open(topPath, ioMode.r);
      var reader = f.reader();
      
      var line: string;
      var resNames, atomNames: list(string);
      var resPointers: list(int);

      while reader.readLine(line) {
        line = line.strip();
        if line.startsWith("%FLAG") then writeln("Processing flag: ", line);
        if line == "%FLAG POINTERS" {
          reader.readLine(line); // Skip %FORMAT
          nAtoms = reader.read(int); // 1: NATOM
          for 2..11 do reader.read(int); // skip to 11
          nRes = reader.read(int); // 12: NRES
          box = reader.read(int); // 13: IFBOX
          writeln("Parsed POINTERS: nAtoms=", nAtoms,
                  " nRes=", nRes, " box=", box);
        } else if line == "%FLAG ATOM_NAME" {
          reader.readLine(line); // Skip %FORMAT
          for 1..nAtoms {
            var s = "";
            for 1..4 {
              var c: string;
              reader.readf("%c", c);
              if c == "\n" {
                reader.readf("%c", c);
              }
              s += c;
            }
            atomNames.pushBack(s);
          }
        }
        else if line == "%FLAG RESIDUE_LABEL" {
          reader.readLine(line); // Skip %FORMAT
          writeln("Reading ", nRes, " residue labels");
          for 1..nRes {
            var s = "";
            for 1..4 {
              var c: string;
              reader.readf("%c", c);
              if c == "\n" {
                reader.readf("%c", c);
              }
              s += c;
            }
            resNames.pushBack(s);
            if resNames.size <= 10 then
              writeln("Read residue label: '", s, "'");
          }
        } else if line == "%FLAG RESIDUE_POINTER" {
          reader.readLine(line); // Skip %FORMAT
          for 1..nRes do resPointers.pushBack(reader.read(int));
        }
      }
      f.close();

      // Identify nucleotide residues and sequence
      var rawSeq: list(string);
      var rawResPointers: list(int);


      for (name, ptr) in zip(resNames, resPointers) {
        var found = false;
        var base = "";
        for n in A_l do if name.strip() == n.strip() {
          base = "A";
          found = true;
        }
        for n in G_l do if name.strip() == n.strip() {
          base = "G";
          found = true;
        }
        for n in C_l do if name.strip() == n.strip() {
          base = "C";
          found = true;
        }
        for n in T_l do if name.strip() == n.strip() {
          base = "T";
          found = true;
        }
        for n in U_l do if name.strip() == n.strip() {
          base = "U";
          found = true;
        }
        if found {
          rawSeq.pushBack(base);
          rawResPointers.pushBack(ptr);
        }
      }

      nbp = rawSeq.size;
      writeln("Found ", nbp, " residues matching DNA labels");
      if strandsType == 2 {
        if nbp % 2 != 0 then
          halt("Incomplete double stranded structure: ", nbp);
        nbp /= 2;
      }
      writeln("Final nbp = ", nbp);

      if nbp <= 4 then
        halt("DNA fragment must be larger than 4 bp (nbp=", nbp, ")");

      // Extract ring atoms indices
      var ringIndices: list(int);
      // For XYZ files from WrLINE, we only want the first nbp sequence elements
      // if we are treating the trajectory as having one point per bp step.
      var finalSeqList: list(string);
      if strandsType == 2 {
        for i in 0..nbp-1 do finalSeqList.pushBack(rawSeq[i]);
      } else {
        for i in 0..rawSeq.size-1 do finalSeqList.pushBack(rawSeq[i]);
      }

      for i in 0..rawSeq.size-1 {
        var s1 = rawResPointers[i];
        var s2 = if i < rawSeq.size-1 then rawResPointers[i+1] - 1 else nAtoms;
        
        var base = rawSeq[i];
        if base == "A" || base == "G" {
          // Purine: N9, C8, N7, C5, C6, N1, C2, N3, C4
          var targets = [N9, C8, N7, C5, C6, N1, C2, N3, C4];
          for target in targets {
            var foundIdx = 0;
            for j in s1..s2 {
              if atomNames[j-1].strip() == target.strip() {
                foundIdx = j;
                break;
              }
            }
            ringIndices.pushBack(foundIdx);
          }
        } else {
          // Pyrimidine: N1, C2, N3, C4, C5, C6
          var targets = [N1, C2, N3, C4, C5, C6];
          for target in targets {
            var foundIdx = 0;
            for j in s1..s2 {
              if atomNames[j-1].strip() == target.strip() {
                foundIdx = j;
                break;
              }
            }
            ringIndices.pushBack(foundIdx);
          }
        }
      }

      return (nbp, nAtoms, box, finalSeqList.toArray(), ringIndices.toArray());
    }
  }

  /*
   * Reads atom coordinates from AMBER trajectory format (.crd / .x).
   * Supports both open and closed structures.
   */
  proc coordinatesAmberCrd(trajPath: string, nAtoms: int, nbp: int,
                             boxPresent: int, strandsType: int,
                             isCircular: bool, ringIndices: [] int,
                             seq: [] string) {
    try! {
      var f = open(trajPath, ioMode.r);
      var reader = f.reader();
      
      var header: string;
      reader.readLine(header);

      // Determine frames
      var framesList: list([1..3, 1..nAtoms] real);
      
      try {
        while true {
          var frame: [1..3, 1..nAtoms] real;
          for i in 1..nAtoms {
            for j in 1..3 do frame[j, i] = reader.read(real);
          }
          if boxPresent != 0 {
            for 1..3 do reader.read(real); // Skip box dimensions
          }
          framesList.pushBack(frame);
        }
      } catch {
        // End of file
      }
      f.close();

      var numFrames = framesList.size;
      
      // Filter coordinates to ring atoms and handle open structure ends
      var startBp = 1,
          endBp = nbp,
          actualNbp = nbp;

      if !isCircular {
        actualNbp = nbp - 4;
        startBp = 3;
        endBp = nbp - 2;
      }

      var rAtomsPerStrand = 0;
      for i in startBp..endBp {
        if seq[i] == "A" || seq[i] == "G" then rAtomsPerStrand += 9;
        else rAtomsPerStrand += 6;
      }

      var totalRAtoms = if strandsType == 2 
                          then 2 * rAtomsPerStrand 
                          else rAtomsPerStrand;
      var filteredCoords: [1..3, 1..totalRAtoms, 1..numFrames] real;

      forall k in 1..numFrames {
        var l = 0;
        // Strand 1
        var ringPtr = 1;
        // Skip first 2 bp if linear
        if !isCircular {
          for i in 1..2 {
            if seq[i] == "A" || seq[i] == "G" then ringPtr += 9;
            else ringPtr += 6;
          }
        }
        
        for i in startBp..endBp {
          var atomsInBase = if seq[i] == "A" || seq[i] == "G" then 9 else 6;
          for 1..atomsInBase {
            l += 1;
            var atomIdx = ringIndices[ringPtr];
            filteredCoords[1..3, l, k] = framesList[k][1..3, atomIdx];
            ringPtr += 1;
          }
        }

        // Strand 2
        if strandsType == 2 {
          // Re-calculate ringPtr for second strand
          var ringPtr2 = 1;
          for i in 1..nbp {
             if seq[i] == "A" || seq[i] == "G" then ringPtr2 += 9;
             else ringPtr2 += 6;
          }
          // Now at start of strand 2
          if !isCircular {
             for i in 1..2 {
               if seq[nbp+i] == "A" || seq[nbp+i] == "G" then ringPtr2 += 9;
               else ringPtr2 += 6;
             }
          }
          for i in startBp..endBp {
            var idxInSeq = nbp + i;
            var isPurine = seq[idxInSeq] == "A" || seq[idxInSeq] == "G";
            var atomsInBase = if isPurine then 9 else 6;
            for 1..atomsInBase {
              l += 1;
              var atomIdx = ringIndices[ringPtr2];
              filteredCoords[1..3, l, k] = framesList[k][1..3, atomIdx];
              ringPtr2 += 1;
            }
          }
        }
      }

      var finalSeq: list(string);
      for i in startBp..endBp do finalSeq.pushBack(seq[i]);
      if strandsType == 2 {
        for i in startBp..endBp do finalSeq.pushBack(seq[nbp + i]);
      }

      return (filteredCoords, numFrames, actualNbp, finalSeq.toArray());
    }
  }

  /* 
   * Reads 3-column coordinates (.3col) or XYZ format (.xyz) 
   */
  proc coordinatesOther(path: string, nAtoms: int) {
    try! {
      // First pass: count frames
      var f1 = open(path, ioMode.r);
      var r1 = f1.reader();
      var isXyz = path.endsWith(".xyz");
      var numFrames = 0;
      try {
        while true {
          if isXyz {
            r1.read(int);
            r1.readLine();
          }
          for 1..nAtoms {
            if isXyz then r1.read(string);
            for 1..3 do r1.read(real);
          }
          numFrames += 1;
        }
      } catch { }
      f1.close();

      // Second pass: read directly into array
      var coords: [1..3, 1..nAtoms, 1..numFrames] real;
      var f2 = open(path, ioMode.r);
      var r2 = f2.reader();
      for k in 1..numFrames {
        if isXyz {
          r2.read(int);
          r2.readLine();
        }
        for i in 1..nAtoms {
          if isXyz then r2.read(string);
          for j in 1..3 do coords[j, i, k] = r2.read(real);
        }
      }
      f2.close();

      return (coords, numFrames);
    }
  }

  /* Read extract input file and identify its type */
  proc readExtractInputFile(fileIn: string) {
    try! {
      var f = open(fileIn, ioMode.r);
      var reader = f.reader();
      var cdum: string;
      reader.readLine(cdum);
      f.close();

      if cdum.find("BASE-STEP PARAMETERS") != -1 {
        return readBSP(fileIn);
      }
      if cdum.find("ELASTIC PARAMETERS") != -1 {
        return readElasticParms(fileIn);
      }
      if cdum.find("STRUCTURAL PARAMETERS") != -1 {
        return readStructuralParms(fileIn);
      }
      if cdum.find("BASE-PAIR PARAMETERS") != -1 {
        return readBPP(fileIn);
      }
      
      halt("Couldn't identify the type of data file: " + cdum);
    }
  }

  /* Internal helper to read the header/metadata part of parameter files */
  private proc readMetadataHelper(reader) {
    try! {
      var cdum: string;
      var str, strands, nbp, frames: int;

      reader.read(cdum); // CLOSED or LINEAR
      str = if cdum == "CLOSED" then 2 else 1;
      reader.read(cdum); // STRUCTURE
      reader.read(cdum); // DOUBLE-STRANDED or SINGLE-STRANDED
      strands = if cdum == "DOUBLE-STRANDED" then 2 else 1;
      reader.read(cdum); // STRUCTURE
      reader.read(cdum); // ANALYSED
      reader.read(cdum); // BASE-PAIRS:
      reader.read(nbp);
      reader.read(cdum); // FRAMES:
      reader.read(frames);
      reader.read(cdum); // SEQUENCE:
      reader.read(cdum); // STRAND
      reader.read(cdum); // 1:
      
      var seqI: [1..nbp] string;
      for i in 1..nbp do reader.read(seqI[i]);
      
      var seqII: [1..nbp] string;
      reader.read(cdum); // STRAND
      reader.read(cdum); // 2:
      if strands == 2 {
        for i in 1..nbp do reader.read(seqII[i]);
      }
      
      return (nbp, frames, strands, str, seqI, seqII);
    }
  }

  proc readBSP(fileIn: string) {
    try! {
      var f = open(fileIn, ioMode.r);
      var reader = f.reader();
      var cdum: string;
      reader.readLine(cdum); // BASE-STEP PARAMETERS
      
      var (nbp, frames, strands, str, seqI, seqII) = readMetadataHelper(reader);
      var nBsp = if str == 2 then nbp else nbp - 1;
      
      var BSP: [1..2, 1..7, 1..nBsp] real,
          ovBsp: [1..2, 1..7, 1..1] real;
      
      reader.read(cdum); // First column averages...
      reader.readLine(cdum); 
      reader.readLine(cdum); // base-step...
      reader.readLine(cdum); // ----

      for l in 1..nBsp {
        var idum: int;
        var sdum: string;
        reader.read(idum); // i
        reader.read(sdum); // -
        reader.read(idum); // w
        for 1..6 do reader.read(sdum); // base types and separators
        for i in 1..2 do for j in 1..7 do reader.read(BSP[i, j, l]);
      }
      reader.readLine(cdum); // ----
      reader.read(cdum); // AVG-STD=
      for i in 1..2 do for j in 1..7 do reader.read(ovBsp[i, j, 1]);
      
      f.close();
      return (1:int, nbp, frames, strands, str, nBsp, seqI, seqII, BSP, ovBsp, 
              matrix(0, 0), array3D(0, 0, 0), // elasp
              // strucp, avstrp, ovStrucp, ovAvstrp
              array3D(0, 0, 0), matrix(0, 0), 
              array3D(0, 0, 0), array3D(0, 0, 0),
              array3D(0, 0, 0), matrix(0, 0)); // BPP
    }
  }

  // Define placeholders for empty arrays to keep tuple size consistent
  private proc matrix(d1: int, d2: int) {
    var A: [1..d1, 1..d2] real;
    return A;
  }
  private proc array3D(d1: int, d2: int, d3: int) {
    var A: [1..d1, 1..d2, 1..d3] real;
    return A;
  }

  proc readElasticParms(fileIn: string) {
    try! {
      var f = open(fileIn, ioMode.r);
      var reader = f.reader();
      var cdum: string;
      reader.readLine(cdum); // ELASTIC PARAMETERS
      
      var (nbp, frames, strands, str, seqI, seqII) = readMetadataHelper(reader);
      var nBsp = if str == 2 then nbp*(nbp-1) else nbp*(nbp-1)/2;
      
      var elasp: [1..13, 1..nBsp] real,
          ovElasp: [1..2, 1..13, 1..nbp-1] real;
      
      for j in 1..nbp-1 {
        reader.read(cdum); // bp j
        reader.read(cdum); 
        reader.read(cdum); // base-step ...
        reader.read(cdum); // ----
        
        var nSteps = if str == 2 then nbp else nbp - j;
        for i in 1..nSteps {
          var idum: int;
          var sdum: string;
          reader.read(idum); // i
          reader.read(sdum); // -
          reader.read(idum); // w
          for 1..6 do reader.read(sdum); 
          var idx = if str == 2 then (j-1)*nbp + i 
            else (j-1)*nbp - (j-1)*j/2 + i;
          for k in 1..13 do reader.read(elasp[k, idx]);
        }
        reader.read(cdum); // ----
        reader.read(cdum); // AVG-STD=
        for k in 1..2 do for l in 1..13 do reader.read(ovElasp[k, l, j]);
      }
      
      f.close();
      return (2:int, nbp, frames, strands, str, nBsp, seqI, seqII, 
              array3D(0, 0, 0), array3D(0, 0, 0), // BSP
              elasp, ovElasp, 
              // strucp, avstrp, ovStrucp, ovAvstrp
              array3D(0, 0, 0), matrix(0, 0),
              array3D(0, 0, 0), array3D(0, 0, 0),
              array3D(0, 0, 0), matrix(0, 0)); // BPP
    }
  }

  proc readStructuralParms(fileIn: string) {
    try! {
      var f = open(fileIn, ioMode.r);
      var reader = f.reader();
      var cdum: string;
      reader.readLine(cdum); // STRUCTURAL PARAMETERS
      
      var (nbp, frames, strands, str, seqI, seqII) = readMetadataHelper(reader);
      var nBsp = if str == 2 then nbp*(nbp-1) else nbp*(nbp-1)/2;
      
      var strucp: [1..2, 1..11, 1..nBsp] real,
          avstrp: [1..3, 1..nBsp] real,
          ovStrucp: [1..2, 1..11, 1..nbp-1] real,
          ovAvstrp: [1..2, 1..3, 1..nbp-1] real;

      for j in 1..nbp-1 {
        reader.read(cdum); // bp j
        reader.read(cdum);
        reader.read(cdum); // base-step...
        reader.read(cdum); // ----
        
        var nSteps = if str == 2 then nbp else nbp - j;
        for i in 1..nSteps {
          var idum: int;
          var sdum: string;
          reader.read(idum); 
          reader.read(sdum);
          reader.read(idum);
          for 1..6 do reader.read(sdum);
          var idx = if str == 2 then (j-1)*nbp + i 
            else (j-1)*nbp - (j-1)*j/2 + i;
          for k in 1..2 do for l in 1..11 do reader.read(strucp[k, l, idx]);
          for k in 1..3 do reader.read(avstrp[k, idx]);
        }
        reader.read(cdum); // ----
        reader.read(cdum); // AVG-STD=
        for k in 1..2 do for l in 1..11 do reader.read(ovStrucp[k, l, j]);
        for k in 1..2 do for l in 1..3 do reader.read(ovAvstrp[k, l, j]);
      }
      
      f.close();
      return (3:int, nbp, frames, strands, str, nBsp, seqI, seqII, 
              array3D(0, 0, 0), array3D(0, 0, 0), // BSP
              matrix(0, 0), array3D(0, 0, 0), // elasp
              strucp, avstrp, ovStrucp, ovAvstrp,
              array3D(0, 0, 0), matrix(0, 0)); // BPP
    }
  }

  proc readBPP(fileIn: string) {
    try! {
      var f = open(fileIn, ioMode.r);
      var reader = f.reader();
      var cdum: string;
      reader.readLine(cdum); // BASE-PAIR PARAMETERS
      
      var (nbp, frames, strands, str, seqI, seqII) = readMetadataHelper(reader);
      var nBsp = nbp;
      
      var BPP: [1..2, 1..6, 1..nBsp] real,
          ovBpp: [1..2, 1..6] real;

      reader.read(cdum); // base-pair ...
      reader.readLine(cdum);
      reader.readLine(cdum); // ----
      
      for i in 1..nbp {
        var idum: int;
        var sdum: string;
        reader.read(idum); // i
        for 1..4 do reader.read(sdum); // base types and separators
        for k in 1..2 do for l in 1..6 do reader.read(BPP[k, l, i]);
      }
      reader.readLine(cdum); // ----
      reader.read(cdum); // AVG-STD=
      for k in 1..2 do for l in 1..6 do reader.read(ovBpp[k, l]);
      
      f.close();
      return (4:int, nbp, frames, strands, str, nBsp, seqI, seqII, 
              array3D(0, 0, 0), array3D(0, 0, 0), // BSP
              matrix(0, 0), array3D(0, 0, 0), // elasp
              // strucp, avstrp, ovStrucp, ovAvstrp
              array3D(0, 0, 0), matrix(0, 0),
              array3D(0, 0, 0), array3D(0, 0, 0),
              BPP, ovBpp);
    }
  }

  /* Write BPP parameters to file */
  proc writeBPP(ref BPP: [] real, ref ovBpp: [] real, seq: [] string, nbp: int,
                 numFrames: int, strandsType: int, isCircular: bool) {
    try! {
      var f = open("BPP.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln("BASE-PAIR PARAMETERS");
      writer.writeln("");
      writer.writeln(if isCircular then "CLOSED STRUCTURE"
                                    else "LINEAR STRUCTURE");
      writer.writeln(if strandsType == 2
                     then "DOUBLE-STRANDED STRUCTURE ANALYSED"
                     else "SINGLE-STRANDED STRUCTURE ANALYSED");
      writer.writef("BASE-PAIRS: %d\n", nbp);
      writer.writef("FRAMES:     %d\n", numFrames);
      writer.writeln("SEQUENCE:   ");
      writer.write("STRAND 1:   ");
      for i in 1..nbp do writer.write(seq[i]);
      writer.writeln();
      if strandsType == 2 {
        writer.write("STRAND 2:   ");
        for i in 1..nbp do writer.write(seq[2*nbp - i + 1]);
        writer.writeln();
      }
      writer.writeln("");
      writer.writeln(
        "First column averages, second column standard deviations"
      );
      writer.writef("%-10s %20s %20s %20s %20s %20s %20s\n", " base-pair",
                    "Shear", "Stretch", "Stagger", "Buckle", "Propeller",
                    "Opening");
      writer.writeln("-" * 130);
      
      for i in 1..nbp {
        var s2 = if strandsType == 2 then seq[2*nbp - i + 1] else "";
        var sep = if strandsType == 2 then "-" else " ";
        writer.writef("%6d %s %s%s%s", i, " ", seq[i], sep, s2);
        for p in 1..6 {
          writer.writef("%10.3f%10.3f", BPP[1, p, i], BPP[2, p, i]);
        }
        writer.writeln();
      }
      writer.writeln("-" * 130);
      writer.write(" AVG-STD= ");
      for p in 1..6 {
        writer.writef("%10.3f%10.3f", ovBpp[1, p], ovBpp[2, p]);
      }
      writer.writeln();
      writer.close();
      f.close();
    }
  }

  /* Write BSP parameters to file */
  proc writeBSP(ref BSP: [] real, ref ovBsp: [] real, seq: [] string, nbp: int, 
                numFrames: int, strandsType: int, isCircular: bool) {
    try! {
      var f = open("BSP.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln("BASE-STEP PARAMETERS");
      writer.writeln("");
      writer.writeln(if isCircular then "CLOSED STRUCTURE"
                                    else "LINEAR STRUCTURE");
      writer.writeln(if strandsType == 2
                     then "DOUBLE-STRANDED STRUCTURE ANALYSED"
                     else "SINGLE-STRANDED STRUCTURE ANALYSED");
      writer.writeln("BASE-PAIRS: ", nbp);
      writer.writeln("FRAMES:     ", numFrames);
      writer.write("SEQUENCE:   \nSTRAND 1:   ");
      for s in seq[1..nbp] do writer.write(s);
      writer.write("\nSTRAND 2:   ");
      if strandsType == 2 {
          for i in 1..nbp do writer.write(seq[2*nbp-i+1]);
      }
      writer.writeln("\n\nFirst column averages, " + 
                     "second column standard deviations");
      writer.writef(F_BSP_1, " base-step ", "Shift", "Slide", "Rise", "Tilt", 
                    "Roll", "Twist", "Bending");
      writer.writeln("-" * 130);
      
      var l = 0;
      var nSteps = if isCircular then nbp else nbp-1;
      for i in 1..nSteps {
          var w = if isCircular && i == nbp then 1 else i + 1;
          l += 1;
          var s2_1 = if strandsType == 2 then seq[2*nbp-w+1] else "#";
          var s2_2 = if strandsType == 2 then seq[2*nbp-i+1] else "#";
          writer.writef(F_BSP_2, i, "-",
                        w, " ", seq[i],
                        seq[w], "/", s2_1, s2_2);
          for p in 1..7 {
              writer.writef("%10.3f%10.3f", BSP[1, p, l], BSP[2, p, l]);
          }
          writer.writeln();
      }
      writer.writeln("-" * 130);
      writer.write(" AVG-STD= ");
      for p in 1..7 {
          writer.writef("%10.3f%10.3f", ovBsp[1, p, 1], ovBsp[2, p, 1]);
      }
      writer.writeln();
      writer.close();
      f.close();
    }
  }

  /* Write structural parameters to file */
  proc writeStructural(ref strucp: [] real, ref ovStrucp: [] real, 
                       ref avstrp: [] real, ref ovAvstrp: [] real, 
                       seq: [] string, nbp: int, numFrames: int, 
                       strandsType: int, isCircular: bool) {
    try! {
      var f = open("structural.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln("STRUCTURAL PARAMETERS");
      writer.writeln("");
      writer.writeln(if isCircular then "CLOSED STRUCTURE" 
                                    else "LINEAR STRUCTURE");
      writer.writeln(if strandsType == 2 
                     then "DOUBLE-STRANDED STRUCTURE ANALYSED" 
                     else "SINGLE-STRANDED STRUCTURE ANALYSED");
      writer.writeln("BASE-PAIRS: ", nbp);
      writer.writeln("FRAMES:     ", numFrames);
      writer.write("SEQUENCE:   \nSTRAND 1:   ");
      for s in seq[1..nbp] do writer.write(s);
      writer.write("\nSTRAND 2:   ");
      if strandsType == 2 {
          for i in 1..nbp do writer.write(seq[2*nbp-i+1]);
      }
      writer.writeln(
        "\n\nFirst column averages, second column standard deviations"
      );
      writer.writef(F_STRP_1,
                    " base-step ", "Shift", "Slide", "Rise", "Tilt",
                    "Roll", "Twist", "Bending", "Stiffness", "Energy");
      writer.writeln("-" * 130);
      
      var l = 0;
      var nSteps = if isCircular then nbp else nbp-1;
      for i in 1..nSteps {
          var w = if isCircular && i == nbp then 1 else i + 1;
          l += 1;
          var s2_1 = if strandsType == 2 then seq[2*nbp-w+1] else "#";
          var s2_2 = if strandsType == 2 then seq[2*nbp-i+1] else "#";
          writer.writef(F_STRP_2, i, "-", w, " ",
                      seq[i], seq[w], "/", s2_1, s2_2);
          for p in 1..11 {
              writer.writef("%10.3f%10.3f", strucp[1, p, l], strucp[2, p, l]);
          }
          for p in 1..3 {
              writer.writef("%10.3f", avstrp[p, l]);
          }
          writer.writeln();
      }
      writer.writeln("-" * 130);
      writer.write(" AVG-STD= ");
      for p in 1..11 {
          writer.writef("%10.3f%10.3f", ovStrucp[1, p, 1], ovStrucp[2, p, 1]);
      }
      for p in 1..3 {
          writer.writef("%10.3f%10.3f", ovAvstrp[1, p, 1], ovAvstrp[2, p, 1]);
      }
      writer.writeln();
      writer.close();
      f.close();
    }
  }

  /* Write elastic parameters to file */
  proc writeElasticParms(ref elasp: [] real, ref ovElasp: [] real,
                         seq: [] string,
                         nbp: int, numFrames: int, strandsType: int,
                         isCircular: bool) {
    try! {
      var f = open("elastic.out", ioMode.cw);
      var writer = f.writer();
      writer.writeln("ELASTIC PARAMETERS");
      writer.writeln("");
      writer.writeln(if isCircular 
      then "CLOSED STRUCTURE" 
      else "LINEAR STRUCTURE");
      writer.writeln(if strandsType == 2 
      then "DOUBLE-STRANDED STRUCTURE ANALYSED" 
      else "SINGLE-STRANDED STRUCTURE ANALYSED");
      writer.writeln("BASE-PAIRS: ", nbp);
      writer.writeln("FRAMES:     ", numFrames);
      writer.write("SEQUENCE:   \nSTRAND 1:   ");
      for s in seq[1..nbp] do writer.write(s);
      writer.write("\nSTRAND 2:   ");
      if strandsType == 2 {
          for i in 1..nbp do writer.write(seq[2*nbp-i+1]);
      }
      writer.writeln(
        "\n\nFirst column averages, second column standard deviations"
      );
      writer.writef(F_ELAP_1,
      " base-step ", "Shift", "Slide", "Rise", "Tilt", "Roll",
      "Twist", "Bending", "Stiffness", "Energy");
      writer.writeln("-" * 130);

      var l = 0;
      var nSteps = if isCircular then nbp else nbp-1;
      for i in 1..nSteps {
        var w = if isCircular && i == nbp then 1 else i + 1;
        l += 1;
        var s2_1 = if strandsType == 2 then seq[2*nbp-w+1] else "#";
        var s2_2 = if strandsType == 2 then seq[2*nbp-i+1] else "#";
        writer.writef(F_ELAP_2, i, "-", w, " ",
                    seq[i], seq[w], "/", s2_1, s2_2);
        for p in 1..13 {
            writer.writef("%20.3f", elasp[p, l]);
        }
        writer.writeln();
      }
      writer.writeln("-" * 130);
      writer.write(" AVG-STD= ");
      for p in 1..13 {
        writer.writef("%10.3f%10.3f", ovElasp[1, p, 1], ovElasp[2, p, 1]);
      }
      writer.writeln();
      writer.close();
      f.close();
    }
  }
}
