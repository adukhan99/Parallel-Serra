/* This is a re-write of the SerraNA package developed
   by Victor Manuel Velasco Berrelleza and Agnes Noy. All Chapel
   code is written by Alexander Dukhan. The purpose of this project
   is both pedagogical and to permit highly parallel analysis
   of nucleotide structure. */

module SerraNA {

  use Parms;
  use Math;
  use LinearAlgebra;
  use CustomIO;
  use Functions;
  use IO;

  config const top: string = "";
  config const traj: string = "";
  // 1 for single-stranded, 2 for double-stranded
  config const strandsType: int = 2;
  config const structType: int = 1;  // 1 for linear, 2 for circular

  proc elasticF(IV: [1..4, 1..4] real, N: int) {
    var elastic_F: [1..10] real;
    elastic_F[1]  = 1.0e+23 * bKTpN * IV[4,4] * N:real * bnm;
    elastic_F[2]  = IV[3,3] * N:real * bnm * rad_to_deg * rad_to_deg;
    elastic_F[3]  = IV[2,2] * N:real * bnm * rad_to_deg * rad_to_deg;
    elastic_F[4]  = IV[1,1] * N:real * bnm * rad_to_deg * rad_to_deg;
    elastic_F[5]  = IV[4,3] * N:real * bnm * rad_to_deg / w_0;
    elastic_F[6]  = IV[4,2] * N:real * bnm * rad_to_deg / w_0;
    elastic_F[7]  = IV[4,1] * N:real * bnm * rad_to_deg / w_0;
    elastic_F[8]  = IV[3,2] * N:real * bnm * rad_to_deg * rad_to_deg;
    elastic_F[9]  = IV[3,1] * N:real * bnm * rad_to_deg * rad_to_deg;
    elastic_F[10] = IV[2,1] * N:real * bnm * rad_to_deg * rad_to_deg;
    return elastic_F;
  }

  proc main() {
    try! {
    if top == "" || traj == "" {
      writeln(
        "Topology and trajectory files must be specified using " +
        "--top=<file> and --traj=<file>");
      return;
    }

    if strandsType == 1 then writeln("Single stranded structure");
    else if strandsType == 2 then writeln("Double stranded structure");
    else 
      halt("Tell me if it is a single or double stranded structure (1 or 2)");

    if structType == 1 then writeln("Linear structure");
    else if structType == 2 then writeln("Circular structure");
    else 
      halt("Tell me if it is a circular or linear structure (1 or 2)");

    writeln("Topology file = ", top);
    writeln("Trajectory file = ", traj);

    writeln("Reading topology file");
    var (origNbp, nAtoms, box, origSeq, ringIndices) =
      topologyAmber(top, strandsType);

    writeln("Reading trajectory file");
    var isCircular = (structType == 2);
    var (coords, frames, nbp, seq) =
      coordinatesAmberCrd(
        traj, nAtoms, origNbp, box, strandsType, isCircular,
        ringIndices, origSeq);

    writeln("Calculating structural parameters");
    var nb = strandsType * nbp;
    var n_bsp =
      if structType == 2 then nbp * (nbp - 1) else nbp * (nbp - 1) / 2;

    var R: [1..3, 1..3, 1..nb] real;
    var O: [1..3, 1..nb] real;
    var Rmbt: [1..3, 1..3, 1..nbp] real;
    var Ombt: [1..3, 1..nbp] real;

    var BSP: [1..9, 1..n_bsp, 1..frames] real;
    var BPP: [1..6, 1..nbp, 1..frames] real;
    var E_E_dist: [1..n_bsp, 1..frames] real;
    var C_length: [1..n_bsp, 1..frames] real;
    var added_bsp: [1..3, 1..n_bsp, 1..frames] real;

    var V: [1..4, 1..4, 1..n_bsp] real;
    var V__1: [1..4, 1..4, 1..n_bsp] real;
    var F: [1..10, 1..n_bsp] real;
    var Ad2: [1..n_bsp] real;

    var i_case = 0;

    select (structType, strandsType) {
      when (1, 2) do i_case = 1;
      when (1, 1) do i_case = 2;
      when (2, 2) do i_case = 3;
      when (2, 1) do i_case = 4;
      otherwise do halt("Error in determining case");
    }

    for k in 1..frames {
      var l = 1;
      for i in 1..nb {
        var nAuth: int;
        if seq[i] == "G" || seq[i] == "A" then nAuth = n_R;
        else nAuth = n_Y;

        var Ex: [1..3, 1..nAuth] real = coords[1..3, l..l+nAuth-1, k];
        var R_O: ( [1..3, 1..3] real, [1..3] real );

        select seq[i] {
          when "G" do R_O = getRotationROriginO(Ex, G_b);
          when "A" do R_O = getRotationROriginO(Ex, A_b);
          when "C" do R_O = getRotationROriginO(Ex, C_b);
          when "T" do R_O = getRotationROriginO(Ex, T_b);
          when "U" do R_O = getRotationROriginO(Ex, U_b);
        }

        R[1..3, 1..3, i] = R_O(0);
        O[1..3, i] = R_O(1);
        l += nAuth;
      }

      if strandsType == 2 {
        for i in nbp+1..nb {
          R[1..3, 2, i] = -R[1..3, 2, i];
          R[1..3, 3, i] = -R[1..3, 3, i];
        }

        for i in 1..nbp {
          var BPP_T_O =
            basepairParameters(
              O[1..3, i], O[1..3, nb+1-i],
              R[1..3, 1..3, i], R[1..3, 1..3, nb+1-i]);
          BPP[1..6, i, k] = BPP_T_O(0);
          Rmbt[1..3, 1..3, i] = BPP_T_O(1);
          Ombt[1..3, i] = BPP_T_O(2);
        }
      }

      var bspIdx = 0;
      var t_twist = 0.0;

      if i_case == 1 {
        for i in 1..nbp-1 {
          bspIdx += 1;
          var res =
            basestepParameters(
              Ombt[1..3, i], Ombt[1..3, i+1],
              Rmbt[1..3, 1..3, i], Rmbt[1..3, 1..3, i+1], t_twist);
          BSP[1..9, bspIdx, k] = res;
          E_E_dist[bspIdx, k] = absv(Ombt[1..3, i+1] - Ombt[1..3, i]);
          C_length[bspIdx, k] = E_E_dist[bspIdx, k];
          added_bsp[1..3, bspIdx, k] = BSP[1..3, i, k];
        }
        for j in 2..nbp-1 {
          for i in 1..nbp-j {
            bspIdx += 1;
            t_twist = BSP[6, bspIdx-nbp+j-1, k];
            var res =
              basestepParameters(
                Ombt[1..3, i], Ombt[1..3, i+j],
                Rmbt[1..3, 1..3, i], Rmbt[1..3, 1..3, i+j], t_twist);
            BSP[1..9, bspIdx, k] = res;
            E_E_dist[bspIdx, k] = absv(Ombt[1..3, i+j] - Ombt[1..3, i]);
            // index check? i+j-1 was the 1bp length step
            C_length[bspIdx, k] =
              C_length[bspIdx-nbp+j-1, k] + C_length[i+j-1, k];
            added_bsp[1..3, bspIdx, k] =
              added_bsp[1..3, bspIdx-nbp+j-1, k] + BSP[1..3, i+j-1, k];
          }
        }
      } else if i_case == 2 {
        for i in 1..nbp-1 {
          bspIdx += 1;
          var res =
            basestepParameters(
              O[1..3, i], O[1..3, i+1],
              R[1..3, 1..3, i], R[1..3, 1..3, i+1], t_twist);
          BSP[1..9, bspIdx, k] = res;
          E_E_dist[bspIdx, k] = absv(O[1..3, i+1] - O[1..3, i]);
          C_length[bspIdx, k] = E_E_dist[bspIdx, k];
          added_bsp[1..3, bspIdx, k] = BSP[1..3, i, k];
        }
        for j in 2..nbp-1 {
          for i in 1..nbp-j {
            bspIdx += 1;
            t_twist = BSP[6, bspIdx-nbp+j-1, k];
            var res =
              basestepParameters(
                O[1..3, i], O[1..3, i+j],
                R[1..3, 1..3, i], R[1..3, 1..3, i+j], t_twist);
            BSP[1..9, bspIdx, k] = res;
            E_E_dist[bspIdx, k] = absv(O[1..3, i+j] - O[1..3, i]);
            C_length[bspIdx, k] =
              C_length[bspIdx-nbp+j-1, k] + C_length[i+j-1, k];
            added_bsp[1..3, bspIdx, k] =
              added_bsp[1..3, bspIdx-nbp+j-1, k] + BSP[1..3, i+j-1, k];
          }
        }
      } else 
      if i_case == 3 {
        for i in 1..nbp-1 {
          bspIdx += 1;
          var w = i+1;
          var res =
            basestepParameters(
              Ombt[1..3, i], Ombt[1..3, w],
              Rmbt[1..3, 1..3, i], Rmbt[1..3, 1..3, w], t_twist);
          BSP[1..9, bspIdx, k] = res;
          E_E_dist[bspIdx, k] = absv(Ombt[1..3, w] - Ombt[1..3, i]);
          C_length[bspIdx, k] = E_E_dist[bspIdx, k];
          added_bsp[1..3, bspIdx, k] = BSP[1..3, i, k];
        }
        for i in nbp..nbp {
          var w = i+1-nbp;
          bspIdx += 1;
          var res =
            basestepParameters(
              Ombt[1..3, i], Ombt[1..3, w],
              Rmbt[1..3, 1..3, i], Rmbt[1..3, 1..3, w], t_twist);
          BSP[1..9, bspIdx, k] = res;
          E_E_dist[bspIdx, k] = absv(Ombt[1..3, w] - Ombt[1..3, i]);
          C_length[bspIdx, k] = E_E_dist[bspIdx, k];
          added_bsp[1..3, bspIdx, k] = BSP[1..3, i, k];
        }
        for j in 2..nbp-1 {
          for i in 1..nbp-j {
            var w = i+j;
            bspIdx += 1;
            t_twist = BSP[6, bspIdx-nbp, k];
            var res =
              basestepParameters(
                Ombt[1..3, i], Ombt[1..3, w],
                Rmbt[1..3, 1..3, i], Rmbt[1..3, 1..3, w], t_twist);
            BSP[1..9, bspIdx, k] = res;
            E_E_dist[bspIdx, k] = absv(Ombt[1..3, w] - Ombt[1..3, i]);
            C_length[bspIdx, k] =
              C_length[bspIdx-nbp, k] + C_length[w-1, k];
            added_bsp[1..3, bspIdx, k] =
              added_bsp[1..3, bspIdx-nbp, k] + BSP[1..3, w-1, k];
          }
          var i2 = nbp-j+1;
          var w2 = 1;
          bspIdx += 1;
          t_twist = BSP[6, bspIdx-nbp, k];
          var res2 =
            basestepParameters(
              Ombt[1..3, i2], Ombt[1..3, w2],
              Rmbt[1..3, 1..3, i2], Rmbt[1..3, 1..3, w2], t_twist);
          BSP[1..9, bspIdx, k] = res2;
          E_E_dist[bspIdx, k] = absv(Ombt[1..3, w2] - Ombt[1..3, i2]);
          C_length[bspIdx, k] = C_length[bspIdx-nbp, k] + C_length[nbp, k];
          added_bsp[1..3, bspIdx, k] =
            added_bsp[1..3, bspIdx-nbp, k] + BSP[1..3, nbp, k];
          for i3 in nbp-j+2..nbp {
            var w3 = i3+j-nbp;
            bspIdx += 1;
            t_twist = BSP[6, bspIdx-nbp, k];
            var res3 =
              basestepParameters(
                Ombt[1..3, i3], Ombt[1..3, w3],
                Rmbt[1..3, 1..3, i3], Rmbt[1..3, 1..3, w3], t_twist);
            BSP[1..9, bspIdx, k] = res3;
            E_E_dist[bspIdx, k] = absv(Ombt[1..3, w3] - Ombt[1..3, i3]);
            C_length[bspIdx, k] =
              C_length[bspIdx-nbp, k] + C_length[w3-1, k];
            added_bsp[1..3, bspIdx, k] =
              added_bsp[1..3, bspIdx-nbp, k] + BSP[1..3, w3-1, k];
          }
        }
      } else if i_case == 4 {
        for i in 1..nbp-1 {
          bspIdx += 1;
          var w = i+1;
          var res =
            basestepParameters(
              O[1..3, i], O[1..3, w],
              R[1..3, 1..3, i], R[1..3, 1..3, w], t_twist);
          BSP[1..9, bspIdx, k] = res;
          E_E_dist[bspIdx, k] = absv(O[1..3, w] - O[1..3, i]);
          C_length[bspIdx, k] = E_E_dist[bspIdx, k];
          added_bsp[1..3, bspIdx, k] = BSP[1..3, i, k];
        }
        for i in nbp..nbp {
          var w = i+1-nbp;
          bspIdx += 1;
          var res =
            basestepParameters(
              O[1..3, i], O[1..3, w],
              R[1..3, 1..3, i], R[1..3, 1..3, w], t_twist);
          BSP[1..9, bspIdx, k] = res;
          E_E_dist[bspIdx, k] = absv(O[1..3, w] - O[1..3, i]);
          C_length[bspIdx, k] = E_E_dist[bspIdx, k];
          added_bsp[1..3, bspIdx, k] = BSP[1..3, i, k];
        }
        for j in 2..nbp-1 {
          for i in 1..nbp-j {
            var w = i+j;
            bspIdx += 1;
            t_twist = BSP[6, bspIdx-nbp, k];
            var res =
              basestepParameters(
                O[1..3, i], O[1..3, w],
                R[1..3, 1..3, i], R[1..3, 1..3, w], t_twist);
            BSP[1..9, bspIdx, k] = res;
            E_E_dist[bspIdx, k] = absv(O[1..3, w] - O[1..3, i]);
            C_length[bspIdx, k] =
              C_length[bspIdx-nbp, k] + C_length[w-1, k];
            added_bsp[1..3, bspIdx, k] =
              added_bsp[1..3, bspIdx-nbp, k] + BSP[1..3, w-1, k];
          }
          var i2 = nbp-j+1;
          var w2 = 1;
          bspIdx += 1;
          t_twist = BSP[6, bspIdx-nbp, k];
          var res2 =
            basestepParameters(
              O[1..3, i2], O[1..3, w2],
              R[1..3, 1..3, i2], R[1..3, 1..3, w2], t_twist);
          BSP[1..9, bspIdx, k] = res2;
          E_E_dist[bspIdx, k] = absv(O[1..3, w2] - O[1..3, i2]);
          C_length[bspIdx, k] = C_length[bspIdx-nbp, k] + C_length[nbp, k];
          added_bsp[1..3, bspIdx, k] =
            added_bsp[1..3, bspIdx-nbp, k] + BSP[1..3, nbp, k];
          for i3 in nbp-j+2..nbp {
            var w3 = i3+j-nbp;
            bspIdx += 1;
            t_twist = BSP[6, bspIdx-nbp, k];
            var res3 =
              basestepParameters(
                O[1..3, i3], O[1..3, w3],
                R[1..3, 1..3, i3], R[1..3, 1..3, w3], t_twist);
            BSP[1..9, bspIdx, k] = res3;
            E_E_dist[bspIdx, k] = absv(O[1..3, w3] - O[1..3, i3]);
            C_length[bspIdx, k] =
              C_length[bspIdx-nbp, k] + C_length[w3-1, k];
            added_bsp[1..3, bspIdx, k] =
              added_bsp[1..3, bspIdx-nbp, k] + BSP[1..3, w3-1, k];
          }
        }
      }
    }

    if frames > 1 {
      writeln("Calculating elastic parameters");
      for j in 1..nbp-1 {
        var w =
          if structType == 2 then (if j < nbp - j then j else nbp - j) else j;
        var limit = if structType == 2 then nbp else nbp - j;
        for l in 1..limit {
          V[1..4, 1..4, l] = deformationCovariance(
                               BSP[5, l, 1..frames],  // Roll
                               BSP[4, l, 1..frames],  // Tilt
                               BSP[6, l, 1..frames],  // Twist
                               E_E_dist[l, 1..frames] // Stretch
                             );
          V__1[1..4, 1..4, l] = inverseMatrixAnalytic4X4(V[1..4, 1..4, l]);
          F[1..10, l] = elasticF(V__1[1..4, 1..4, l], w);
          Ad2[l] = dynamicPersistenceLength2(F[4, l], F[3, l]);
        }
      }
    } else {
      writeln("Only one frame read");
      writeln("Cannot calculate elastic parameters");
    }

    writeln("Calculating averages and standard deviations");

    var avstd_BPP: [1..2, 1..6, 1..nbp] real;
    var avstd_BSP: [1..2, 1..9, 1..n_bsp] real;
    var avstd_E_E: [1..2, 1..n_bsp] real;
    var avstd_C_l: [1..2, 1..n_bsp] real;
    var avstd_added: [1..2, 1..3, 1..n_bsp] real;

    if strandsType == 2 {
      for i in 1..nbp {
        for j in 1..6 {
          var res = averageStd(BPP[j, i, 1..frames]);
          avstd_BPP[1, j, i] = res[1];
          avstd_BPP[2, j, i] = res[2];
        }
      }
    }

    for i in 1..n_bsp {
      for j in 1..9 {
        var res = averageStd(BSP[j, i, 1..frames]);
        avstd_BSP[1, j, i] = res[1];
        avstd_BSP[2, j, i] = res[2];
      }
      var resE = averageStd(E_E_dist[i, 1..frames]);
      avstd_E_E[1, i] = resE[1]; 
      avstd_E_E[2, i] = resE[2];
      var resC = averageStd(C_length[i, 1..frames]);
      avstd_C_l[1, i] = resC[1]; 
      avstd_C_l[2, i] = resC[2];
      for j in 1..3 {
        var resA = averageStd(added_bsp[j, i, 1..frames]);
        avstd_added[1, j, i] = resA[1];
        avstd_added[2, j, i] = resA[2];
      }
    }

    writeln("Calculating average structure parameters");
    var avstr_BSP: [1..9, 1..n_bsp] real;
    var avRmbt: [1..3, 1..3, 1..nbp] real;
    var avOmbt: [1..3, 1..nbp] real;

    var Tg_rg = reverseAlgorithmBaseTriads(avstd_BSP[1, 1..6, 1..nbp]);
    avRmbt = Tg_rg(0);
    avOmbt = Tg_rg(1);

    var bspIdx = 0;
    var t_twist = 0.0;
    if structType != 2 {
      for i in 1..nbp-1 {
        bspIdx += 1;
        var res =
          basestepParameters(
            avOmbt[1..3, i], avOmbt[1..3, i+1],
            avRmbt[1..3, 1..3, i], avRmbt[1..3, 1..3, i+1], t_twist);
        avstr_BSP[1..9, bspIdx] = res;
      }
      for j in 2..nbp-1 {
        for i in 1..nbp-j {
          bspIdx += 1;
          t_twist = avstr_BSP[6, bspIdx-nbp+j-1];
          var res =
            basestepParameters(
              avOmbt[1..3, i], avOmbt[1..3, i+j],
              avRmbt[1..3, 1..3, i], avRmbt[1..3, 1..3, i+j], t_twist);
          avstr_BSP[1..9, bspIdx] = res;
        }
      }
    } else {
      for i in 1..nbp-1 {
        bspIdx += 1;
        var res =
          basestepParameters(
            avOmbt[1..3, i], avOmbt[1..3, i+1],
            avRmbt[1..3, 1..3, i], avRmbt[1..3, 1..3, i+1], t_twist);
        avstr_BSP[1..9, bspIdx] = res;
      }
      var resLoop =
        basestepParameters(
          avOmbt[1..3, nbp], avOmbt[1..3, 1],
          avRmbt[1..3, 1..3, nbp], avRmbt[1..3, 1..3, 1], t_twist);
      bspIdx += 1;
      avstr_BSP[1..9, bspIdx] = resLoop;

      for j in 2..nbp-1 {
        for i in 1..nbp-j {
          var w = i+j;
          bspIdx += 1;
          t_twist = avstr_BSP[6, bspIdx-nbp];
          var res =
            basestepParameters(
              avOmbt[1..3, i], avOmbt[1..3, w],
              avRmbt[1..3, 1..3, i], avRmbt[1..3, 1..3, w], t_twist);
          avstr_BSP[1..9, bspIdx] = res;
        }
        var i2 = nbp-j+1;
        var w2 = 1;
        bspIdx += 1;
        t_twist = avstr_BSP[6, bspIdx-nbp];
        var res2 =
          basestepParameters(
            avOmbt[1..3, i2], avOmbt[1..3, w2],
            avRmbt[1..3, 1..3, i2], avRmbt[1..3, 1..3, w2], t_twist);
        avstr_BSP[1..9, bspIdx] = res2;

        for i3 in nbp-j+2..nbp {
          var w3 = i3+j-nbp;
          bspIdx += 1;
          t_twist = avstr_BSP[6, bspIdx-nbp];
          var res3 =
            basestepParameters(
              avOmbt[1..3, i3], avOmbt[1..3, w3],
              avRmbt[1..3, 1..3, i3], avRmbt[1..3, 1..3, w3], t_twist);
          avstr_BSP[1..9, bspIdx] = res3;
        }
      }
    }

    writeln("Calculating overall parameters");
    var OV_BPP: [1..2, 1..6] real;
    var OV_BSP: [1..2, 1..9, 1..nbp-1] real;
    var OV_E_E: [1..2, 1..nbp-1] real;
    var OV_C_l: [1..2, 1..nbp-1] real;
    var OV_added: [1..2, 1..3, 1..nbp-1] real;
    var OV_V_E_E: [1..2, 1..nbp-1] real;
    var OV_pV_E_E: [1..2, 1..nbp-1] real;
    var OV_F: [1..2, 1..10, 1..nbp-1] real;
    var OV_Ad2: [1..2, 1..nbp-1] real;
    var OV_BSP_avstr: [1..2, 1..9, 1..nbp-1] real;

    if strandsType == 2 {
      for i in 1..6 {
        var res = averageStd(avstd_BPP[1, i, 1..nbp]);
        OV_BPP[1, i] = res[1];
        OV_BPP[2, i] = res[2];
      }
    }

    var idxL = 1;
    for j in 1..nbp-1 {
      var limit = if structType == 2 then nbp else nbp - j;
      for i in 1..9 {
        var res = averageStd(avstd_BSP[1, i, idxL..idxL+limit-1]);
        OV_BSP[1, i, j] = res[1];
        OV_BSP[2, i, j] = res[2];
      }
      var resE = averageStd(avstd_E_E[1, idxL..idxL+limit-1]);
      OV_E_E[1, j] = resE[1];
      OV_E_E[2, j] = resE[2];
      var resC = averageStd(avstd_C_l[1, idxL..idxL+limit-1]);
      OV_C_l[1, j] = resC[1];
      OV_C_l[2, j] = resC[2];
      for i in 1..3 {
        var resA = averageStd(avstd_added[1, i, idxL..idxL+limit-1]);
        OV_added[1, i, j] = resA[1];
        OV_added[2, i, j] = resA[2];
      }
      idxL += limit;
    }

    if frames > 1 {
      idxL = 1;
      for j in 1..nbp-1 {
        var limit = if structType == 2 then nbp else nbp - j;
        var resV = averageStd(V[4, 4, idxL..idxL+limit-1]);
        OV_V_E_E[1, j] = resV[1];
        OV_V_E_E[2, j] = resV[2];

        var v1Arr: [1..limit] real;
        for i in 1..limit do v1Arr[i] = 1.0 / V__1[4, 4, idxL+i-1];
        var resPv = averageStd(v1Arr);
        OV_pV_E_E[1, j] = resPv[1];
        OV_pV_E_E[2, j] = resPv[2];

        for i in 1..10 {
          var resF = averageStd(F[i, idxL..idxL+limit-1]);
          OV_F[1, i, j] = resF[1];
          OV_F[2, i, j] = resF[2];
        }
        var resAd2 = averageStd(Ad2[idxL..idxL+limit-1]);
        OV_Ad2[1, j] = resAd2[1];
        OV_Ad2[2, j] = resAd2[2];
        idxL += limit;
      }
    }

    idxL = 1;
    for j in 1..nbp-1 {
      var limit = if structType == 2 then nbp else nbp - j;
      for i in 1..9 {
        var res = averageStd(avstr_BSP[i, idxL..idxL+limit-1]);
        OV_BSP_avstr[1, i, j] = res[1];
        OV_BSP_avstr[2, i, j] = res[2];
      }
      idxL += limit;
    }

    writeln("Writing output files");
    if strandsType == 2 {
      writeBPP(
        avstd_BPP, OV_BPP, seq, nbp, frames, strandsType, structType==2);
    }
    writeBSP(
      avstd_BSP, OV_BSP, seq, nbp, frames, strandsType, structType==2);

    // Build strucp [1..2, 1..11, 1..n_bsp]:
    //   params 1..9  = averaged base-step parameters (Shift..Bending²)
    //   param  10    = end-to-end distance
    //   param  11    = contour length
    var strucp: [1..2, 1..11, 1..n_bsp] real;
    for i in 1..n_bsp {
      for j in 1..9 {
        strucp[1, j, i] = avstd_BSP[1, j, i];
        strucp[2, j, i] = avstd_BSP[2, j, i];
      }
      strucp[1, 10, i] = avstd_E_E[1, i];
      strucp[2, 10, i] = avstd_E_E[2, i];
      strucp[1, 11, i] = avstd_C_l[1, i];
      strucp[2, 11, i] = avstd_C_l[2, i];
    }
    var ovStrucp: [1..2, 1..11, 1..nbp-1] real;
    for i in 1..nbp-1 {
      for j in 1..9 {
        ovStrucp[1, j, i] = OV_BSP[1, j, i];
        ovStrucp[2, j, i] = OV_BSP[2, j, i];
      }
      ovStrucp[1, 10, i] = OV_E_E[1, i];
      ovStrucp[2, 10, i] = OV_E_E[2, i];
      ovStrucp[1, 11, i] = OV_C_l[1, i];
      ovStrucp[2, 11, i] = OV_C_l[2, i];
    }
    // avstrp[1..3, 1..n_bsp]: Shift, Slide, Rise from the average-structure BSP
    // (rows 1-3 of avstr_BSP: Shift, Slide, Rise).

    var avstrp: [1..3, 1..n_bsp] real;
    avstrp[1..3, 1..n_bsp] = avstr_BSP[1..3, 1..n_bsp];

    // ovAvstrp[1..2, 1..3, 1..nbp-1]: mean and std of avstrp over each
    // separation length.
    var ovAvstrp: [1..2, 1..3, 1..nbp-1] real;
    for i in 1..3 do
      for j in 1..nbp-1 do
        ovAvstrp[1, i, j] = OV_BSP_avstr[1, i, j];
    writeStructural(
      strucp, ovStrucp, avstrp, ovAvstrp, seq, nbp, frames,
      strandsType, structType==2);

    if frames > 1 {
      var elasp: [1..13, 1..n_bsp] real;
      for i in 1..n_bsp {
        for j in 1..10 do elasp[j, i] = F[j, i];
        elasp[11, i] = Ad2[i];
      }
      var ovElasp: [1..2, 1..13, 1..nbp-1] real;
      writeElasticParms(
        elasp, ovElasp, seq, nbp, frames, strandsType, structType==2);
    }
    }
  }
}
