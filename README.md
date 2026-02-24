___________________________________________________________________
                        \       _
			'       '
                \       ---------        /
		|--------SERRA-DUO--------|
                        ---------
                        |       |
                        -       -
  by Victor Velasco with the guidance of Agnes Noy & parallelised by Alexander Dukhan
___________________________________________________________________
___________________________________________________________________
                       OVERVIEW
                       --------

SerraNA and SerraLINE are composed of one executable each.

The whole program calculates flexibility parameters from simulations of 
nucleic acid structures (circular or linear) at different length scales, 
by following the Length-Dependent Elastic Model (LDEM).

The methodology implemented, allows to study how global elastic
properties arise from the base-pair (bp) level, and to infer elastic
constants that describe the fragment's overall flexibility.
These elastic constants correspond to torsional, stretch modulus,
and the persistence length.

For more information check the manual and the paper:
  Velasco-Berrelleza, V. et al. SerraNA: a program to determine nucleic acids elasticity from simulation data. Phys. Chem. Chem. Phys. 22, 19254–19266 (2020).
This is the paper to cite. For more information, please visit:
https://agnesnoylab.wordpress.com/

___________________________________________________________________
                       REQUIREMENTS
                       ------------
To build from scratch:
- The chapel language compiler (chpl)
- The mason build tool
- CBLAS

___________________________________________________________________

                       COMPILATION
                       -----------

`mason build`

___________________________________________________________________

                       SerraNA
                       -------

First main process that carries the calculation of structural and 
flexible parameters flexibility for every possible pair of bp.
For running SerraNA type:

./SerraNA < s_NA.in

A trajectory file and a topology file in AMBER style format is needed 
(10F8.3 for the trajectory). The files can contain ions or other residues 
and SerraNA will ignore them.

s_NA.in is the input file that indicates: 
  1.- The path for topology and trajectory
  2.- If the structure is double-stranded ("2") or if it is 
      single-stranded ("1").
  3.- If the structure is linear ("1") or closed ("2"). 
      For linear NA, SerraNA ignores the two base-pairs at each end for
      avoiding end effects.

SerraNA creates 4 ouputs:

- BPP.out which have the base-pair parameters as they are calculated in 3DNA
  
- BSP.out which have the base-step parameters as they are calculated in 3DNA
  plus the total bending

- structural_parameters.out which have variables describing the geometry of the
  DNA molecule for all possible sub-fragments:
  1.- Added-rise, added-slide, added-rise, roll, tilt and twist as an extension
      of CEHS algorithm to sequences > 2 bp.
  2.- End-to-end distance and contour length
  3.- Bending, bending**2, and directional correlation (D correlation) 
  4.- From averaged structure: Bending (AVSTR B), bending**2 (AVSTR B**2) 
      and directional correlation (AVSTR D C) 

- elastic_parameters.out have elastic constants for Stretch (pN), Twist (nm), 
  Roll (nm), Tilt (nm), as well as their couplings (nm), together with the
  variance and partial variance  of the end-to-end distance (angstroms).

___________________________________________________________________
                       SerraLINE
                       --------

SerraLINE is a software that calculates bending angles, width, height,
aspect ratio and deviatin from planarity of DNA molecules using the 
global molecular contour WrLINE. 

The molecular contour WrLINE defines a point for each bp.

SerraLINE can process closed (circular) and opened (linear) DNAs.

Bending angles are measured between two tangent vectors that can be separated 
by different number of bps (or points).

Tangent vectors are constructed by the vectors joining two points that 
can be consecutive or arbitrarily separated.
This latter option is included to further approach microscopy imaging and 
analysis with a measurable limited resolution. 

SerraLINE can project the WrLINE contour to the plane that best fits the 
molecule or to the given region specified by the user.
This projection method mimics experiments where structures are 
visualized on a two dimensional plane.
Then, quantities such as width, height, aspect ratio (width/height)
and deviation from planarity can be calculated.

For both options (projection or not), bending angles are calculated 
following the same criteria.

___________________________________________________________________

                        LICENSE
                       ---------

    SerraDUO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    any later version.

    SerraDUO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

___________________________________________________________________
