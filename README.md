___________________________________________________________________
                        \       _
			'       '
                \       ---------        /
		|----------SERRA----------|
                        ---------
                        |       |
                        -       -
  Chapel by: Alex Dukhan; Original implementation by: Victor Velasco with the guidance of Agnes Noy
___________________________________________________________________
___________________________________________________________________
                       OVERVIEW
                       --------

Serra is an analysis suite for determining the flexibility and structural 
parameters of nucleic acids from simulation data. It is composed of five 
main programs: SerraNA, SerraLINE, Analysis, Extract, and WrLINE.

The whole program calculates flexibility parameters from simulations of 
nucleic acid structures (circular or linear) at different length scales, 
by following the Length-Dependent Elastic Model (LDEM).

The methodology implemented allows the study of how global elastic
properties arise from the base-pair (bp) level, and the inference of elastic
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

To compile the entire suite:

`make all`

Individual components can be built using `mason build` within their 
respective directories or via the root Makefile (e.g., `make SerraNA`).
The compiled binaries will be located in `<Component>/target/debug/<Component>`.

___________________________________________________________________

                       CONFIG FILE SUPPORT
                       -------------------

All tools (SerraNA, SerraLINE, Extract, and Analysis) support an
optional `--configFile=<path>` flag that reads parameters from a plain
text file instead of (or in addition to) command-line flags.

**File format** (`key = value`, one pair per line):
- Lines beginning with `#` or `c` are comments and are ignored.
- Blank lines are ignored.
- Unknown keys are silently ignored.
- **CLI flags always take priority** over config file values, so you can
  mix both: put common settings in the file and override individual
  parameters on the command line as needed.

Example invocation:

```
./SerraNA --configFile=SerraNA.in
./SerraNA --configFile=SerraNA.in --structType=2   # override one value
```

Example `.in` files are provided alongside each tool's source code:

| Tool      | Example config               |
|-----------|------------------------------|
| SerraNA   | `SerraNA/src/SerraNA.in`     |
| SerraLINE | `SerraLINE/src/SerraLINE.in` |
| Extract   | `Extract/src/Extract.in`     |
| Analysis  | `Analysis/src/Analysis.in`   |


___________________________________________________________________

                       SerraNA
                       -------

First main process that carries the calculation of structural and 
flexible parameters for every possible pair of bp.
For running SerraNA, use the following flags:

```
./SerraNA --top=<topology_file> --traj=<trajectory_file> [options]
```

Or via a config file:

```
./SerraNA --configFile=SerraNA.in
```

Options:
  --strandsType=<1|2>   1 for single-stranded, 2 for double-stranded (default: 2)
  --structType=<1|2>    1 for linear, 2 for circular (default: 1)
  --configFile=<path>   Path to a key = value config file (see Config File Support)

SerraNA expects the topology and trajectory to be provided via command-line flags
or via a config file. The files can contain ions or other residues and SerraNA
will ignore them.
For linear NA, SerraNA ignores the two base-pairs at each end for avoiding end effects.

SerraNA creates 4 outputs:

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
  variance and partial variance of the end-to-end distance (angstroms).

___________________________________________________________________
                       Analysis
                       --------

Analysis is the second main process which analyzes the parameters 
calculated by SerraNA and estimates global elastic constants. 

For running Analysis, use the following flags:

```
./Analysis --elasFile=<elastic_file> --strucFile=<structural_file> [options]
```

Or via a config file:

```
./Analysis --configFile=Analysis.in
```

Options (Ranges for subfragment [a..b] and sublength [l1..l2]):
  --rA_a, --rA_b, --rA_l1, --rA_l2   (Persistence length)
  --rT_a, --rT_b, --rT_l1, --rT_l2   (Twist)
  --rS_a, --rS_b, --rS_l1, --rS_l2   (Stretch)
  --configFile=<path>                  Path to a key = value config file

It calculates global stretch modulus, torsion and persistence length 
for specified regions of the molecule.

___________________________________________________________________
                       Extract
                       -------

The Extract program processes SerraNA outputs (elastic and structural 
parameters), creating simple files ready to plot. You can filter 
particular sub-lengths or extract averages and standard deviations 
as a function of length.

For running Extract, use the following flags:

```
./Extract --fileIn=<input_file> --typeExt=<0|1> [options]
```

Or via a config file:

```
./Extract --configFile=Extract.in
```

Options:
  --typeExt=<0|1>   0 for sublength extraction, 1 for subfragment overall (default: 0)
  --sublength=<int> Length to extract (for typeExt=0)
  --a, --b          Range to extract (for typeExt=1)
  --configFile=<path>  Path to a key = value config file

___________________________________________________________________
                       SerraLINE
                       ---------

SerraLINE calculates bending angles, width, height, aspect ratio 
and deviation from planarity of DNA molecules using the global 
molecular contour WrLINE. 

The molecular contour WrLINE defines a point for each bp. SerraLINE 
can process closed (circular) and opened (linear) DNAs.

Bending angles are measured between two tangent vectors that can be 
separated by different number of bps (or points). Tangent vectors 
are constructed by the vectors joining two points that can be 
consecutive or arbitrarily separated.

SerraLINE can project the WrLINE contour to the plane that best fits 
the molecule or to the given region specified by the user. This 
projection method mimics experiments where structures are 
visualized on a 2D plane.

For running SerraLINE, use the following flags:

```
./SerraLINE --traj=<trajectory_file> [options]
```

Or via a config file:

```
./SerraLINE --configFile=SerraLINE.in
```

Options:
  --top=<topology_file>    Amber topology file
  --isCircular=<bool>      Set to true for circular DNA (default: false)
  --strandsType=<0|1|2>    0 for none, 1 for single, 2 for double (default: 1)
  --nbpArg=<int>           Number of base-pairs if no topology provided
  --fitplane=<bool>        Fit a plane to the molecule (default: true)
  --bpFitting=<indices>    Comma-separated indices for plane fitting
  --tLength=<int>          Tangent length (default: 1)
  --printProj=<bool>       Print projected coordinates (default: false)
  --xyzFormat=<bool>       Use xyz format for projected output (default: false)
  --configFile=<path>      Path to a key = value config file

___________________________________________________________________
                       WrLINE
                       ------

WrLINE extracts the helical axis and calculates 'writhe' from an 
AMBER trajectory of DNA. It provides the xyz trajectory 
of the helical axis, calculated writhe time series, and twist 
at each base-pair.

For running WrLINE, use the following flags:

./WrLINE --name=<jobname> --nbp=<n_basepairs> --nstep=<n_timesteps> [options]

Options:
  --top=<prmtop>    Topology file
  --traj=<mdcrd>    Trajectory file
  --doStrip=<bool>  Run stripC.sh (default: false)

___________________________________________________________________

                        LICENSE
                       ---------

    Serra is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    any later version.

    Serra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

___________________________________________________________________

