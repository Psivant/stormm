#!/bin/bash

cat > ffld.in << EOF
&files
  -sys { -p ${STORMM_SOURCE}/test/Topology/sulfonamide.top
         -c ${STORMM_SOURCE}/test/MoleculeFormat/sulfonamide_rots.sdf
         -label sulfonamide frame_end -1 -r sulfonamide_min.sdf r_kind SDF }

  -sys { -p ligand_11_5707.parm7
         -c ligand_11_5707.sdf frame_end -1
         -label ligand11 -r ligand11_min.sdf r_kind SDF }

&end

&ffrefine
  
&end

&minimize
  ncyc 50,  maxcyc 500,
&end

&restraint
  system sulfonamide,
  ensemble heavy_dihedrals,
  mask '@O1,S1,C2,C3',
  penalty 50.0, fbhw 0.0,
&end

&restraint
  system sulfonamide
  ensemble heavy_dihedrals
  penalty 20.0, fbhw 0.1,
&end

&restraint
  system sulfonamide
  ensemble prevent_hbonds
  penalty 16.0, proximity 3.2,
&end

&report
  sdf_item { -title etot  -label ALL -energy TOTAL_ENERGY }
  sdf_item { -title eptot -label ALL -energy POTENTIAL_ENERGY }
  sdf_item { -title edihe -label ALL -energy PROPER_DIHEDRAL }
  sdf_item { -title note  -label ALL -message "Hello world" }
  sdf_item { -title "S_N_bond" -parameter BOND -label ALL -typeI ab -typeJ ad }
&end
EOF

${STORMM_HOME}/apps/bin/ffrefine.stormm -O -i ffld.in -warn
