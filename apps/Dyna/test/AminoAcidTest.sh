#!/bin/bash

cat > md.in << EOF
&files
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd
         -x grp_one.crd
         -label GlyArg -n 8 }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd
         -label GlyArg_II -n 20 }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd
         -label GlyArg_III -n 20 }
  -x mdna.crd
  x_kind AMBER_CRD
  -o mdn.out
&end

&minimize
  cdcyc 20,  ncyc 40,  maxcyc 60,
  ntpr 1,
&end

&dynamics
  nstlim = 100000,  ntpr = 2500,  ntwx = 50000, dt = 1.0,
  ntt = 3, % nscm = 0,
  rigid_geom on,
  temperature = { tempi 100.0, temp0 300.0, -label GlyArg },
  temperature = { tempi 100.0, temp0 400.0, -label GlyArg_II },
  temperature = { tempi 300.0, temp0 200.0, -label GlyArg_III },
  tevo_start = 250, tevo_end = 750,
  tcache_depth 1,
&end

&solvent
  igb = 8,
&end

&restraint
  ensemble positions,
  mask '@N,CA,C,O & :1-2000',
  r1 = 0.0, r2 = 0.0, r3 = 0.0, r4 = 10.0, rk2 = 100.0, rk3 = 100.0,
  system all_possible,
&end

&report
  syntax = Matlab,
  energy total,
&end
EOF

${STORMM_BUILD}/apps/Dyna/dynamics.stormm.cuda -O -i md.in -except warn
