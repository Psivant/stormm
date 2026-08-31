#!/bin/bash

cat > md.in << EOF
&files
  -sys { -p ${STORMM_SOURCE}/test/Topology/jac.top
         -c ${STORMM_SOURCE}/test/Trajectory/jac.inpcrd
         -label JAC -n 2 }
  -o md_jac.m
  -a hb_jac.m
&end

&dynamics
  nstlim = 5000000,  ntpr = 500,  ntwx = 10000, dt = 1.0,
  cut = 9.0,
  ntt = 0,
  rigid_geom on,
!  temperature = { tempi 100.0, temp0 300.0, -label JAC },
  tevo_start = 25, tevo_end = 75,
&end

&precision
  nonbonded single,
  valence single,
&end
&analysis

&end

&report
  syntax = Matlab,
  energy total,
  energy bond,
  energy angle,
  energy dihedral,
  energy electrostatic,
  energy vdw,
  energy elec_14,
  energy vdw_14,
  energy kinetic,
  ascii_salvage STARS,
&end
EOF

if [ ${1} ] ; then
  if [ ${1} == "CPU" ] ; then
    ${STORMM_BUILD}/apps/Dyna/dynamics.stormm -O -i md.in -except warn
  else
    ${STORMM_BUILD}/apps/Dyna/dynamics.stormm.cuda -O -i md.in -except warn
  fi
else
  ${STORMM_BUILD}/apps/Dyna/dynamics.stormm.cuda -O -i md.in -except warn
fi
