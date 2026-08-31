#!/bin/bash

cat > md.in << EOF
&files
#  -sys { -p ${STORMM_SOURCE}/test/Topology/tip3p.top
#         -c ${STORMM_SOURCE}/test/Trajectory/tip3p.inpcrd
#         -label TIP3P -n 1 }
#  -sys { -p ${STORMM_SOURCE}/test/Topology/tip4p.top
#         -c ${STORMM_SOURCE}/test/Trajectory/tip4p.inpcrd
#         -label TIP4P -n 1 }
#  -sys { -p ${STORMM_SOURCE}/test/Topology/ubiquitin.top
#         -c ${STORMM_SOURCE}/test/Trajectory/ubiquitin.inpcrd
#         -label UBI -n 2 }
#  -sys { -p ${STORMM_SOURCE}/test/Topology/drug_example.top
#         -c ${STORMM_SOURCE}/test/Trajectory/drug_example.inpcrd
#         -label DRUG -n 3 }
  -sys { -p ${STORMM_SOURCE}/test/Topology/jac.top
         -c ${STORMM_SOURCE}/test/Trajectory/jac.inpcrd
         -label JAC -n 8 }
  -o mdn.out
&end

&dynamics
  nstlim = 250000,  ntpr = 5000,  ntwx = 25000, dt = 1.0,
  cut = 9.0,
  ntt = 0,
  rigid_geom off,
#  temperature = { tempi 100.0, temp0 300.0, -label TIP3P },
#  temperature = { tempi 100.0, temp0 400.0, -label TIP4P },
#  temperature = { tempi 300.0, temp0 300.0, -label UBI },
#  temperature = { tempi 300.0, temp0 300.0, -label DRUG },
#  tevo_start = 25, tevo_end = 75,
#  tcache_depth 1,
&end

&precision
#  globalpos_bits 48, velocity_bits 60, force_bits 48
  nonbonded single,
  valence single,
&end

!&debug
!  max_reports = 1000,
!  forces,
!  force_threshold = 250.0,
!  force_trigger, interval_trigger 100,
!  ngbr_placement
!&end

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
!  progbar none,
&end
EOF

if [ ${1} ] ; then
  if [ ${1} == "CPU" ] ; then
    ${STORMM_BUILD}/apps/Dyna/dynamics.stormm -O -i md.in -except warn
  else
    ${STORMM_BUILD}/apps/Dyna/dynamics.stormm.cuda -O -i md.in -except warn
  fi
else
  #compute-sanitizer
  ${STORMM_BUILD}/apps/Dyna/dynamics.stormm.cuda -O -i md.in -except warn
fi
