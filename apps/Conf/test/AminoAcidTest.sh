#!/bin/bash

echo "&files" > cgen.in
#for SYS in gly_lys gly_gly ala arg gly lys phe pro trp tyr gly_ala gly_arg gly_phe gly_pro \
#           gly_trp gly_tyr ; do
for SYS in ala arg gly phe pro trp tyr ; do
  echo "  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/${SYS}.top" >> cgen.in
  echo "         -c ${STORMM_SOURCE}/test/Namelists/coord/${SYS}.inpcrd " >> cgen.in
  echo "         -x conf_${SYS}.crd x_kind SDF -label ${SYS} }" >> cgen.in
done
for SYS in gly_lys gly_gly ; do
  echo "  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/${SYS}.top" >> cgen.in
  echo "         -c ${STORMM_SOURCE}/test/Namelists/coord/${SYS}.inpcrd " >> cgen.in
  echo "         -x conf_${SYS}.crd x_kind SDF -label ${SYS}_part2 }" >> cgen.in
done
cat >> cgen.in << EOF
  -x conf.crd
  x_kind amber_crd
&end

&conformer
  rotation_sample_count 3,
  trial_limit 100, local_trial_limit 200, final_states 6,
  core_mask { atoms "@N,CA,C & !(:ACE,NME)", rk3 4.0, grace 1.25 }
  effort LIGHT
&end

&solvent
  igb 8
&end

&minimize
  ncyc 25, cdcyc 0, maxcyc 25, ntpr = 1,
  clash_vdw_ratio 0.65,
&end

&random
  igseed 9183025
&end

&report
  sdf_item { -title total_energy    -label ALL -energy TOTAL_ENERGY }
  sdf_item { -title dihedral_energy -label ALL -energy PROPER_DIHEDRAL }
  sdf_item { -title elec_energy     -label ALL -energy ELECTROSTATIC_NONBONDED }
  sdf_item { -title elec_14_energy  -label ALL -energy ELECTROSTATIC_NEAR }
  report_width 99
  e_precision 2
  scope cluster_outlier
  syntax matlab
&end
EOF

cat >> cgen.in << EOF

EOF

if [ -e ${STORMM_BUILD}/apps/Conf/conformer.stormm.cuda ] ; then
  valgrind ${STORMM_BUILD}/apps/Conf/conformer.stormm.cuda -O -i cgen.in -warn
elif [ -e ${STORMM_BUILD}/apps/Conf/conformer.stormm ] ; then
  ${STORMM_BUILD}/apps/Conf/conformer.stormm -O -i cgen.in -t cgen.out -warn
fi
