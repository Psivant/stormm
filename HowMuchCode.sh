#!/bin/bash

# Remove backup files created by emacs
rm *~ */*~ */*/*~ */*/*/*~

if [ -e "lines of code" ] ; then
  rm "lines of code"
fi
if [ -e "total lines of code" ] ; then
  rm "total lines of code"
fi
for DIR in src test ; do
  touch ${DIR}/UnitTesting/nothing_anyone_wants.h
  touch ${DIR}/UnitTesting/nothing_anyone_wants.cpp
  touch ${DIR}/UnitTesting/nothing_anyone_wants.tpp
  touch ${DIR}/UnitTesting/nothing_anyone_wants.cu
  touch ${DIR}/UnitTesting/nothing_anyone_wants.cuh
  touch ${DIR}/UnitTesting/nothing_anyone_wants.cui
  for EXT in h cpp tpp cu cuh cui ; do
    cat ${DIR}/*/*.${EXT} >> "lines of code"
    cat ${DIR}/*/*.${EXT} >> "total lines of code"
  done
  rm ${DIR}/UnitTesting/nothing_anyone_wants.h
  rm ${DIR}/UnitTesting/nothing_anyone_wants.cpp
  rm ${DIR}/UnitTesting/nothing_anyone_wants.tpp
  rm ${DIR}/UnitTesting/nothing_anyone_wants.cu
  rm ${DIR}/UnitTesting/nothing_anyone_wants.cuh
  rm ${DIR}/UnitTesting/nothing_anyone_wants.cui
  grep -v "//" "lines of code" > code_dump_file.tmp
  grep ";.*//" "lines of code" >> code_dump_file.tmp
  mv code_dump_file.tmp "lines of code"
  TT=`wc -l "lines of code"`
  echo "${TT} in ${DIR}"
  rm "lines of code"
done

for DIR in benchmark ; do
  touch ${DIR}/ForceAccumulation/nothing_anyone_wants.h
  touch ${DIR}/ForceAccumulation/nothing_anyone_wants.cpp
  touch ${DIR}/ForceAccumulation/nothing_anyone_wants.tpp
  touch ${DIR}/ForceAccumulation/nothing_anyone_wants.cu
  touch ${DIR}/ForceAccumulation/nothing_anyone_wants.cuh
  touch ${DIR}/ForceAccumulation/nothing_anyone_wants.cui
  for EXT in h cpp tpp cu cuh cui ; do
    cat ${DIR}/*/*.${EXT} >> "lines of code"
    cat ${DIR}/*/*.${EXT} >> "total lines of code"
  done
  rm ${DIR}/ForceAccumulation/nothing_anyone_wants.h
  rm ${DIR}/ForceAccumulation/nothing_anyone_wants.cpp
  rm ${DIR}/ForceAccumulation/nothing_anyone_wants.tpp
  rm ${DIR}/ForceAccumulation/nothing_anyone_wants.cu
  rm ${DIR}/ForceAccumulation/nothing_anyone_wants.cuh
  rm ${DIR}/ForceAccumulation/nothing_anyone_wants.cui
  grep -v "//" "lines of code" > code_dump_file.tmp
  grep ";.*//" "lines of code" >> code_dump_file.tmp
  mv code_dump_file.tmp "lines of code"
  TT=`wc -l "lines of code"`
  echo "${TT} in ${DIR}"
  rm "lines of code"
done

for DIR in apps ; do
  touch ${DIR}/Conf/src/nothing_anyone_wants.h
  touch ${DIR}/Conf/src/nothing_anyone_wants.cpp
  touch ${DIR}/Conf/src/nothing_anyone_wants.tpp
  touch ${DIR}/Conf/src/nothing_anyone_wants.cu
  touch ${DIR}/Conf/src/nothing_anyone_wants.cuh
  touch ${DIR}/Conf/src/nothing_anyone_wants.cui
  for EXT in h cpp tpp cu cuh cui ; do
    cat ${DIR}/*/*/*.${EXT} >> "lines of code"
    cat ${DIR}/*/*/*.${EXT} >> "total lines of code"
  done
  rm ${DIR}/Conf/src/nothing_anyone_wants.h
  rm ${DIR}/Conf/src/nothing_anyone_wants.cpp
  rm ${DIR}/Conf/src/nothing_anyone_wants.tpp
  rm ${DIR}/Conf/src/nothing_anyone_wants.cu
  rm ${DIR}/Conf/src/nothing_anyone_wants.cuh
  rm ${DIR}/Conf/src/nothing_anyone_wants.cui
  grep -v "//" "lines of code" > code_dump_file.tmp
  grep ";.*//" "lines of code" >> code_dump_file.tmp
  mv code_dump_file.tmp "lines of code"
  TT=`wc -l "lines of code"`
  echo "${TT} in ${DIR}"
  rm "lines of code"
done

TT=`wc -l "total lines of code"`
echo "${TT}"
grep -v "//" "total lines of code" > code_dump_file.tmp
grep ";.*//" "total lines of code" >> code_dump_file.tmp
mv code_dump_file.tmp "lines of code excluding documentation"
TT=`wc -l "lines of code excluding documentation"`
echo "${TT}"
rm "total lines of code"
rm "lines of code excluding documentation"
