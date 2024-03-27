#!/bin/bash

if [ ! ${1} ] || [ ! ${2} ] ; then
  echo "Usage: ./ChangeAllSuchWords <original term> <replacement term> <authorization>"
  echo
  echo "  If the authorization is not set to 'GO' then the script will merely display a"
  echo "  list of all grep results for the original term in all code files.  Look "
  echo "  before you leap!"
  exit
fi

if [ ! ${3} ] ; then
  N=0
elif [ ${3} == "GO" ] ; then
  echo "#!/bin/bash" > meaninglessChangeScript.sh
fi
ALL_EXT="cpp h tpp cu cuh cui"
for EXT in ${ALL_EXT} ; do
  for LEVEL in 0 1 2 3 ; do
    KK="NOT_A_FILE "
    if [ ${LEVEL} -eq 0 ] ; then
      KK="${KK} `ls *.${EXT}`"
    elif [ ${LEVEL} -eq 1 ] ; then
      KK="${KK} `ls */*.${EXT}`"
    elif [ ${LEVEL} -eq 2 ] ; then
      KK="${KK} `ls */*/*.${EXT}`"
    elif [ ${LEVEL} -eq 3 ] ; then
      KK="${KK} `ls */*/*/*.${EXT}`"
    fi
    if [ "${KK}" == "NOT_A_FILE " ] ; then
      echo "No ${BASE}... ls ${BASE}.${EXT}"
    else
      for FI in ${KK} ; do
        if [ ! -e ${FI} ] ; then
          continue
        fi
        TT="NOT_A_RESULT "
        TT="${TT} `grep ${1} ${FI}`"
        if [ ${#TT} -le 15 ] ; then
          N=0
        else
          if [ ! ${3} ] ; then
            echo "Instances in ${FI}:"
            grep "${1}" ${FI}
            echo
          elif [ ${3} == "GO" ] ; then
            echo "sed -i 's/${1}/${2}/g' ${FI}" >> meaninglessChangeScript.sh
          fi
        fi
      done
    fi
  done
done
if [ ! ${3} ] ; then
  N=0
elif [ ${3} == "GO" ] ; then
  chmod +x meaninglessChangeScript.sh
  ./meaninglessChangeScript.sh
  rm meaninglessChangeScript.sh
fi
