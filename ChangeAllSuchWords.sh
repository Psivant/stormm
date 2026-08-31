#!/usr/bin/env bash

# zsh globbing and [ ] test syntax differ from bash; re-invoke under bash when needed.
if [ -n "${ZSH_VERSION:-}" ]; then
  exec env bash "$0" "$@"
fi

# BSD sed (macOS) requires an argument after -i; GNU sed (Linux) does not.
if [ "$(uname -s)" = "Darwin" ]; then
  run_sed_inplace() {
    sed -i '' "$@"
  }
else
  run_sed_inplace() {
    sed -i "$@"
  }
fi

if [ -z "${1:-}" ] || [ -z "${2:-}" ] ; then
  echo "Usage: ./ChangeAllSuchWords <original term> <replacement term> <authorization>"
  echo
  echo "  If the authorization is not set to 'GO' then the script will merely display a"
  echo "  list of all grep results for the original term in all code files.  Look "
  echo "  before you leap!"
  exit
fi

shopt -s nullglob

if [ -z "${3:-}" ] ; then
  N=0
elif [ "${3}" = "GO" ] ; then
  echo "#!/usr/bin/env bash" > meaninglessChangeScript.sh
  if [ "$(uname -s)" = "Darwin" ]; then
    echo "# macOS (BSD sed): in-place edits use sed -i ''" >> meaninglessChangeScript.sh
  fi
fi
ALL_EXT="cpp h tpp cu cuh cui"
for EXT in ${ALL_EXT} ; do
  for LEVEL in 0 1 2 3 ; do
    KK=( )
    if [ ${LEVEL} -eq 0 ] ; then
      KK=( *.${EXT} )
    elif [ ${LEVEL} -eq 1 ] ; then
      KK=( */*.${EXT} )
    elif [ ${LEVEL} -eq 2 ] ; then
      KK=( */*/*.${EXT} )
    elif [ ${LEVEL} -eq 3 ] ; then
      KK=( */*/*/*.${EXT} )
    fi
    if [ ${#KK[@]} -eq 0 ] ; then
      continue
    fi
    for FI in "${KK[@]}" ; do
      if [ ! -e "${FI}" ] ; then
        continue
      fi
      if ! grep -q "${1}" "${FI}" ; then
        continue
      fi
      if [ -z "${3:-}" ] ; then
        echo "Instances in ${FI}:"
        grep "${1}" "${FI}"
        echo
      elif [ "${3}" = "GO" ] ; then
        if [ "$(uname -s)" = "Darwin" ]; then
          echo "sed -i '' 's/${1}/${2}/g' \"${FI}\"" >> meaninglessChangeScript.sh
        else
          echo "sed -i 's/${1}/${2}/g' \"${FI}\"" >> meaninglessChangeScript.sh
        fi
      fi
    done
  done
done
if [ -z "${3:-}" ] ; then
  N=0
elif [ "${3}" = "GO" ] ; then
  chmod +x meaninglessChangeScript.sh
  ./meaninglessChangeScript.sh
  rm meaninglessChangeScript.sh
fi
