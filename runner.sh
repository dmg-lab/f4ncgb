#!/usr/bin/env bash

add_to_d=0
executable="~/sai/f4ncgb/build-release/f4ncgb"

POSITIONAL_ARGS=()

enable_outandproof=false
outandproof=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --add-to-d)
      add_to_d=$2
      shift
      shift
      ;;
    --f4ncgb)
      echo "[f4ncgb-runner] Using executable $2 instead of $executable"
      executable=$2
      shift
      shift
      ;;
    --with-out-and-proof)
      enable_outandproof=true
      shift
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}"

if [[ "$1" =~ ^.*-([0-9]+)\.ms$ ]]; then
  dint=${BASH_REMATCH[1]}
  echo "[f4ncgb-runner] Extracted -d$dint from argument $param, adding $add_to_d to it"
  ((dint+=$add_to_d))
else
  echo "[freebg-runner] Regex did not match!"
  exit 1
fi

if [ "$enable_outandproof" = true ]; then
  outpath="/tmp/f4ncgb-out"
  mkdir -p $outpath
  outandproof="-o $outpath/$(basename $1).out -p $outpath/$(basename $1).proof"
  echo "[f4ncgb-runner] Using additional parameters $outandproof"
fi

"$executable" -d$dint $outandproof $@
