#!/usr/bin/env bash

add_to_d=0
executable="~/sai/kommunopp/build-release/freegb"

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    --add-to-d)
      add_to_d=$2
      shift
      shift
      ;;
    --freegb)
      echo "[freegb-runner] Using executable $2 instead of $executable"
      executable=$2
      shift
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
  echo "[freegb-runner] Extracted -d$dint from argument $param, adding $add_to_d to it"
  ((dint+=$add_to_d))
else
  echo "[freebg-runner] Regex did not match!"
  exit 1
fi

"$executable" -d$dint $@
