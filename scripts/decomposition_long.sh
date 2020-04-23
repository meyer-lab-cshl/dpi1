#!/usr/bin/env bash

function usage()
{
  cat <<EOF
Usage: $0
        [ -t tagcluster_file ]
        [ -c threshold_counts ]
        [ -l threshold_length ]
        [ -o output_file ]"

EOF
  exit 1;
}

while getopts t:c:l:o: opt; do
  case ${opt} in
  t) tc=${OPTARG};;
  c) counts=${OPTARG};;
  l) length=${OPTARG};;
  o) output=${OPTARG};;
  *) usage;;
  esac
done

if [ "${tc}" = "" ]; then usage; fi
if [ "${counts}" = "" ]; then usage; fi
if [ "${length}" = "" ]; then usage; fi
if [ "${output}" = "" ]; then usage; fi

gunzip -c ${tc} | \
awk -v length=${length} -v count=${count} '{
    if((($3 - $2) > length) || ($5 > count)){print}
}' | \
gzip -c > ${output}
