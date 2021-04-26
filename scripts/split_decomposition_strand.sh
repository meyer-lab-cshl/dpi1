#!/usr/bin/env bash

function usage()
{
  cat <<EOF
Usage: $0
        [ -b bed_file ]
        [ -c counts_file ]
        [ -s sign ]
        [ -C color ]
        [ -o output_file ]"

EOF
  exit 1;
}

while getopts b:c:s:C:o: opt; do
  case ${opt} in
  b) bed=${OPTARG};;
  c) counts=${OPTARG};;
  s) sign=${OPTARG};;
  C) color=${OPTARG};;
  o) output=${OPTARG};;
  *) usage;;
  esac
done

if [ "${bed}" = "" ]; then usage; fi
if [ "${counts}" = "" ]; then usage; fi
if [ "${sign}" = "" ]; then usage; fi
if [ "${color}" = "" ]; then usage; fi
if [ "${output}" = "" ]; then usage; fi

gunzip -c ${bed} |
awk -v strand=${sign} 'BEGIN{OFS="\t"} {if($6 == strand) \
    {print $1,$2,$3,$4,$5,$6}}' |
mergeBed -i stdin | \
awk -v strand=${sign} 'BEGIN{OFS="\t"} \
  {print $1,$2,$3, $1":"$2".."$3","strand, 1000,strand }' | \
scripts/bed2peakBed9.sh \
  -b ${counts} \
  -c ${color} \
> ${output}
