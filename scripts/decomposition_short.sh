#!/usr/bin/env bash

function usage()
{
  cat <<EOF
Usage: $0
        [ -b ctssTotalCounts_fwd|rev ]
        [ -t tagcluster_file ]
        [ -g genome file ]
        [ -c threshold_counts ]
        [ -l threshold_length ]
        [ -n noise_subtraction_ratio ]
        [ -C color ]
        [ -s sign ]
        [ -o output_file ]"

EOF
  exit 1;
}

while getopts b:t:g:c:l:n:C:s:o: opt; do
  case ${opt} in
  b) bw=${OPTARG};;
  t) tc=${OPTARG};;
  g) genome=${OPTARG};;
  c) counts=${OPTARG};;
  l) length=${OPTARG};;
  n) noise=${OPTARG};;
  C) color=${OPTARG};;
  s) sign=${OPTARG};;
  o) output=${OPTARG};;
  *) usage;;
  esac
done

if [ "${bw}" = "" ]; then usage; fi
if [ "${tc}" = "" ]; then usage; fi
if [ "${genome}" = "" ]; then usage; fi
if [ "${counts}" = "" ]; then usage; fi
if [ "${length}" = "" ]; then usage; fi
if [ "${noise}" = "" ]; then usage; fi
if [ "${color}" = "" ]; then usage; fi
if [ "${sign}" = "" ]; then usage; fi
if [ "${output}" = "" ]; then usage; fi
signstr="${sign}$"

echo "b (ctssTotalCounts): ${b}"
echo "t (tagcluster): ${tc}"
echo "g (genome): ${genome}"
echo "c (threshold_counts): ${counts}"
echo "l (threshold_length): ${length}"
echo "n (noise_subtraction_ratio): ${noise}"
echo "C (color): ${color}"
echo "s (sign): ${sign}"
echo "o (output): ${output}"


gunzip -c ${tc} | \
awk -v clength=${length} -v ccount=${counts} '{
    if((($3 - $2) <= clength) || ($5 <= ccount)){print}
}' | \
grep -- $signstr | \
scripts/bed2peakBed9_with_boundaryTrimming.sh \
    -r ${noise} \
    -b ${bw} \
    -c ${color} > ${output}
