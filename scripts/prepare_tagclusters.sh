#!/usr/bin/env bash

function usage()
{
  cat <<EOF
Usage: $0
        [ -m ctssMaxCounts ]
        [ -t ctssTotalCounts ]
        [ -c threshold_counts ]
        [ -d distance ]
        [ -s sign ]
        [ -o output_file ]

EOF
  exit 1;
}

while getopts m:t:c:d:s:o: opt; do
  case ${opt} in
  m) max=${OPTARG};;
  t) total=${OPTARG};;
  c) thr=${OPTARG};;
  d) dist=${OPTARG};;
  s) sign=${OPTARG};;
  o) output=${OPTARG};;
  *) usage;;
  esac
done

if [ "${max}" = "" ]; then usage; fi
if [ "${total}" = "" ]; then usage; fi
if [ "${thr}" = "" ]; then usage; fi
if [ "${dist}" = "" ]; then usage; fi
if [ "${sign}" = "" ]; then usage; fi
if [ "${output}" = "" ]; then usage; fi

bigWigToBedGraph ${max} /dev/stdout | \
awk -v thr=${thr} 'BEGIN{OFS="\t"}{
    if($4 >=thr){print $1,$2,$3,$4}
}' | \
sort -k1,1 -k2,2n | \
mergeBed -d $dist | \
awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,$1":"$2".."$3
}' | \
awk -v sign=${sign} 'BEGIN{OFS="\t"} {
    print $1,$2,$3,$4","sign,1000,sign
}' | \
grep -v '_' | \
bigWigAverageOverBed ${total} /dev/stdin /dev/stdout | \
awk '{printf "%s\t%i\n",$1,$4}' \
> ${output}

