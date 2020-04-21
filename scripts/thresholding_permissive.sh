#!/usr/bin/env bash

function usage()
{
  cat <<EOF
Usage: $0
        [ -f ctssMaxCounts_fwd ]
        [ -r ctssMaxCounts_rev ]
        [ -c threshold_counts ]
        [ -g smoothed_tag_cluster_file ]
        [ -o output_file ]"

EOF
  exit 1;
}

while getopts f:r:c:t:o:g: opt; do
  case ${opt} in
  f) counts_fwd=${OPTARG};;
  r) counts_rev=${OPTARG};;
  c) counts=${OPTARG};;
  o) output=${OPTARG};;
  g) tc=${OPTARG};;
  *) usage;;
  esac
done

if [ "${counts_fwd}" = "" ]; then usage; fi
if [ "${counts_rev}" = "" ]; then usage; fi
if [ "${counts}" = "" ]; then usage; fi
if [ "${output}" = "" ]; then usage; fi
if [ "${tc}" = "" ]; then usage; fi

bigWigToBedGraph ${counts_fwd} /dev/stdout | \
awk -v counts=${counts} 'BEGIN{OFS="\t"}{
    if($4 >= counts){print $1,$2,$3,".",$4,"+"}
}' > ${output}.tmp

bigWigToBedGraph ${counts_rev} /dev/stdout | \
awk -v counts=${counts} 'BEGIN{OFS="\t"}{
 if($4 >= counts){print $1,$2,$3,".",$4,"-"}
}' >> ${output}.tmp

gunzip -c ${tc} | \
intersectBed \
    -s \
    -wa \
    -u \
    -a stdin \
    -b ${output}.tmp | \
gzip -c > ${output}
rm ${output}.tmp
