#!/usr/bin/env bash

function usage()
{
  cat <<EOF
  Usage: $0
          [ -f ctssMaxCounts_fwd ]
          [ -r ctssMaxCounts_rev ]
          [ -F ctssMaxTPM_fwd ]
          [ -R ctssMaxTPM_rev ]
          [ -c threshold_counts ]
          [ -t threshold_tpm ]
          [ -g smoothed_tag_cluster_file ]
          [ -o output_file ]"

EOF
  exit 1;
}

while getopts f:r:F:R:c:t:o:g: opt; do
  case ${opt} in
  f) counts_fwd=${OPTARG};;
  r) counts_rev=${OPTARG};;
  F) tpm_fwd=${OPTARG};;
  R) tpm_rev=${OPTARG};;
  c) counts=${OPTARG};;
  t) tpm=${OPTARG};;
  o) output=${OPTARG};;
  g) tc=${OPTARG};;
  *) usage;;
  esac
done

if [ "${counts_fwd}" = "" ]; then usage; fi
if [ "${counts_rev}" = "" ]; then usage; fi
if [ "${tpm_fwd}" = "" ]; then usage; fi
if [ "${tpm_rev}" = "" ]; then usage; fi
if [ "${counts}" = "" ]; then usage; fi
if [ "${tpm}" = "" ]; then usage; fi
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
gzip -c > ${output}.tmp.gz

bigWigToBedGraph ${tpm_fwd} /dev/stdout | \
awk -v counts=${tpm} 'BEGIN{OFS="\t"}{
    if($4 >= counts){print $1,$2,$3,".",$4,"+"}
}' > ${output}.tmp2

bigWigToBedGraph ${tpm_rev} /dev/stdout | \
awk -v counts=${tpm} 'BEGIN{OFS="\t"}{
 if($4 >= counts){print $1,$2,$3,".",$4,"-"}
}' >> ${output}.tmp2

intersectBed \
    -s \
    -wa \
    -u \
    -a ${output}.tmp.gz \
    -b ${output}.tmp2 | \
gzip -c > ${output}
rm ${output}.tmp.gz ${output}.tmp ${output}.tmp2
