#!/usr/bin/env bash

function usage()
{
  cat <<EOF
Usage: $0
        [ -s tc_short ]
        [ -f tc.long.spi.fwd.bed.gz ]
        [ -r tc.long.spi.rev.bed.gz ]
        [ -o output ]

EOF
  exit 1;
}

while getopts s:f:r:o: opt; do
  case ${opt} in
  s) short=${OPTARG};;
  f) fwd=${OPTARG};;
  r) rev=${OPTARG};;
  o) output=${OPTARG};;
  *) usage;;
  esac
done

if [ "${short}" = "" ]; then usage; fi
if [ "${fwd}" = "" ]; then usage; fi
if [ "${rev}" = "" ]; then usage; fi
if [ "${output}" = "" ]; then usage; fi

# combined strands short
gunzip -c ${fwd} ${rev} ${short} | \
sort -k1,1 -k2,2n | \
awk 'BEGIN{OFS="\t"}{
    print $0,1,$3-$2",",0","
}' | \
gzip -c > ${output}
