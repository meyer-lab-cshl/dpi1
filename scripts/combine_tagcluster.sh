#!/usr/bin/env bash

function usage()
{
  cat <<EOF
Usage: $0
        [ -f tc.fwd.bed ]
        [ -r tc.rev.bed ]
        [ -o tc.bed.gz ]
EOF
  exit 1;
}

while getopts m:t:c:d:s:o: opt; do
  case ${opt} in
  f) fwd=${OPTARG};;
  r) rev=${OPTARG};;
  o) output=${OPTARG};;
  *) usage;;
  esac
done

if [ "${fwd}" = "" ]; then usage; fi
if [ "${rev}" = "" ]; then usage; fi
if [ "${output}" = "" ]; then usage; fi

cat ${fwd} ${rev} | \
sed -e 's/[:|,|..]/\t/g' | \
awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,$1":"$2".."$3","$4,$5,$4
}' | \
sort -k1,1 -k2,2n | \
gzip -c > ${output}
