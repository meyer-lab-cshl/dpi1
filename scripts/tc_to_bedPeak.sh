#!/usr/bin/env bash

function usage()
{
  cat <<EOF
Usage: $0
        [ -l tc_long ]
        [ -b ctssTotalCounts.fwd|rev.bw ]
        [ -c color ]
        [ -s sign ]
        [ -o output ]

EOF
  exit 1;
}

while getopts l:b:c:s:o: opt; do
  case ${opt} in
  l) long=${OPTARG};;
  b) bw=${OPTARG};;
  c) color=${OPTARG};;
  s) sign=${OPTARG};;
  o) output=${OPTARG};;
  *) usage;;
  esac
done

if [ "${long}" = "" ]; then usage; fi
if [ "${bw}" = "" ]; then usage; fi
if [ "${color}" = "" ]; then usage; fi
if [ "${sign}" = "" ]; then usage; fi
if [ "${output}" = "" ]; then usage; fi

gunzip -c ${long} | \
grep -v ^# | \
awk -v sign=${sign} '{if($6 == sign){print}}' | \
awk 'BEGIN{OFS="\t"}{
    print $1,$2,$3, $1":"$2".."$3","$4, 1000,$6
}' | \
awk -v bw=${bw} -v color=${color} 'BEGIN{OFS="\t"}{
  chrom = $1
  start = $2
  stop  = $3
  name  = $4
  score  = $5
  strand = $6

  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout",
  chrom, start, stop, bw )

  rep_start = 0
  rep_stop  = 0
  rep_max   = 0
  total = 0
  while ((command | getline ) > 0)
  {
    total = total + ($4 * ($3 - $2))
    if ($4 > rep_max){
      rep_start = $2
      rep_stop  = $3
      rep_max   = $4
    }
  }
  close(command)

  score = total
  if (score > 0)
  {
    print chrom, start, stop, name, score, strand,
          rep_start, rep_stop, color
  }
}' |  gzip -c > ${output}
