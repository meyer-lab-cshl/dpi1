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

while getopts b:t:g:c:l:n:s:o: opt; do
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

gunzip -c ${tc} | \
awk -v length=${length} -v count=${count} '{
    if((($3 - $2) <= length) || ($5 <= count)){print}
}' | \
grep -- $signstr | \
awk -v bw=${bw} -v color=${color} -v noise_subtraction_ratio=${noise} \
'BEGIN{OFS="\t"}{

  chrom = $1; start = $2; stop = $3; name = $4; score = $5; strand = $6;

  ### for representative position ###
  rep_start = 0; rep_stop = 0; rep_max = 0;
  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout",
  chrom, start, stop, bw )
  while ((command | getline ) > 0)
  {
    #total = total + ($4 * ($3 - $2))
    if ($4 > rep_max){
      rep_start = $2
      rep_stop  = $3
      rep_max   = $4
    }
  }
  close(command)

  ### for boundary trimming by noise subtraction ###
  boundary_start = 0; boundary_stop = 0 ;
  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout | sort -k1,1 -k2,2n ",
  chrom, start, stop, bw )
  while ((command | getline ) > 0)
  {
    if ( ( $4 - (rep_max * noise_subtraction_ratio) ) > 0  ){
      boundary_start = $2
      break
    }
  }
  close(command)

  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout | sort -k1,1 -k2,2nr ",
  chrom, start, stop, bw )
  while ((command | getline ) > 0)
  {
    if ( ( $4 - (rep_max * noise_subtraction_ratio) ) > 0  ){
      boundary_stop = $3
      break
    }
  }
  close(command)

  ### get total counts on the new boundary ###
  total = 0
  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout",
  chrom, boundary_start, boundary_stop, bw )
  while ((command | getline ) > 0)
  {
    total = total + ($4 * ($3 - $2))
  }
  close(command)

  if (total > 0)
  {
    # added this line at 13th Jan, 2012
    name=chrom":"boundary_start".."boundary_stop","strand

    print chrom, boundary_start, boundary_stop, name, total, strand,
          rep_start, rep_stop, color
  }
}' | \
> ${output}
