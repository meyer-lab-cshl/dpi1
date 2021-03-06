#!/usr/bin/env bash


### prep
function usage()
{
  cat <<EOF
usage: $0
          [ -i INFILE.ctss.bed|INFILE.ctss.bed.gz;
                full path to ctss files in bed format, optionally gzipped]
          [ -o OUT_PREFIX;
                full path prefix to output files ]
          [ -g GENOME;
                full path the genome file as obtained in bedtools/genomes ]
EOF
  exit 1;
}

while getopts i:o:g: opt
do
  case ${opt} in
  i) infile=${OPTARG};;
  o) outprefix=${OPTARG};;
  g) genome=${OPTARG};;
  *) usage;;
  esac
done

if [ "${infile}" = "" ]; then usage; fi
if [ "${outprefix}" = "" ]; then usage; fi
if [ "${genome}" = "" ]; then usage; fi

echo "infile: $infile"
echo "outprefix: $outprefix"
echo "genome: $genome"


fwd=${outprefix}.fwd.bw
rev=${outprefix}.rev.bw
tmpfile=${outprefix}.tmp.bg
tmpfile_g=${outprefix}.tmp_g

sort ${genome} > ${tmpfile_g}

### fwd
if file --mime-type "$infile" | grep -q gzip; then
    echo "input file is gz"
    gunzip -c ${infile} \
    | grep ^chr \
    | grep +$ \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5}' \
    | sort -k1,1 -k2,2n \
    > ${tmpfile}

    bedGraphToBigWig ${tmpfile} ${tmpfile_g} ${fwd}

    ### rev
    gunzip -c ${infile} \
    | grep ^chr \
    | grep -v +$ \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5}' \
    | sort -k1,1 -k2,2n \
    > ${tmpfile}
    bedGraphToBigWig ${tmpfile} ${tmpfile_g} ${rev}
else
    echo "input file is unzipped"
    cat ${infile} \
    | grep ^chr \
    | grep +$ \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5}' \
    | sort -k1,1 -k2,2n \
    > ${tmpfile}

    bedGraphToBigWig ${tmpfile} ${tmpfile_g} ${fwd}

    ### rev
    cat ${infile} \
    | grep ^chr \
    | grep -v +$ \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5}' \
    | sort -k1,1 -k2,2n \
    > ${tmpfile}
    bedGraphToBigWig ${tmpfile} ${tmpfile_g} ${rev}
fi

rm -f ${tmpfile}
rm -f ${tmpfile_g}

