rule merge:
    input:
        spi="{dir}/outPooled/tc.long.spi.bed.gz",
        short="{dir}/outPooled/tc.short.bed.gz",
        fwd="{dir}/outPooled/ctssTotalCounts.fwd.bw",
        rev="{dir}/outPooled/ctssTotalCounts.rev.bw",
    output:
        spi="{dir}/outPooled/tc.spi.merged.bed.gz"
    shell:
        """
        # forward strand long
        gunzip -c {input.spi} | \
        grep -v ^# | \
        awk '{{if($6=="+"){{print}}}}' | \
        awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3, $1":"$2".."$3","$4, 1000,$6}}' | \
        scripts/bed2peakBed9.sh \
            -b {input.fwd} \
            -c "255,0,0" \
        > {output.spi}.tmp
        # reverse strand long
        gunzip -c {input.spi} | \
        grep -v ^# | \
        awk '{{if($6=="-"){{print}}}}' | \
        awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3, $1":"$2".."$3","$4, 1000,$6}}' | \
        scripts/bed2peakBed9.sh \
            -b {input.rev} \
            -c "255,0,0" \
        >> {output.spi}.tmp
        # combined strands short
        gunzip -c {input.short} >> {output.spi}.tmp
        sort -k1,1 -k2,2n {output.spi}.tmp | \
        awk 'BEGIN{{OFS="\\t"}}{{print $0,1,$3-$2",",0","}}' | \
        gzip -c > {output.spi}
        rm -f {output.spi}.tmp
        """
