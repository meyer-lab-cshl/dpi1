rule permissive:
    input:
        fwd_counts="{dir}/outPooled/ctssMaxCounts.fwd.bw",
        rev_counts="{dir}/outPooled/ctssMaxCounts.rev.bw",
        spi="{dir}/outPooled/tc.spi.merged.bed.gz"
    output:
        bedgraph="{dir}/outPooled/ctssMaxCounts{counts}.all.bg",
        spi="{dir}/outPooled/tc.spi_merged.ctssMaxCounts{counts}.bed.gz"
    shell:
        """
        bigWigToBedGraph {input.fwd_counts} /dev/stdout | \
        awk -v cutoff={wildcards.counts} \
            'BEGIN{{OFS="\\t"}}{{if($4 >= cutoff ){{print $1,$2,$3,".",$4,"+"}}}}' \
        > {output.spi}.tmp
        bigWigToBedGraph {input.rev_counts} /dev/stdout | \
        awk -v cutoff={wildcards.counts} \
            'BEGIN{{OFS="\t"}}{{if($4 >= cutoff ){{print $1,$2,$3,".",$4,"-"}}}}' \
        >> {output.spi}.tmp

        gunzip -c {input.spi} | \
        intersectBed \
            -s \
            -wa \
            -u \
            -a stdin \
            -b {output.spi}.tmp | \
        gzip -c > {output.spi}
        rm {output.spi}.tmp
        """

rule robust:
    input:
        fwd_tpm="{dir}/outPooled/ctssMaxTpm.fwd.bw",
        rev_tpm="{dir}/outPooled/ctssMaxTpm.rev.bw",
        spi="{dir}/outPooled/tc.spi_merged.bed.gz",
    output:
        spi="{dir}/outPooled/tc.spi_merged.ctssMaxCounts{counts}_ctssMaxTpm{tpm}.bed.gz"
    shell:
        """
        bigWigToBedGraph {input.fwd_counts} /dev/stdout | \
        awk -v cutoff={wildcards.counts} \
            'BEGIN{{OFS="\\t"}}{{if($4 >= cutoff ){{print $1,$2,$3,".",$4,"+"}}}}' \
        > {output.spi}.tmp
        bigWigToBedGraph {input.rev_counts} /dev/stdout | \
        awk -v cutoff={wildcards.counts} \
            'BEGIN{{OFS="\\t"}}{{if($4 >= cutoff ){print $1,$2,$3,".",$4,"-"}}}}' \
        >> {output.spi}.tmp

        gunzip -c {input.spi} | \
        intersectBed \
            -s \
            -wa \
            -u \
            -a stdin \
            -b {output.spi}.tmp | \
        gzip -c > {output.spi}.tmp.gz

        bigWigToBedGraph {input.fwd_tpm} /dev/stdout | \
        awk -v cutoff={wildcards.tpm} \
            'BEGIN{{OFS="\\t"}}{{if($4 >= cutoff ){{print $1,$2,$3,".",$4,"+"}}}}' \
        > {output.spi}.tmp2
        bigWigToBedGraph {input.rev_tpm} /dev/stdout | \
        awk -v cutoff={wildcards.tpm} \
            'BEGIN{{OFS="\\t"}}{{if($4 >= cutoff ){{print $1,$2,$3,".",$4,"-"}}}}' \
        >> {output.spi}.tmp2

        intersectBed \
            -s \
            -wa \
            -u \
            -a {output.spi}.tmp.gz \
            -b {output.spi}.tmp2 | \
        gzip -c > {output.spi}
        rm {output.spi}.tmp.gz {output.spi}.tmp {output.spi}.tmp2
        """
