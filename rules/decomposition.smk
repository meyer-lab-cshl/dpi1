rule decompositions_short:
    input:
        tc="{dir}/outPooled/tc.bed.gz",
        fwd="{dir}/outPooled/ctssTotalCounts.fwd.bw",
        rev="{dir}/outPooled/ctssTotalCounts.rev.bw",
        genome=config["genome"]
    output:
        tc="{dir}/outPooled/tc.short.bed.gz"
    params:
        length=config["decomposition"]["length"],
        count=config["decomposition"]["count"],
        ratio=config["decomposition"]["noise_subtraction_ratio"],
    shell:
        """
        gunzip -c {input.tc} | \
        awk '{{if((($3 - $2) <= {params.length}) || ($5 <= {params.count})){{print}}}}' | \
        grep -- '+$' | \
        scripts/bed2peakBed9_with_boundaryTrimming.sh \
            -r {params.ratio} \
            -b {wildcards.dir}/outPooled/ctssTotalCounts.fwd.bw \
            -c '255,0,0' \
        > {output.tc}.tmp
        gunzip -c {input.tc} | \
        awk '{{if((($3 - $2) <= {params.length}) || ($5 <= {params.count})){{print}}}}' | \
        grep -- '-$' | \
        scripts/bed2peakBed9_with_boundaryTrimming.sh \
            -r {params.ratio} \
            -b {wildcards.dir}/outPooled/ctssTotalCounts.rev.bw \
            -c '0,0,255' \
        >> {output.tc}.tmp
        sort -k1,1 -k2,2n {output.tc}.tmp | \
        gzip -c > {output.tc}
        rm -f {output.tc}.tmp
        """

rule decompositions_long:
    input:
        tc="{dir}/outPooled/tc.bed.gz",
        genome=config["genome"]
    output:
        #dynamic("{dir}/outPooled/tc.long/tc.{suffix}"),
        alltc="{dir}/outPooled/tc.long.bed.gz"
    params:
        length=config["decomposition"]["length"],
        count=config["decomposition"]["count"],
        ratio=config["decomposition"]["noise_subtraction_ratio"],
    shell:
        """
        # Read tag clusters; filter by > length and count; zip to file
        gunzip -c {input.tc} | \
        awk '{{if((($3 - $2) > {params.length}) || ($5 > {params.count})){{print}}}}' | \
        gzip -c > {output.alltc}
        # unzip and split file into small chunks; named by `suffix_length` letters
        # to form the suffix of the file name; `lines` per file; prefix `tc`
        gunzip -c {output.alltc} | \
        gsplit \
            --suffix-length=5 \
            --lines=1000  - {wildcards.dir}/outPooled/tc.long/
        """
