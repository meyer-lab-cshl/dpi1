rule decompositions_short_stranded:
    input:
        tc="{dir}/outPooled/tc.bed.gz",
        bw="{dir}/outPooled/ctssTotalCounts.{strand}.bw",
        genome=config["genome"]
    output:
        bw="{dir}/outPooled/tc.short.{strand}.bed.gz",
    params:
        length=config["decomposition"]["length"],
        count=config["decomposition"]["count"],
        ratio=config["decomposition"]["noise_subtraction_ratio"],
        color=lambda wildcards: '255,0,0' if strand == "fwd" else '0,0,255',
        sign=lambda wildcards: '+' if strand == "fwd" else '-'
    shell:
        """
        scripts/decomposition_short.sh \
            -b {input.bw} \
            -t {input.tc} \
            -g {input.genome} \
            -c {params.count} \
            -l {params.length} \
            -n {params.ratio} \
            -C {params.color} \
            -s {params.sign} \
            -o {output.bw}
        """
rule decomposition_short_merge:
    input:
        fwd="{dir}/outPooled/tc.short.fwd.bed.gz",
        rev="{dir}/outPooled/tc.short.rev.bed.gz",
    output:
        tc="{dir}/outPooled/tc.short.bed.gz"
    shell:
        """
        gunzip -c {input.fwd} {input.rev} | \
        sort -k1,1 -k2,2n | \
        gzip -c > {output.tc}
        """

rule decompositions_long:
    input:
        tc="{dir}/outPooled/tc.bed.gz",
    output:
        tc="{dir}/outPooled/tc.long.bed.gz"
    params:
        length=config["decomposition"]["length"],
        count=config["decomposition"]["count"],
    shell:
        """
        scripts/decompositions_long.sh \
            -t {input.tc} \
            -l {params.length} \
            -c {params.count} \
            -o {output.tc}
        """
