rule split_short:
    input:
        tc="{dir}/outPooled/tc.bed.gz",
        bw="{dir}/outPooled/ctssTotalCounts.{strand}.bw",
        genome=config["genome"]
    output:
        bw="{dir}/outPooled/tc.short.{strand}.bed"
    conda:
        "../envs/split.yaml"
    params:
        length=config["split"]["length"],
        ct=config["split"]["count"],
        ratio=config["split"]["noise_subtraction_ratio"],
        color=lambda wildcards: '255,0,0' if wildcards.strand == "fwd" else '0,0,255',
        sign=lambda wildcards: '+' if wildcards.strand == "fwd" else '-'
    shell:
        """
        scripts/split_short.sh \
            -b {input.bw} \
            -t {input.tc} \
            -g {input.genome} \
            -c {params.ct} \
            -l {params.length} \
            -n {params.ratio} \
            -C {params.color} \
            -s {params.sign} \
            -o {output.bw}
        """
rule merge_short:
    input:
        fwd="{dir}/outPooled/tc.short.fwd.bed",
        rev="{dir}/outPooled/tc.short.rev.bed",
    output:
        tc="{dir}/outPooled/tc.short.bed.gz"
    shell:
        """
        cat {input.fwd} {input.rev} | \
        sort -k1,1 -k2,2n | \
        gzip -c > {output.tc}
        """

rule split_long:
    input:
        tc="{dir}/outPooled/tc.bed.gz",
    output:
        tc="{dir}/outPooled/tc.long.bed.gz"
    params:
        length=config["split"]["length"],
        ct=config["split"]["count"],
    shell:
        """
        scripts/split_long.sh \
            -t {input.tc} \
            -l {params.length} \
            -c {params.ct} \
            -o {output.tc}
        """
