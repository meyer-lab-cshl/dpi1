rule smoothing:
    input:
        tc_long="{dir}/outPooled/tc.long.bed.gz",
        fwd="{dir}/outPooled/ctssTotalCounts.fwd.bw",
        rev="{dir}/outPooled/ctssTotalCounts.rev.bw",
    output:
        spi="{dir}/outPooled/tc.long.spi.bed.gz"
    conda:
        "../envs/smooth.yaml"
    params:
        length=config["split"]["length"],
        ratio=config["split"]["noise_subtraction_ratio"],
        window=config["smoothing"]["gaussian_window_size_half"],
    shell:
        """
        Rscript scripts/decompose_peakclustering_hvm.R \
            --tagclusters {input.tc_long} \
            --ctssprefix {wildcards.dir}/outPooled/ctssTotalCounts \
            --window {params.window} \
            --length {params.length} \
            --ratio {params.ratio} \
            --outfile {wildcards.dir}/outPooled/tc.long.spi.bed \
            --analysis spi
        gzip {wildcards.dir}/outPooled/tc.long.spi.bed
        """

rule tc_long_stranded:
    input:
        tc="{dir}/outPooled/tc.long.spi.bed.gz",
        bw="{dir}/outPooled/ctssTotalCounts.{strand}.bw",
    output:
        bed=temp("{dir}/outPooled/tc.long.spi.{strand}.bed.gz")
    conda:
        "../envs/merge_simple.yaml"
    params:
        color=lambda wildcards: '255,0,0' if wildcards.strand == "fwd" else '0,0,255',
        sign=lambda wildcards: '+' if wildcards.strand == "fwd" else '-'
    shell:
        """
        scripts/tc_to_bedPeak.sh \
            -l {input.tc} \
            -b {input.bw} \
            -c {params.color} \
            -s {params.sign} \
            -o {output.bed}
        """

rule tc_merge:
    input:
        fwd="{dir}/outPooled/tc.long.spi.fwd.bed.gz",
        rev="{dir}/outPooled/tc.long.spi.rev.bed.gz",
        short="{dir}/outPooled/tc.short.bed.gz",
    output:
        merged="{dir}/outPooled/tc.spi.merged.bed.gz"
    shell:
        """
        scripts/merge_simple.sh \
            -s {input.short} \
            -f {input.fwd} \
            -r {input.rev} \
            -o {output.merged}
        """

