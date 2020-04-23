rule tc_long_stranded:
    input:
        tc="{dir}/outPooled/tc.long.spi.bed.gz",
        bw="{dir}/outPooled/ctssTotalCounts.{strand}.bw",
    params:
        color=lambda wildcards: '255,0,0' if wildcards.strand == "fwd" else '0,0,255',
        sign=lambda wildcards: '+' if wildcards.strand == "fwd" else '-'
    output:
        bed="{dir}/outPooled/tc.long.spi.{strand}.bed.gz"
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

