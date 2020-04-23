rule permissive:
    input:
        fwd_counts="{dir}/outPooled/ctssMaxCounts.fwd.bw",
        rev_counts="{dir}/outPooled/ctssMaxCounts.rev.bw",
        spi="{dir}/outPooled/tc.spi.merged.bed.gz"
    output:
        spi="{dir}/outPooled/tc.spi.merged.ctssMaxCounts{counts}.bed.gz"
    conda:
        "../envs/thresholding_simple.yaml"
    wildcard_constraints:
        counts="\d+"
    shell:
        """
        scripts/thresholding_permissive.sh \
            -f {input.fwd_counts} \
            -r {input.rev_counts} \
            -c {wildcards.counts} \
            -o {output.spi} \
            -g {input.spi}
        """

rule robust:
    input:
        fwd_counts="{dir}/outPooled/ctssMaxCounts.fwd.bw",
        rev_counts="{dir}/outPooled/ctssMaxCounts.rev.bw",
        fwd_tpm="{dir}/outPooled/ctssMaxTpm.fwd.bw",
        rev_tpm="{dir}/outPooled/ctssMaxTpm.rev.bw",
        spi="{dir}/outPooled/tc.spi.merged.bed.gz",
    output:
        spi="{dir}/outPooled/tc.spi.merged.ctssMaxCounts{counts}.ctssMaxTpm{tpm}.bed.gz"
    conda:
        "../envs/thresholding_simple.yaml"
    wildcard_constraints:
        counts="\d+",
        tpm="\d+"
    shell:
        """
        scripts/thresholding_robust.sh \
            -f {input.fwd_counts} \
            -r {input.rev_counts} \
            -F {input.fwd_tpm} \
            -R {input.rev_tpm} \
            -c {wildcards.counts} \
            -t {wildcards.tpm} \
            -o {output.spi} \
            -g {input.spi}
        """
