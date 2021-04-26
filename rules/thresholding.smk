rule permissive:
    input:
        fwd_counts="{dir}/outPooled/ctssMaxCounts.fwd.bw",
        rev_counts="{dir}/outPooled/ctssMaxCounts.rev.bw",
        dpi="{dir}/outPooled/tc.decompose_smoothing_merged.bed.gz"
    output:
        dpi="{dir}/outPooled/tc.decompose_smoothing_merged.ctssMaxCounts{counts}.bed.gz"
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
            -o {output.dpi} \
            -g {input.dpi}
        """

rule robust:
    input:
        fwd_counts="{dir}/outPooled/ctssMaxCounts.fwd.bw",
        rev_counts="{dir}/outPooled/ctssMaxCounts.rev.bw",
        fwd_tpm="{dir}/outPooled/ctssMaxTpm.fwd.bw",
        rev_tpm="{dir}/outPooled/ctssMaxTpm.rev.bw",
        dpi="{dir}/outPooled/tc.decompose_smoothing_merged.bed.gz"
    output:
        dpi="{dir}/outPooled/tc.decompose_smoothing_merged.ctssMaxCounts{counts}.ctssMaxTpm{tpm}.bed.gz"
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
            -o {output.dpi} \
            -g {input.dpi}
        """
