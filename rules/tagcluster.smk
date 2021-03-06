rule create_stranded_tagcluster:
    input:
        total_counts="{dir}/outPooled/ctssTotalCounts.{strand}.bw",
        max_counts="{dir}/outPooled/ctssMaxCounts.{strand}.bw",
    output:
        tc="{dir}/outPooled/tc.{strand}.bed"
    resources:
        mem_mb = 4000
    conda:
        "../envs/tagcluster.yaml"
    wildcard_constraints:
        strand="fwd|rev"
    params:
        sign=lambda wildcards: '+' if wildcards.strand == "fwd" else '-',
        thr=config["cluster"]["ctss_max_counts"],
        dist=config["cluster"]["dist"]
    shell:
        """
        scripts/prepare_tagclusters.sh \
            -m {input.max_counts} \
            -t {input.total_counts} \
            -c {params.thr} \
            -d {params.dist} \
            -s {params.sign} \
            -o {output.tc}
        """

rule combine_stranded_tagcluster:
    input:
        fwd="{dir}/outPooled/tc.fwd.bed",
        rev="{dir}/outPooled/tc.rev.bed"
    output:
        tc="{dir}/outPooled/tc.bed.gz"
    resources:
        mem_mb = 4000
    shell:
        """
        scripts/combine_tagclusters.sh \
            -f {input.fwd} \
            -r {input.rev} \
            -o {output.tc}
        """
