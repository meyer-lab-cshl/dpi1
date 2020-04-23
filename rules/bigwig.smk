rule prepare_bigwig_counts:
    input:
        sample="{dir}/{sample}.ctss.bed.gz",
        genome=config["genome"]
    output:
        "{dir}/outCounts/{sample}.ctss.fwd.bw",
        "{dir}/outCounts/{sample}.ctss.rev.bw"
    conda:
        "../envs/bigwig.yaml"
    shell:
        """
        scripts/ctssBed2bigWig.sh \
            -i {input.sample} \
            -o {wildcards.dir}/outCounts/{wildcards.sample}.ctss \
            -g {input.genome}
        """

rule prepare_bigwig_tpm:
    input:
        sample="{dir}/{sample}.ctss.bed.gz",
        genome=config["genome"]
    output:
        "{dir}/outTpm/{sample}.ctss.fwd.bw",
        "{dir}/outTpm/{sample}.ctss.rev.bw"
    conda:
        "../envs/bigwig.yaml"
    shell:
        """
        scripts/ctssBed2TpmBigWig.sh \
            -i {input.sample} \
            -o {wildcards.dir}/outTpm/{wildcards.sample}.ctss \
            -g {input.genome}
        """
