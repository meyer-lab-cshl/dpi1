rule pool_ctss_total_counts:
    input:
        expand("{{dir}}/outCounts/{sample}.ctss.{strand}.bw",
          sample=samples['sample'],
          strand=['fwd', 'rev']),
        genome=config["genome"]
    output:
        total="{dir}/outPooled/ctssTotalCounts.{strand}.bw"
    resources:
        mem_mb = 10000
    conda:
        "../envs/pool.yaml"
    shell:
        """
        bigWigMerge {wildcards.dir}/outCounts/*{wildcards.strand}.bw /dev/stdout | \
            sort -k1,1 -k2,2n > \
            {output.total}.tmp
        bedGraphToBigWig {output.total}.tmp {input.genome} {output.total}
        rm {output.total}.tmp
        """

rule pool_ctss_max_counts:
    input:
        expand("{{dir}}/outCounts/{sample}.ctss.{strand}.bw",
          sample=samples['sample'],
          strand=['fwd', 'rev']),
        genome=config["genome"]
    output:
        max="{dir}/outPooled/ctssMaxCounts.{strand}.bw"
    resources:
        mem_mb = 10000
    conda:
        "../envs/pool.yaml"
    shell:
        """
        bigWigMerge -max {wildcards.dir}/outCounts/*{wildcards.strand}.bw /dev/stdout | \
            sort -k1,1 -k2,2n > \
            {output.max}.tmp
        bedGraphToBigWig {output.max}.tmp {input.genome} {output.max}
        rm {output.max}.tmp
        """

rule pool_ctss_max_tpm:
    input:
        expand("{{dir}}/outTpm/{sample}.ctss.{strand}.bw",
          sample=samples['sample'],
          strand=['fwd', 'rev']),
        genome=config["genome"]
    output:
        max="{dir}/outPooled/ctssMaxTpm.{strand}.bw"
    resources:
        mem_mb = 10000
    conda:
        "../envs/pool.yaml"
    shell:
        """
        bigWigMerge -max {wildcards.dir}/outTpm/*{wildcards.strand}.bw /dev/stdout | \
            sort -k1,1 -k2,2n > \
            {output.max}.tmp
        bedGraphToBigWig {output.max}.tmp {input.genome} {output.max}
        rm {output.max}.tmp
        """
