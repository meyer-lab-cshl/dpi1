rule prepare_tag_cluster:
    input:
        total_counts="{dir}/outPooled/ctssTotalCounts.{strand}.bw",
        max_counts="{dir}/outPooled/ctssMaxCounts.{strand}.bw",
        genome=config["genome"]
    output:
        tc="{dir}/outPooled/tc.{strand}.bed"
    params:
        sign=lambda wildcards: '+' if wildcards.strand == "fwd" else '-',
        thr=config["cluster"]["ctss_max_counts"],
        dist=config["cluster"]["dist"]
    shell:
        """
        bigWigToBedGraph {input.max_counts} /dev/stdout | \
        awk -v thr={params.thr} 'BEGIN{{OFS="\t"}} \
            {{if($4 >={params.thr}){{print $1,$2,$3,$4}}}}' | \
        scripts/largeMergeBed.sh \
            -g {input.genome} \
            -d {params.dist} | \
        awk -v sign={params.sign} 'BEGIN{{OFS="\t"}} \
            {{print $1,$2,$3,$4","sign,1000,sign}}' | \
        grep -v '_' | \
        bigWigAverageOverBed {input.total_counts} /dev/stdin  /dev/stdout | \
        awk '{{printf "%s\t%i\\n",$1,$4}}' \
        > {output.tc}
        """

rule combine_tag_cluster:
    input:
        fwd="{dir}/outPooled/tc.fwd.bed",
        rev="{dir}/outPooled/tc.rev.bed"
    output:
        tc="{dir}/outPooled/tc.bed.gz"
    shell:
        """
        cat {input.fwd} {input.rev} | \
        sed -e 's/[:|,|..]/\t/g' | \
        awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$1":"$2".."$3","$4,$5,$4}}' | \
        sort -k1,1 -k2,2n | \
        gzip -c > {output.tc}
        """
