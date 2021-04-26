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
