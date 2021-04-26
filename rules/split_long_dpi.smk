rule split_long:
    input:
        tc="{dir}/outPooled/tc.bed.gz",
        genome=config["genome"]
    output:
        #dynamic("{dir}/outPooled/tc.long/tc.{suffix}"),
        alltc="{dir}/outPooled/tc.long.bed.gz"
    params:
        length=config["split"]["length"],
        ct=config["split"]["count"],
        ratio=config["split"]["noise_subtraction_ratio"],
        chunk_length=config["decomposition"]["chunk"],
        suffix_length=config["decomposition"]["suffix"],
    shell:
        """
        # Read tag clusters; filter by > length and count; zip to file
        gunzip -c {input.tc} | \
        awk '{{if((($3 - $2) > {params.length}) || ($5 > {params.ct})){{print}}}}' | \
        gzip -c > {output.alltc}
        # unzip and split file into small chunks; named by `suffix_length` letters
        # to form the suffix of the file name; `lines` per file; prefix `tc`
        gunzip -c {output.alltc} | \
        gsplit \
            --suffix-length={params.suffix_length} \
            --lines={params.chunk_length} \
            - {wildcards.dir}/outPooled/tc.long/
        """
