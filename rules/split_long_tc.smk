rule split_long:
    input:
        alltc="{dir}/outPooled/tc.long.bed.gz"
    output:
        dynamic("{dir}/outPooled/tc.long/tc.{suffix}"),
    shell:
        """
        # unzip and split file into small chunks; named by `suffix_length` letters
        # to form the suffix of the file name; `lines` per file; prefix `tc`
        gunzip -c {output.alltc} | \
        gsplit \
            --suffix-length=5 \
            --lines=1000  - {wildcards.dir/outPooled/tc.long/
        """
