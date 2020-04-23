rule smoothing:
    input:
        tc_long="{dir}/outPooled/tc.long.bed.gz",
        fwd="{dir}/outPooled/ctssTotalCounts.fwd.bw",
        rev="{dir}/outPooled/ctssTotalCounts.rev.bw",
    output:
        spi="{dir}/outPooled/tc.long.spi.bed.gz"
    params:
        length=config["decomposition"]["length"],
        ratio=config["decomposition"]["noise_subtraction_ratio"],
        window=config["smoothing"]["gaussian_window_size_half"],
    shell:
        """
        Rscript scripts/decompose_peakclustering_hvm.R \
            --tagclusters {input.tc_long} \
            --ctssprefix {wildcards.dir}/outPooled/ctssTotalCounts \
            --window {params.window} \
            --length {params.length} \
            --ratio {params.ratio} \
            --outfile {wildcards.dir}/outPooled/tc.long.spi.bed \
            --analysis spi
        gzip {wildcards.dir}/outPooled/tc.long.spi.bed
        """
