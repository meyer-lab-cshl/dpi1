# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint split:
    input:
        alltc="{dir}/outPooled/tc.long.bed.gz"
    output:
        directory("{dir}/outPooled/tc.long.files")
    shell:
        """
        mkdir {wildcards.dir}/outPooled/tc.long.files;
        # unzip and split file into small chunks; named by `suffix_length` letters
        # to form the suffix of the file name; `lines` per file; prefix `tc`
        gunzip -c {input.alltc} | \
        gsplit \
            --suffix-length=5 \
            --lines=1000  - {wildcards.dir}/outPooled/tc.long.files/
        """


rule decompose:
    input:
        tc="{dir}/outPooled/tc.long.files/{suffix}",
    output:
        bed="{dir}/outPooled/tc.long.decompose/{suffix}.decompose_smoothing.bed",
        ica="{dir}/outPooled/tc.long.decompose/{suffix}.ica.txt"
    params:
        pattern=".bw$",
        length=config["split"]["length"],
        ratio=config["split"]["noise_subtraction_ratio"],
        bound=config["split"]["n_comp_upper_bound"],
        window=config["smoothing"]["gaussian_window_size_half"],
    shell:
        """
        Rscript scripts/decompose_peakclustering_hvm.R \
            --icafile={output.ica} \
            --outfile={output.bed} \
            --tagclusters={input.tc} \
            --path={wildcards.dir}/outCounts \
            --pattern={params.pattern} \
            --exclude=xxxxxxxxxxxx\. \
            --window={params.window} \
            --bound={params.bound} \
            --length={params.length} \
            --ratio={params.ratio} \
            --analysis=dpi
        """

def preparesmooth_input(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    suffix=glob_wildcards(os.path.join(checkpoint_output, "{suffix}")).suffix
    suffix=[a for a in suffix if not a.startswith(".")]
    return expand("{pdir}/outPooled/tc.long.decompose/{suffix}.decompose_smoothing.bed",
           pdir=wildcards.dir,
           suffix=suffix)

rule preparesmooth:
    input:
        preparesmooth_input
    output:
        bed=temp("{dir}/outPooled/tc.long.decompose_smoothing.bed.tmp")
    shell:
        """
        cat {input} | \
        sort -k 1,1 -k 2,2n  > {output.bed}
        """

rule smooth:
    input:
        bed="{dir}/outPooled/tc.long.decompose_smoothing.bed.tmp"
    output:
        bed="{dir}/outPooled/tc.long.decompose_smoothing.bed.gz"
    run:
        import gzip
        import re
        colors = {"+": "255,0,0", "-": "0,0,255"}
        with gzip.open(output.bed, "wt") as op:
            with open(input.bed) as ip:
                lines = ip.readlines()
                for line in lines:
                    #chr1 566854 566953 chr1:566727..567119,+;peak:chr1:566874..566875,+;peakCounts:348 1000 +
                    #chr1 909815 910037 chr1:909815..910037,+ 1000 +
                    cols = line.rstrip().split("\t")
                    p = re.compile('peak:chr.*:(\d+)\.\.(\d+),')
                    m = p.search(cols[3])
                    if m:
                        cols.append(m.groups()[0])
                        cols.append(m.groups()[1])
                    else:
                        cols.append(cols[1])
                        cols.append(cols[2])
                    cols.append(colors[cols[5]])
                    op.writelines("\t".join(cols) + '\n')

rule strand:
    input:
        bed="{dir}/outPooled/tc.long.decompose_smoothing.bed.gz",
        counts="{dir}/outPooled/ctssTotalCounts.{strand}.bw"
    output:
        temp("{dir}/outPooled/tc.long.decompose_smoothing_merged.{strand}.bed")
    params:
        sign = lambda wildcards:  "+" if wildcards.strand == "fwd" else "-",
        color = lambda wildcards:  "255,0,0" if wildcards.strand == "fwd" else "0,0,255",
    shell:
        """
        scripts/split_decomposition_strand.sh \
            -b {input.bed} \
            -c {input.counts} \
            -s {params.sign} \
            -C {params.color} \
            -o {output}
        """

rule merge_strand:
    input:
        fwd="{dir}/outPooled/tc.long.decompose_smoothing_merged.fwd.bed",
        rev="{dir}/outPooled/tc.long.decompose_smoothing_merged.rev.bed",
    output:
        "{dir}/outPooled/tc.long.decompose_smoothing_merged.bed.gz"
    shell:
        """
        cat {input.fwd} {input.rev} | \
        sort -k1,1 -k2,2n | \
        gzip -c > {output}
        """

rule merge_all_bed:
    input:
        bedlong="{dir}/outPooled/tc.long.decompose_smoothing_merged.bed.gz",
        bedshort="{dir}/outPooled/tc.short.bed.gz"
    output:
        "{dir}/outPooled/tc.decompose_smoothing_merged.bed.gz",
    shell:
        """
        gunzip -c {input.bedlong} {input.bedshort} | \
        sort -k1,1 -k2,2n -k6,6 | \
        awk 'BEGIN{{OFS="\t"}}{{print $0,1,$3-$2",",0","}}' | \
        gzip -c > {output}
        """


def ica_input(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    suffix=glob_wildcards(os.path.join(checkpoint_output, "{suffix}")).suffix
    suffix=[a for a in suffix if not a.startswith(".")]
    return expand("{pdir}/outPooled/tc.long.decompose/{suffix}.ica.txt",
           pdir=wildcards.dir,
           suffix=suffix)

rule ica:
    input:
        ica_input
    output:
        bedgraph = "{dir}/outPooled/tc.long.decompose_smoothing.component{n}_ctss.{strand}.bedGraph.gz"
    params:
        strand = lambda wildcards:  "+" if wildcards.strand == "fwd" else "-",
        nplus1 = lambda wildcards: int(wildcards.n) + 1
    log:
        "{dir}/logs/ica.component{n}.{strand}.log"
    shell:
        """
        touch {output.bedgraph}
        cat {input} | \
        awk -v cols={wildcards.n} '(NF > cols) {{print}}'| \
        cut -f 1,{params.nplus1} | \
        sed 's/:/ /;s/,/ /;s/\\.\\./ /' | \
        tr ' ' '\t' |
        awk -v strand={params.strand} '{{if (($4 == strand) && ($5 != "NA")) \
            {{printf "%s\\t%i\\t%i\\t%i\\n", $1,$2,$3,$5}} }}' | \
        sort -k 1,1 -k 2,2n | \
        gzip -c > {output.bedgraph} 2>{log}
        """
