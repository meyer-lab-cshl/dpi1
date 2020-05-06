# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint split:
    input:
        alltc="{dir}/outPooled/tc.long.bed.gz"
    output:
        files=directory("{dir}/outPooled/tc.long")
    shell:
        """
        # unzip and split file into small chunks; named by `suffix_length` letters
        # to form the suffix of the file name; `lines` per file; prefix `tc`
        gunzip -c {input.alltc} | \
        gsplit \
            --suffix-length=5 \
            --lines=1000  - {wildcards.dir}/outPooled/tc.long/
        """


rule decompose:
    input:
        tc="{dir}/outPooled/tc.long/tc.{suffix}",
    output:
        bed="{dir}/outPooled/tc.long/{suffix}.decompose_smoothing.bed"
    wildcard_constraints:
        suffix="!gz$|!bed$"
    params:
        pattern=".bw$",
        length=config["decomp_length"],
        ratio=config["noise_subtraction_ratio"],
        window=config["gaussian_window_size_half"],
        bound=config["n_comp_upper_bound"],
    shell:
        """
        Rscript scripts/decompose_peakclustering_hvm.R \
            --tagclusters={input.tc} \
            --path={wildcards.dir}/outCounts \
            --pattern={params.pattern} \
            --exclude_prefix = xxxxxxxxxxxx\. \
            --gaussian_window_size_half={params.window} \
            --n.comp.upper_bound={params.bound} \
            --length_to_decompose={params.length} \
            --noise_subtraction_ratio={params.ratio} \
            --analysis=dpi
        """

def preparesmooth_input(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    return expand("{pdir}/outPooled/tc.long/{suffix}.decompose_smoothing.bed",
           pdir=wildcards.dir,
           suffix=glob_wildcards(os.path.join(checkpoint_output, "{suffix}")).suffix)


def ica_input(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    return expand("{pdir}/outPooled/tc.long/{suffix}.ica.txt",
           pdir=wildcards.dir,
           suffix=glob_wildcards(os.path.join(checkpoint_output, "{suffix}")).suffix)

rule ica:
    input:
        ica_input
    output:
        bedgraph = "{dir}/outPooled/tc.long.decompose_smoothing.component{n}_ctss.{strand}.bedGraph.gz"
    params:
        strand = lambda wildcards:  "+" if wildcards.strand == "fwd" else "-",
        nplus1 = lambda wildcards: int(wildcards.n) + 1
    run:
        """
        cat {input} | \
        awk -v n={wildcards.n} '(NF > n) {print}'| \
        cut -f 1,{params.nplus1} | \
        sed 's/:/\t/;s/,/\t/;s/\\.\\./\t/' \
        awk -v strand={params.strand} '{if (($4 == strand) && ($5 \!= "NA")) \
            {printf "%s\t%i\t%i\t%i\n", $1,$2,$3,$5} }' | \
        sort -k 1,1 -k 2,2n | \
        gzip -c > {output.bedgraph}
        """

rule preparesmooth:
    input:
        preparesmooth_input
    output:
    #   temp(ica="{dir}/outPooled/tc.long.{suffix}.ica.tmp")
        temp(bed="{dir}/outPooled/tc.long.decompose_smoothing.bed.tmp")
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
    params:
        colors = lambda wildcards:  "255,0,0" if wildcards.strand == "+" else "0,0,255",
    run:
        import re
        colors = {"+": "255,0,0", "-": "0,0,255")
        with open({output}, "w") as op:
        with open({input.bed}) as ip:
            lines = ip.readlines()
        for line in lines:
            #chr1 566854 566953 chr1:566727..567119,+;peak:chr1:566874..566875,+;peakCounts:348 1000 +
            #chr1 909815 910037 chr1:909815..910037,+ 1000 +
            cols = line.split("\t")
            p = re.compile('peak:chr.*:(\d+)\.\.(\d+),')
            if m = p.search(cols[3]):
                cols.append(m.groups()[0])
                cols.append(m.groups()[1])
            else:
                cols.append(cols[1])
                cols.append(cols[2])
            cols.append(colors[cols[5])
            op.writelines("\t".join(cols))

rule strand:
    input:
        bed="{dir}/outPooled/tc.long.decompose_smoothing.bed.gz",
        counts="{dir}/outPooled/ctssTotalCounts.{strand}.bw"
    output:
        temp("{dir}/outPooled/tc.long.decompose_smoothing_merged.{strand}.bed")
    params:
        sign = lambda wildcards:  "+" if wildcards.strand == "fwd" else "-",
    shell:
        """
        gunzip -c {input.bed} |
        awk '{if($6 == "+"){print $1,$2,$3,$4,$5,$6}}'| \
        mergeBed -i stdin | \
        awk -v strand={params.sign} 'BEGIN{OFS="\t"} \
            {print $1,$2,$3, $1":"$2".."$3","strand, 1000,strand }' | \
        scripts/bed2peakBed9.sh \
            -b {input.counts} \
            -c "255,0,0" \
          > {output}
        """

rule merge_strand:
    input:
        fwd="{dir}/outPooled/tc.long.decompose_smoothing_merged.fwd.bed",
        rev="{dir}/outPooled/tc.long.decompose_smoothing_merged.rev.bed",
    output:
        "{dir}/outPooled/tc.long.decompose_smoothing_merged.bed"
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
        awk 'BEGIN{OFS="\t"}{print $0,1,$3-$2",",0","}' | \
        gzip -c > {output}
        """


