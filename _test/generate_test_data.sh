dir=/Users/hannah/data/tss/mouse/fantom/bed/GRCm38
f1=$dir/ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep1.CNhs14104.14357-155I1.mm10.nobarcode.ctss.bed.gz
f2=$dir/ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep2.CNhs14109.14362-155I6.mm10.nobarcode.ctss.bed.gz

lines=2000
gunzip -c $f1 |head -n $lines |
    gzip -c > rake/test.1.ctss.bed.gz
gunzip -c $f1 |tail -n $lines |
    gzip -c >> rake/test.1.ctss.bed.gz

gunzip -c $f2 |head -n $lines |
    gzip -c > rake/test.2.ctss.bed.gz
gunzip -c $f2 |tail -n $lines |
    gzip -c >> rake/test.2.ctss.bed.gz

gunzip -c $f1 |head -n $lines |
    gzip -c > snakemake/test.1.ctss.bed.gz
gunzip -c $f1 |tail -n $lines |
    gzip -c >> snakemake/test.1.ctss.bed.gz

gunzip -c $f2 |head -n $lines |
    gzip -c > snakemake/test.2.ctss.bed.gz
gunzip -c $f2 |tail -n $lines |
    gzip -c >> snakemake/test.2.ctss.bed.gz
