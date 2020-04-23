from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.14.0")


##### load config and sample sheets #####
configfile: "config/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

SAMPLE, = glob_wildcards("/Users/hannah/data/tss/mouse/fantom/bed/GRCm38/{sample}.ctss.bed.gz")

##### target rules #####
def get_input(wildcards):
    if config["analysis"] == "spi":
        robust = "{pdir}/outPooled/tc.spi.merged.ctssMaxCounts{robust}.ctssMaxTpm{tpm}.bed.gz".format(
            pdir=config['directory'],
            robust=config['cutoff']['robust'],
            tpm=config['cutoff']['tpm'])
        permissive = "{pdir}/outPooled/tc.spi.merged.ctssMaxCounts{permissive}.bed.gz".format(
            pdir=config['directory'],
            permissive=config['cutoff']['permissive'])
    return [robust, permissive]

rule all:
    input:
        get_input

##### setup report #####
#report: "report/workflow.rst"

##### load rules #####
### prepare bigWig files for individual samples
include: "rules/bigwig.smk"
### pool ctss across samples: total counts, max counts and tpm
include: "rules/pool.smk"
### prepare tag clusters
include: "rules/tagcluster.smk"
### split tag cluster into long and short for decomposition
include: "rules/decomposition.smk"

if config["analysis"] == 'spi':
    ### simple smoothing
    include: "rules/smooth_simple.smk"
    ### combined long and short tag clusters
    include: "rules/merge_simple.smk"
    ### thresholding without decomposition
    include: "rules/thresholding_simple.smk"

