#!/usr/bin/env python3

import multiprocessing
import pandas


def get_min_cutoff(wildcards):
    cutoff_file = checkpoints.genotype.get(**wildcards).output['cutoffs']
    cutoffs = pandas.read_csv(cutoff_file,
                              header=None,
                              index_col=0)
    return cutoffs.loc['min_depth', 1]


def get_max_cutoff(wildcards):
    cutoff_file = checkpoints.genotype.get(**wildcards).output['cutoffs']
    cutoffs = pandas.read_csv(cutoff_file,
                              header=None,
                              index_col=0)
    return cutoffs.loc['max_depth', 1]


###########
# GLOBALS #
###########

bbduk_ref = '/phix174_ill.ref.fa.gz'
bbduk_adaptors = '/adapters.fa'
honeybee_ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'

# containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
bioconductor_container = 'shub://TomHarrop/singularity-containers:bioconductor_3.9'
biopython_container = 'shub://TomHarrop/singularity-containers:biopython_1.73'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
clustalo = 'shub://TomHarrop/singularity-containers:clustalo_1.2.4'
freebayes_container = 'shub://TomHarrop/singularity-containers:freebayes_1.2.0'
pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
sambamba_container = 'shub://TomHarrop/singularity-containers:sambamba_0.6.9'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
vcflib_container = 'shub://TomHarrop/singularity-containers:vcflib_1.0.0-rc2'
honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.2'
    '@0726020591078c43cd45c1e0d138928d3dc197ce')
vcftools = ('shub://TomHarrop/variant-utils:vcftools_0.1.16'
            '@d64cc5a37951760be575c43024c66e69b2563166')
whatshap = 'shub://TomHarrop/variant-utils:whatshap_0.18'

########
# MAIN #
########

all_indivs = sorted(set(
    glob_wildcards("data/raw/{indiv}/{indiv1}_R1.fq.gz").indiv))


#########
# RULES #
#########

rule target:
    input:
        # 'output/050_variant-annotation/csd.vcf',
        # 'output/030_process-aln/merged.bam',
        # 'output/060_derived-alleles/all-indivs_aa.fa',
        # 'output/060_derived-alleles/all-indivs_aa.faa',
        # 'output/060_derived-alleles/drones_aa.f
        'output/030_phased/phased.vcf.gz'

rule align_consensus:
    input:
        'output/060_derived-alleles/all-indivs_aa.fa'
    output:
        aln = 'output/060_derived-alleles/all-indivs_aa.faa',
        dist = 'output/060_derived-alleles/all-indivs_aa.dist'
    log:
        'output/logs/060_derived-alleles/clustalo.fa'
    threads:
        multiprocessing.cpu_count()
    singularity:
        clustalo
    shell:
        'clustalo '
        '-i {input} '
        '--threads {threads} '
        '--dealign '
        '--full '
        '--out {output.aln} '
        '--distmat-out {output.dist} '
        '&> {log}'

rule drones_only:
    input:
        'output/060_derived-alleles/all-indivs_aa.fa'
    output:
        'output/060_derived-alleles/drones_aa.fa'
    log:
        'output/logs/060_derived-alleles/filterbyname.log'
    threads:
        1
    singularity:
        bbduk_container
    shell:
        'filterbyname.sh '
        'in={input} '
        'out={output} '
        'names=DR '
        'substring=name '
        'include=t '
        'ignorejunk=t '
        '2> {log}'

rule translate_consensus:
    input:
        'output/060_derived-alleles/consensus_all-indivs.fa'
    output:
        'output/060_derived-alleles/all-indivs_aa.fa'
    singularity:
        biopython_container
    script:
        'src/translate_consensus.py'

rule combine_cds:
    input:
        expand('output/060_derived-alleles/{indiv}_consensus_condensed.fa',
               indiv=all_indivs)
    output:
        'output/060_derived-alleles/consensus_all-indivs.fa'
    singularity:
        samtools_container
    shell:
        'cat {input} > {output}'

rule condense_cds:
    input:
        'output/060_derived-alleles/{indiv}-consensus.fa'
    output:
        temp('output/060_derived-alleles/{indiv}_consensus_condensed.fa')
    params:
        header = '>{indiv}'
    singularity:
        samtools_container
    shell:
        'echo "{params.header}" > {output} ; '
        'grep -v "^>" {input} >> {output}'


rule extract_derived_cds:
    input:
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
        regions = 'output/050_variant-annotation/hvr_dt.txt',
        vcf = 'output/050_variant-annotation/csd_reheadered.vcf.gz'
    output:
        temp('output/060_derived-alleles/{indiv}-consensus.fa')
    log:
        'output/logs/060_derived-alleles/{indiv}-consensus.log'
    singularity:
        samtools_container
    shell:
        'samtools faidx '
        '{input.fa} '
        '$(cat {input.regions}) '
        '2> {log} '
        '| '
        'bcftools consensus '
        '-s {wildcards.indiv} '
        '-H 1 '
        '{input.vcf} '
        '> {output} '
        '2>> {log}'

rule filter_csd_variants:
    input:
        gff = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.gff',
        vcf = 'output/040_freebayes/variants_csd_filtered.vcf.gz',
        tbi = 'output/040_freebayes/variants_csd_filtered.vcf.gz.tbi',
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
    output:
        coding = 'output/050_variant-annotation/coding.Rds',
        csd = 'output/050_variant-annotation/csd.vcf',
    log:
        'output/logs/050_variant-annotation/filter_csd_variants.log'
    singularity:
        bioconductor_container
    script:
        'src/filter_csd_variants.R'

rule filter_vcf_csd:
    input:
        'output/040_freebayes/variants_csd.vcf'
    output:
        'output/040_freebayes/variants_csd_filtered.vcf'
    params:
        filter = 'QUAL > 20'
    log:
        'output/logs/040_freebayes/freebayes_filter_csd.log'
    singularity:
        vcflib_container
    shell:
        'vcffilter -f \'{params.filter}\' {input} > {output} 2> {log}'


rule phase:
    input:
        vcf = 'output/020_filtered-genotypes/filtered.vcf.gz',
        bam = 'output/010_genotypes/merged.bam',
        ref = 'output/010_genotypes/015_ref/ref.fasta',
        fai = 'output/010_genotypes/015_ref/ref.fasta.fai',
    output:
        'output/030_phased/phased.vcf'
    log:
        'output/logs/030_phase.log'
    singularity:
        whatshap
    shell:
        'whatshap phase '
        '--reference {input.ref} '
        '-o {output} '
        '--indels '
        '{input.vcf} '
        '{input.bam} '
        '&> {log}'


rule filter:
    input:
        cutoffs = 'output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/calls.vcf.gz'
    output:
        'output/020_filtered-genotypes/filtered.vcf'
    params:
        min_depth = get_min_cutoff,
        max_depth = get_max_cutoff,
        maf = 0.1,
        max_missing = 0.9,
        qual = 30
    log:
        'output/logs/020_filter.log'
    singularity:
        vcftools
    shell:
        'vcftools '
        '--gzvcf {input.vcf} '
        '--maf {params.maf} '
        '--max-missing {params.max_missing} '
        '--minQ {params.qual} '
        '--min-meanDP {params.min_depth} '
        '--max-meanDP {params.max_depth} '
        '--minDP {params.min_depth} '
        '--maxDP {params.max_depth} '
        '--recode '
        '--stdout '
        '> {output} '
        '2> {log}'


# try the pipeline
checkpoint genotype:
    input:
        csv = 'data/ty_samples.csv',
        ref = honeybee_ref
    output:
        cutoffs = 'output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/calls.vcf.gz',
        bam = 'output/010_genotypes/merged.bam',
        ref = 'output/010_genotypes/015_ref/ref.fasta',
        fai = 'output/010_genotypes/015_ref/ref.fasta.fai',
    params:
        wd = 'output/010_genotypes',
        ploidy = '2'
    log:
        'output/logs/genotype.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        honeybee_genotype_pipeline
    shell:
        'honeybee_genotype_pipeline '
        '--ref {input.ref} '
        '--samples_csv {input.csv} '
        '--outdir {params.wd} '
        '--ploidy {params.ploidy} '
        '--threads {threads} '
        '--restart_times 1 '
        '&> {log}'

# generic csd rules
rule extract_hvr_exon:
    input:
        gff = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.gff'
    output:
        regions = 'output/050_variant-annotation/hvr_dt.txt'
    log:
        'output/logs/050_variant-annotation/extract_hvr_exon.log'
    singularity:
        bioconductor_container
    script:
        'src/extract_hvr_exon.R'

# generic index rule
rule index_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        gz = 'output/{folder}/{file}.vcf.gz',
        tbi = 'output/{folder}/{file}.vcf.gz.tbi'
    log:
        'output/logs/{folder}/{file}_index-vcf.log'
    singularity:
        samtools_container
    shell:
        'bgzip -c {input} > {output.gz} 2> {log} '
        '; '
        'tabix -p vcf {output.gz} 2>> {log}'

# generic reheader rule
rule reheader_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        'output/{folder}/{file}_reheadered.vcf'
    singularity:
        samtools_container
    shell:
        'grep -v "^##FILTER=All filters passed" {input} > {output}'

