#!/usr/bin/env python3

import multiprocessing


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
        'pipeline-output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv'

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

rule freebayes_csd:
    input:
        bam = expand('output/030_process-aln/{indiv}_marked.bam',
                     indiv=all_indivs),
        bai = expand('output/030_process-aln/{indiv}_marked.bam.bai',
                     indiv=all_indivs),
        fa = honeybee_ref
    output:
        vcf = 'output/040_freebayes/variants_csd.vcf'
    params:
        ploidy = '2',
        region = 'NC_037640.1:11771679-11781139'
    log:
        'output/logs/040_freebayes/freebayes.log'
    singularity:
        freebayes_container
    shell:
        'freebayes '
        '--ploidy {params.ploidy} '
        '--region {params.region} '
        '-f {input.fa} '
        '{input.bam} '
        '> {output} '
        '2> {log}'


rule merge_bam: # for visualisation
    input:
        expand('output/030_process-aln/{indiv}_marked.bam',
               indiv=all_indivs)
    output:
        'output/030_process-aln/merged.bam'
    log:
        'output/logs/030_process-aln/merge.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        sambamba_container
    shell:
        'sambamba merge '
        '--nthreads={threads} '
        '--compression-level=9 '
        '{output} '
        '{input} '
        '2> {log} '
        '; '
        'sambamba index {output} 2>> {log}'


rule index_bamfile:
    input:
        'output/030_process-aln/{indiv}_marked.bam'
    output:
        'output/030_process-aln/{indiv}_marked.bam.bai'
    log:
        'output/logs/030_process-aln/{indiv}_index.log',
    threads:
        2
    singularity:
        samtools_container
    shell:
        'samtools index -@ {threads} {input} 2> {log}'

rule markdup:
    input:
        'output/020_bwa/{indiv}.sam'
    output:
        sorted = temp('output/030_process-aln/{indiv}_sorted.bam'),
        marked = 'output/030_process-aln/{indiv}_marked.bam'
    threads:
        multiprocessing.cpu_count()
    log:
        s = 'output/logs/030_process-aln/{indiv}_sort.log',
        m = 'output/logs/030_process-aln/{indiv}_markdup.log'
    singularity:
        samtools_container
    shell:
        'samtools fixmate '
        '-m '
        '-O BAM '
        '-@ {threads} '
        '{input} '
        '- '
        '2> {log.s} '
        '| '
        'samtools sort '
        '-o {output.sorted} '
        '-O BAM '
        '-l 0 '
        '-@ {threads} '
        '- '
        '2>> {log.s} '
        '; '
        'samtools markdup '
        '-@ {threads} '
        '-s '
        '{output.sorted} '
        '{output.marked} '
        '2> {log.m}'

# map individuals
rule bwa:
    input:
        fq = 'output/010_trim-decon/{indiv}.fq.gz',
        index = expand('output/020_bwa/honeybee_ref.fasta.{suffix}',
                       suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        temp('output/020_bwa/{indiv}.sam')
    params:
        prefix = 'output/020_bwa/honeybee_ref.fasta',
        rg = '\'@RG\\tID:{indiv}\\tSM:{indiv}\''
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/020_bwa/{indiv}.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '-p '
        '-R {params.rg} '
        '{params.prefix} '
        '{input.fq} '
        '> {output} '
        '2> {log}'

# prepare ref
rule index:
    input:
        honeybee_ref
    output:
        expand('output/020_bwa/honeybee_ref.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/020_bwa/honeybee_ref.fasta'
    threads:
        1
    log:
        'output/logs/020_bwa/index.log'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input} '
        '2> {log}'

rule trim_decon:
    input:
        r1 = 'data/raw/{indiv}/{indiv}_R1.fq.gz',
        r2 = 'data/raw/{indiv}/{indiv}_R2.fq.gz'
    output:
        fq = 'output/010_trim-decon/{indiv}.fq.gz',
        f_stats = 'output/010_trim-decon/{indiv}_filter-stats.txt',
        t_stats = 'output/010_trim-decon/{indiv}_trim-stats.txt'
    log:
        filter = 'output/logs/010_trim-decon/{indiv}_filter.log',
        trim = 'output/logs/010_trim-decon/{indiv}_trim.log',
        repair1 = 'output/logs/010_trim-decon/{indiv}_repair1.log',
        repair2 = 'output/logs/010_trim-decon/{indiv}_repair2.log'
    params:
        filter = bbduk_ref,
        trim = bbduk_adaptors
    threads:
        multiprocessing.cpu_count()
    singularity:
        bbduk_container
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '2> {log.repair1} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.f_stats} '
        '2> {log.filter} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '
        '| '
        'repair.sh '
        'in=stdin.fastq '
        'out={output.fq} '
        '2> {log.repair2} '

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


# try the pipeline
rule genotype:
    input:
        csv = 'data/ty_samples.csv',
        ref = honeybee_ref
    output:
        'pipeline-output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv'
    params:
        wd = 'pipeline-output/010_genotypes',
        ploidy = '2'
    log:
        'pipeline-output/logs/genotype.log'
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


