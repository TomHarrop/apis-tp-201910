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
freebayes_container = 'shub://TomHarrop/singularity-containers:freebayes_1.2.0'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
vcflib_container = 'shub://TomHarrop/singularity-containers:vcflib_1.0.0-rc2'
bioconductor_container = 'shub://TomHarrop/singularity-containers:bioconductor_3.9'
pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'


########
# MAIN #
########

all_indivs = sorted(set(glob_wildcards("data/raw/{indiv}/{indiv1}_R1.fq.gz").indiv))


#########
# RULES #
#########


rule target:
    input:
        'output/040_freebayes/variants_filtered.vcf'

rule filter_vcf:
    input:
        'output/040_freebayes/variants.vcf'
    output:
        'output/040_freebayes/variants_filtered.vcf'
    params:
        filter = 'QUAL > 20'
    log:
        'output/logs/040_freebayes/freebayes_filter.log'
    singularity:
        vcflib_container
    shell:
        'vcffilter -f \'{params.filter}\' {input} > {output} 2> {log}'

rule freebayes:
    input:
        bam = expand('output/030_process-aln/{indiv}_marked.bam',
                     indiv=all_indivs),
        bai = expand('output/030_process-aln/{indiv}_marked.bam.bai',
                     indiv=all_indivs),
        fa = honeybee_ref
    output:
        vcf = 'output/040_freebayes/variants.vcf'
    params:
        ploidy = '2'
    log:
        'output/logs/040_freebayes/freebayes.log'
    singularity:
        freebayes_container
    shell:
        'freebayes '
        '--ploidy {params.ploidy} '
        '-f {input.fa} '
        '{input.bam} '
        '> {output} '
        '2> {log}'


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
        'output/020_bwa/{indiv}.sam'
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
