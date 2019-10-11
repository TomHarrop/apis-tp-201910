import multiprocessing

bbduk_ref = '/phix174_ill.ref.fa.gz'
bbduk_adaptors = '/adapters.fa'

# containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

all_indivs = sorted(set(glob_wildcards("data/raw/{indiv}/{indiv1}_R1.fq.gz").indiv))

rule target:
    input:
        expand('output/010_trim-decon/{indiv}.fq.gz',
               indiv=all_indivs)

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
