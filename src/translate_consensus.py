#!/usr/bin/env python3

from Bio import Seq
from Bio import SeqIO

consensus_file = snakemake.input[0]
# consensus_file = 'output/070_regions/all_indivs_consensus.fa'
# consensus_file = 'output/050_derived-alleles/flongle/BB34_17_consensus.fa'
aa_file = snakemake.output[0]

cds_list = [x for x in SeqIO.parse(consensus_file, 'fasta')]
aa_list = [x[0:-1].reverse_complement().translate(id=x.id,
                                                  name=x.name,
                                                  description='')
           for x in cds_list]
SeqIO.write(aa_list, aa_file, 'fasta')
