import generic as gen
import collections
import re
import os
import numpy as np
import pandas as pd
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table
from mmsplice.utils import max_varEff

description = "File to ."
arguments = ["vcf", "gtf", "fasta", "output_file"]
args = gen.parse_arguments(description, arguments, flags = [], ints=[])

vcf, gtf, fasta, output_file = args.vcf, args.gtf, args.fasta, args.output_file

# dataloader to load variants from vcf
dl = SplicingVCFDataloader(gtf, fasta, vcf, split_seq=False)

entries = [i for i in gen.read_many_fields(vcf, "\t") if not i[0].startswith("#")]
print("\n{0} PTCS to predict.\n".format(len(entries)))

# Specify model
model = MMSplice(
    exon_cut_l=0,
    exon_cut_r=0,
    acceptor_intron_cut=6,
    donor_intron_cut=6,
    acceptor_intron_len=50,
    acceptor_exon_len=3,
    donor_exon_len=5,
    donor_intron_len=13)

# Do prediction
predictions = predict_all_table(model, dl, batch_size=1024, split_seq=False, assembly=False)

# Summerize with maximum effect size
predictionsMax = max_varEff(predictions)
# write to output file
predictionsMax.to_csv(path_or_buf = output_file)
