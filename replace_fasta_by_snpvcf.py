#-*- coding:utf-8 -*-
import re
import vcf
import pandas as pd
from Bio import SeqIO
vcf_reader = vcf\
        .Reader(open('dataset/ID-8.Q30trimQ30filMQ25.vcf.snp.eff.pass.imp', 'r'))
tmp = {i: [record.CHROM, record.POS, record.REF, record.ALT] \
       for i, record in enumerate(vcf_reader)}
mutation_list = pd.DataFrame.from_dict(tmp, orient='index')
mutation_list.rename(columns={0:'chr',1:'pos',2:'ref',3:'alt'},inplace=True)
mutation_list['alt'] = mutation_list.alt.apply(lambda x: x[0])
# altが2文字になっている箇所がある
mutation_list = mutation_list[mutation_list.alt.apply(lambda x: len(x) == 1)]

records = []
for record in SeqIO.parse("dataset/reference.fa", "fasta"):
    if not record.id.startswith('chr'):
        print(record.id, 'None')
    else:
        t_mut = mutation_list[mutation_list.chr == \
                              re.match("chr(.*)", record.id).group(1)]
        mutable_seq = record.seq.tomutable()
        count = 0
        for i, mut in t_mut.iterrows():
            if mutable_seq[mut.pos-1] == mut.ref:
                mutable_seq[mut.pos-1] = str(mut.alt)
                count += 1
            else:
                print("Ref Base doesn'tmatch at "+str(mut.pos-1))
                raise
        print(record.id, len(t_mut), count)
        record.seq = mutable_seq.toseq()
    records.append(record)
with open("dataset/reference_replace_snp.fa","w") as f:
    SeqIO.write(records, f, "fasta")
