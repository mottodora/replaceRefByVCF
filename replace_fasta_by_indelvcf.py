#-*- coding:utf-8 -*-
import re
import vcf
import pandas as pd
from Bio import SeqIO
vcf_reader = vcf\
        .Reader(open('dataset/ID-8.Q30trimQ30filMQ25.vcf.indel.eff.pass.imp', 'r'))
tmp = {i: [record.CHROM, record.POS, record.REF, record.ALT] \
       for i, record in enumerate(vcf_reader)}
mutation_list = pd.DataFrame.from_dict(tmp, orient='index')
mutation_list.rename(columns={0:'chr',1:'pos',2:'ref',3:'alt'},inplace=True)
mutation_list['alt'] = mutation_list.alt.apply(lambda x: x[0])
# altが2文字になっている箇所がある
insertion_list = mutation_list[mutation_list.alt.apply(lambda x: len(x) > 1)]
deletion_list = mutation_list[mutation_list.ref.apply(lambda x: len(x) > 1)]

records = []
for record in SeqIO.parse("dataset/reference.fa", "fasta"):
    if not record.id.startswith('chr'):
        print(record.id, 'None')
    else:
        print(record.id)
        t_ins = insertion_list.loc[insertion_list.chr == \
                              re.match("chr(.*)", record.id).group(1)]
        t_del = deletion_list.loc[deletion_list.chr == \
                              re.match("chr(.*)", record.id).group(1)]
        mutable_seq = record.seq.tomutable()
        del_num = t_del.ref.apply(lambda x: len(x)-1).sum()
        del_target = []
        for _,t_d in t_del.iterrows():
            if mutable_seq[t_d.pos-1] != t_d.ref[0]:
                print("Ref Base doesn't match at "+str(t_d.pos-1))
                raise
            mutable_seq[t_d.pos:t_d.pos+len(t_d.ref)-1] = 'D'*(len(t_d.ref)-1)
            del_target.extend(list(range(t_d.pos,t_d.pos+len(t_d.ref)-1)))
        print(del_target[0])
        for _,t_i in t_ins.sort_values('pos', ascending=False).iterrows():
            if mutable_seq[t_i.pos-1] != t_i.ref:
                print("Ref Base doesn't match at "+str(t_i.pos-1))
                raise
            print(mutable_seq[t_i.pos-1:t_i.pos+1],end=" ")
            for i, c in enumerate(str(t_i.alt)):
                mutable_seq.insert(t_i.pos+i,c)
            print(mutable_seq[t_i.pos-1:t_i.pos+len(t_i.alt)+1],t_i.alt)
            del_target = list(map(lambda x: x+len(t_i.alt) if t_i.pos < x else x,\
                             del_target))
            #print(len(t_i.alt),del_target[-1])
        print(len(mutable_seq), del_num,end=" ")
        #del_target = list(del_target)
        #print(list(del_target))
        for i in del_target[::-1]:
            if mutable_seq[i] == 'D':
                del mutable_seq[i]
            else:
                raise ValueError("not D "+mutable_seq[i-100:i+100])
        print(len(mutable_seq))
        try:
            mutable_seq.remove("D")
            raise AssertionError("something wrong")
        except ValueError:
            continue
        print("fin")
with open("dataset/reference_replace_indel.fa","w") as f:
    SeqIO.write(records, f, "fasta")
