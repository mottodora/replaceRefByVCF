#-*- coding:utf-8 -*-
import re
import argparse

import pandas as pd
from Bio import SeqIO
import vcf


def process_vcf(file):
    vcf_reader = vcf.Reader(open(file, 'r'))
    tmp = {i: [record.CHROM, record.POS, record.REF, record.ALT]\
           for i, record in enumerate(vcf_reader)}
    mutation_list = pd.DataFrame.from_dict(tmp, orient='index')
    mutation_list.rename(columns={0:'chr',1:'pos',2:'ref',3:'alt'}\
                         ,inplace=True)
    mutation_list['alt'] = mutation_list.alt.apply(lambda x: x[0])
    return mutation_list

def main():
    parser = argparse.ArgumentParser(description="Replace reference genome by\
                                     vcf files.")
    parser.add_argument('--snp', type=str, default=None)
    parser.add_argument('--indel', type=str, default=None)
    parser.add_argument('--reference', '-r', type=str, required=True)
    parser.add_argument('--output', '-o', type=str, required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--debug', action='store_true')
    group.add_argument('--no-debug', action='store_false')
    parser.set_defaults(debug=False)
    args = parser.parse_args()

    if args.snp is None and args.indel is None:
        raise ValueError("At least one vcf file is necessary")

    if args.snp is not None:
        snp_list = process_vcf(args.snp)
        snp_list = snp_list[snp_list.alt.apply(lambda x: len(x)==1)]

    if args.indel is not None:
        indel_list = process_vcf(args.indel)
        insertion_list = indel_list[indel_list.alt.apply(lambda x: len(x) > 1)]
        deletion_list = indel_list[indel_list.ref.apply(lambda x: len(x) > 1)]

    records = []
    for record in SeqIO.parse(args.reference, "fasta"):
        print(record.id, end=" ")
        if record.id.startswith('chr'):
            t_chr = re.match("chr(.*)", record.id).group(1)
            mutable_seq = record.seq.tomutable()
            if args.snp is not None:
                t_mut = snp_list[snp_list.chr == t_chr]
                for i, mut in t_mut.iterrows():
                    if mutable_seq[mut.pos-1] == mut.ref:
                        mutable_seq[mut.pos-1] = str(mut.alt)
                    else:
                        raise ValueError("Ref base doesn't match at "\
                                         +str(mut.pos-1))
                print("snp:",len(t_mut), end=" ")
            if args.indel is not None:
                t_ins = insertion_list.loc[insertion_list.chr == t_chr]
                t_del = deletion_list.loc[deletion_list.chr == t_chr]
                del_target = []
                for _,t_d in t_del.iterrows():
                    if mutable_seq[t_d.pos-1] != t_d.ref[0]:
                        raise ValueError("Ref base doesn't match at "\
                                         +str(t_d.pos-1))
                    mutable_seq[t_d.pos:t_d.pos+len(t_d.ref)-1] = \
                            'D'*(len(t_d.ref)-1)
                    del_target\
                            .extend(list(range(t_d.pos,t_d.pos+len(t_d.ref)-1)))
                for _,t_i in t_ins.sort_values('pos', ascending=False).iterrows():
                    if mutable_seq[t_i.pos-1] != t_i.ref:
                        raise ValueError("Ref base doesn't match at "\
                                         +str(t_d.pos-1))
                    for i, c in enumerate(str(t_i.alt)):
                        mutable_seq.insert(t_i.pos+i,c)
                    del_target = list(map(lambda x: x+len(t_i.alt) \
                                          if t_i.pos < x else x, del_target))
                for i in del_target[::-1]:
                    if mutable_seq[i] == 'D':
                        del mutable_seq[i]
                    else:
                        raise ValueError("not D "+mutable_seq[i-100:i+100])
                if args.debug:
                    try:
                        mutable_seq.remove("D")
                        raise AssertionError("something wrong")
                    except ValueError:
                        continue

                print("ins:",len(t_ins), end=" ")
                print("del:",len(t_del), end="")
            record.seq = mutable_seq.toseq()
        records.append(record)
        print('')
    with open(args.output,"w") as f:
        SeqIO.write(records, f, "fasta")
if __name__ == '__main__':
    main()
