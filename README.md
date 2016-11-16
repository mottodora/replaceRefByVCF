editFasta
===============
* Reference配列とVCFファイルを与えると、ReferenceファイルをVCFファイルのSNPやINDELで置換してくれます。

## dependency
* python3
* pandas
* biopython
* PyVCF

## how to use
```
$ python replace_reference_by_vcf.py -r reference.fa --snp snp.vcf --indel indel.vcf -o edited.fa
```

