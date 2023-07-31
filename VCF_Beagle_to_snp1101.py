import pandas as pd 
import numpy as np
import sys


snp = pd.read_csv("/gpfs/fs0/scratch/y/ymiar/milad/P/vcf/hwe_QQ_imp_vcf_plink/myvcf.vcf",delimiter=r"\s+",header=0,skiprows=9)
snp2= snp.drop(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'],axis=1)
snp3=snp2.transpose()
snp3.iloc[:, 0:] = snp3.iloc[:, 0:].apply(lambda s: s.str[:3])

snp4=snp3.replace(['0|0','0|1','1|0','1|1','.|.'],['0','1','1','2','5'])

snp5 = snp4[snp4.columns[0:]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
snp6= pd.DataFrame(snp5).reset_index(drop=False)

snp6.columns=["ID","calls..."]
snp7=snp6.reset_index(drop=True)
snp7.to_csv('geno.txt', sep='\t',index=False)


MAP=pd.DataFrame(snp[['ID','#CHROM','POS']])
MAP.columns=["SNPID","Chr",'Pos']
MAP.to_csv('Map.txt', sep='\t',index=False)