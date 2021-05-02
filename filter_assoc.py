import argparse
import sys
import gzip

p=argparse.ArgumentParser()
p.add_argument('--pheno', type = str) #list of phenotypes to be tested
p.add_argument('--path', type = str) #path to files with significant associations from BOLT-LMM, formatted as below
args = p.parse_args()

#filenames of files with significant BOLT-LMM associations from the WES file and imputation files
file1name = args.path + args.pheno + ".imp_v3.stats.gz" ##or append with HRC imputation variants
file2name = args.path + args.pheno + ".WES.stats.gz" ## or append with WES imputation variants
outfilename = args.path + "/sig." + args.pheno + ".FILTER.stats.gz"

#open files with significant variants
out = gzip.open(outfilename,'wb')
lines = gzip.open(file1name, 'rb').readlines()
lines2 = gzip.open(file2name, 'rb').readlines()

#make list of variants in both files
WES_variants = {}
for line in lines2:
        data = line.split()
        chrom  = data[1]
        bp = data[2]
        pos = (chrom,bp)
        alt = data[5]
        WES_variants[pos] = alt
        
for line in lines:
        data = line.split()
        chrom  = data[1]
        bp = data[2]
        pos = (chrom,bp)
        alt = data[5]
        ref = data[4]
        if pos in WES_variants:
                if WES_variants[pos] != alt and  WES_variants[pos] != ref:
                        out.write(line)
        else:
                out.write(line)

out.close()
