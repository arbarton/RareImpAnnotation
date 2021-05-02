#!/bin/bash

PATH1=$1
FILE1=$2 #WES imputed significant variants output from BOLT-LMM
FILE2=$3 #file prefix of imputed WES variants
FILE3=$4 #file prefix of imputed HRC variants
FILE4=$5 #file of individuals to keep for this phenotype
PHENO=$6 #phenotype to be tested
THREADS=$7 #number of threads to use
mkdir -p $PATH1

python filter_assoc.py --pheno "$PHENO" --path "$PATH1"

#Format list of SNPs produced by python program
zcat "$PATH1"/sig."$PHENO".FILTER.stats.gz | awk '{print $1}' > "$PATH1"/sig."$PHENO".FILTER.stats.gz.assoc.snps

zcat "$PATH1"/ | awk '{print $1}' > "$PATH1"/sig."$PHENO".WES.stats.gz.assoc.snps

for CHR in {1..22}
do

### either the old or new bgenix should work for extracting variants
BGENIX= #Insert path to BGENIX
CATBGEN= #Insert path to CATBGEN
### plink2 binaries are available for simple download
PLINK2= #Insert path to plink2


### use bgenix to extract significant variants
$BGENIX \
    -g "$FILE2".bgen  \
    -incl-rsids  "$PATH1"/sig."$PHENO".WES.stats.gz.assoc.snps \
    > "$PATH1"/"$PHENO"_WES_chr"$CHR"_v12.bgen

##list variants to check if empty
$BGENIX -g "$PATH1"/"$PHENO"_WES_chr"$CHR"_v12.bgen -index -clobber > "$PATH1"/"$PHENO"_WES_chr"$CHR"_v12.bgen.bgi
$BGENIX -g "$PATH1"/"$PHENO"_WES_chr"$CHR"_v12.bgen -list > "$PATH1"/"$PHENO"_WES_chr"$CHR"_v12.list

### use plink2 to convert bgen v1.2 to bgen v1.1 and exclude duplicate rsIDs
$PLINK2 --bgen "$PATH1"/${PHENO}_WES_chr"$CHR"_v12.bgen \
    --sample "$FILE2".sample \
    --export bgen-1.1 \
    --threads $THREADS \
    --memory 8000 \
    --out "$PATH1"/"$PHENO"_WES_chr"$CHR"_v11_no_dup_rsIDs

### fix missing sex in v11 sample file
sed -i 's/NA/0/g' "$PATH1"/"$PHENO"_WES_chr"$CHR"_v11_no_dup_rsIDs.sample

##do again for other variants
$BGENIX \
    -g "FILE3".bgen  \
    -incl-rsids   "$PATH1"/sig."$PHENO".FILTER.stats.gz.assoc.snps\
    > "$PATH1"/"$PHENO"_chr"$CHR"_v12.bgen


### list variants to find duplicate rsIDs
$BGENIX -g "$PATH1"/"$PHENO"_chr"$CHR"_v12.bgen -index -clobber > "$PATH1"/"$PHENO"_chr"$CHR"_v12.bgen.bgi
$BGENIX -g "$PATH1"/"$PHENO"_chr"$CHR"_v12.bgen -list > "$PATH1"/"$PHENO"_chr"$CHR"_v12.list
cut -f2 "$PATH1"/"$PHENO"_chr"$CHR"_v12.list | uniq -c | awk '$1>1 {print $2}' > "$PATH1"/"$PHENO"_chr"$CHR"_v12_dup_rsIDs.txt

### use plink2 to convert bgen v1.2 to bgen v1.1 and exclude duplicate rsIDs
$PLINK2 --bgen "$PATH1"/"$PHENO"_chr"$CHR"_v12.bgen \
    --sample "FILE3".sample \
    --export bgen-1.1 \
    --exclude "$PATH1"/"$PHENO"_chr"$CHR"_v12_dup_rsIDs.txt \
    --threads $THREADS \
    --memory 8000 \
    --out "$PATH1"/"$PHENO"_chr"$CHR"_v11_no_dup_rsIDs

### fix missing sex in v11 sample file
sed -i 's/NA/0/g' "$PATH1"/"$PHENO"_chr"$CHR"_v11_no_dup_rsIDs.sample

lines_imp=$( wc -l "$PATH1"/"$PHENO"_chr"$CHR"_v12.list | awk '{print $1}' )
lines_WES=$( wc -l "$PATH1"/"$PHENO"_WES_chr"$CHR"_v12.list | awk '{print $1}' )

if [$lines_imp -eq 3]
then
    PREFIX="$PATH1"/"$PHENO"_chr"$CHR"_v11_no_dup_rsIDs


elif [$lines_WES -eq 3]
then
    PREFIX="$PATH1"/"$PHENO"_WES_chr"$CHR"_v11_no_dup_rsIDs

else
    $CATBGEN  -g "$PATH1"/"$PHENO"_chr"$CHR"_v11_no_dup_rsIDs.bgen "$PATH1"/"$PHENO"_WES_chr"$CHR"_v11_no_dup_rsIDs.bgen -og  "$PATH1"/"$PHENO"_chr"$CHR"_merged.bgen -clobber
    PREFIX="$PATH1"/"$PHENO"_chr"$CHR"_merged
    cp "$PATH1"/"$PHENO"_WES_chr"$CHR"_v11_no_dup_rsIDs.sample $PREFIX.sample
fi

plink --bgen $PREFIX.bgen \
    --sample $PREFIX.sample \
    --keep-allele-order \
    --hard-call-threshold 0.25 \
    --make-bed \
    --memory 8000 \
    --out $PREFIX

plink --bfile $PREFIX \
    --keep "FILE4" \
    --r \
    --memory 8000 \
    --threads $THREADS \
    --ld-window-kb 3000 \
    --ld-window 1000000 \
    --keep-allele-order \
    --out "$PATH1"/"$PHENO"_chr"$CHR"_v11

awk '$7*$7>=1e-4 {print $3,$6,$7}' "$PATH1"/"$PHENO"_chr"$CHR"_v11.ld | gzip > "$PATH1"/"$PHENO"_chr"$CHR"_v11.ld.gz

rm "$PATH1"/"$PHENO"_chr"$CHR"_v11.ld

rm "$PATH1"/"$PHENO"_chr"$CHR"_merged.{bgen,sample}
rm $PREFIX.{bed,bim,fam}

rm "$PATH1"/"$PHENO"{_WES,}_chr"$CHR"_v11_no_dup_rsIDs.{bgen,sample}
rm "$PATH1"/"$PHENO"{_WES,}_chr"$CHR"_v12.bgen{.bgi,}
rm "$PATH1"/"$PHENO"_chr"$CHR"_{v11,merged}.nosex
done

zcat "FILE1" "$PATH1"/sig."$PHENO".FILTER.stats.gz | bgzip > "$PATH1"/sig."$PHENO".ALL.stats.gz

echo Successfully Completed
