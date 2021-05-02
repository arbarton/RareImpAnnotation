import sys
import argparse
import gzip
import math

p=argparse.ArgumentParser()
p.add_argument('--r2', type = str) #plink LD file listing r2 for variants
p.add_argument('--cadd', type = str) #CADD file containing CADD annotations
p.add_argument('--assoc', type = str) #assoc file from output of first script
p.add_argument('--gene', type = str) #annotation file from VEP database
p.add_argument('--dups', type = str) #dup_snps file from first script
p.add_argument('--list', type = str) #list of imputed variants with rs numbers
p.add_argument('--chrom', type = str) #input chromosome number
p.add_argument('--splice', type = str) #spliceAI annotation file
args = p.parse_args()


##create extra dictionaries to store data
genes = {}
chi_stats = {}
caddscores = {}
betas = {}
effects = {}
mafs = {}
p_values = {}
dup_snps = []
lost_snp = {}
win_snp = {}
win_chi = {}
win_num = {}
win_R2 = {}
win_chi_adj = {}
lost_chi_adj = {}
lost_num = {}
lost_chi = {}
lost_R2 = {}
Allele1 = {}
Allele0 = {}
INFO = {}
largest_chi_neighbor = {}
largest_chi_value = {}
largest_chi_R2 = {}
rs_to_name = {}
name_to_rsid = {}

##read in the list of rs numbers that correspond to the variants
with open(args.list) as lines_variant_list:
    for line in lines_variant_list:
            data = line.split()
            if data[0] != '#' and data[0] != 'alternate_ids':
                ids = data[1]
                chrom = int(data[2])
                bp = data[3]
                ref = data[5]
                alt = data[6]
                name = str(chrom) + ":" + str(bp) + ":" + ref + ":" + alt
                rs_to_name[ids] = name
                name_to_rsid[name] = ids

##read in list of variants that are duplicated between the two association files
with open(args.dups) as lines_dups:
    for line in lines_dups:
            data = line.split()
            dup_snps.append(data[0])

all_names = set()
##read in the assoc file and capture the relevant information
with gzip.open(args.assoc,'rb') as lines_assoc:
    for line in lines_assoc:
            data = line.split()
            if line[0] != 'SNP':
                    rs = data[0]
                    rs.strip()
                    if rs in rs_to_name:
                            name = rs_to_name[rs]
                    else:
                            name = rs
                    all_names.add(name)
                    afs = data[6]
                    mafs[name] = afs
                    chi_stats[name] = data[14]
                    p_values[name] = data[15]
                    betas[name] = data[10]
                    Allele1[name] = data[4]
                    Allele0[name] = data[5]
                    INFO[name] = data[7]

ds_ags = {}
ds_als = {}
ds_dgs = {}
ds_dls = {}
splice_genes = {}

##read in the splice AI annotations
with gzip.open(args.splice,'rb') as lines_spliceai:
    for line in lines_spliceai:
        data = line.split()
        ch = data[0]
        chrom = ch.split('r')[1]
        bp = data[2]
        ref = data[3]
        alt = data[4]
        name = str(chrom) + ":" + str(bp) + ":" + ref  + ":" + alt
        if name in all_names:
            splice_annots = data[5]
            annot = splice_annots.split(';')
            gene = annot[0]
            splice_gene = gene.split('=')[1]
            splice_genes[name] = splice_gene
            ag = annot[4]
            ds_ag = ag.split('=')[1]
            ds_ags[name] = ds_ag
            al = annot[5]
            ds_al = al.split('=')[1]
            ds_als[name] = ds_al
            dg = annot[6]
            ds_dg = dg.split('=')[1]
            ds_dgs[name] = ds_dg
            dl = annot[7]
            ds_dl = dl.split('=')[1]
            ds_dls[name] = ds_dl

##read in the CADD file and save values

with gzip.open(args.cadd,'rb') as lines_cadd :
    for line in lines_cadd:
            data = line.split()
            cadds = data[5]
            chrom = data[0]
            bp = data[1]
            ref = data[2]
            alt = data[3]
            name = str(chrom) + ":" + str(bp) + ":" + ref  + ":" + alt
            if name in all_names:
                    caddscores[name] = cadds

##read in the list of genes and save genes and variant effect
with gzip.open(args.gene,'rb') as lines_genes :
    for line in lines_genes:
            data = line.split('\t')
            name = str(data[0]) + ":" + str(data[1]) + ":" + data[3] + ":"  + data[4]
            if name in all_names and name not in genes:
                    gene = data[7]
                    effect = data[6]
                    genes[name] = gene
                    effects[name] = effect

a = "rsID\tChr\tPos\tAllele1\tAllele0\tA1FREQ\tINFO\tCHISQ\tp_value\tBeta\tCADD\tGene\tEffect\tLead\tWin_SNP\tWin_value\tWin_chi\tWin_R2\tWin_Beta\tAdjusted_chi\tLargest_chi_SNP\tLargest_chi\tlargest_chi_r2\tSplice_Gene\tDS_Acceptor_Gain\tDS_Acceptor_Loss\tDS_Donor_Gain\tDS_Donor_Loss"

all_snps = []
not_lead = []
lead = []


##read in the plink ld file and compute the relative values of the variants
with gzip.open(args.r2,'rb') as lines_r2 :
    for line in lines_r2:
            data2 = line.split()
            if data2[0] != 'CHR_A':
                    snp1 = data2[0]
                    snp1.strip
                    if snp1 in rs_to_name:
                        name1 = rs_to_name[snp1]
                    else:
                        name1 = snp1
                    all_snps.append(name1)
                    snp2 = data2[1]
                    snp2.strip
                    if snp2 in rs_to_name:
                        name2 = rs_to_name[snp2]
                    else:
                        name2 = snp2
                    all_snps.append(name2)
                    r2s = (data2[2])
                    r2s.strip
                    R = float(r2s)
                    R2 = R**2
                    if name1 in chi_stats and name2 in chi_stats:
                        if R*float(betas[name1])*float(betas[name2]) > 0:
                        chi1 = float(chi_stats[name1])
                        chi2 = float(chi_stats[name2])
                        if chi1 < chi2:
                            chiL = chi1
                            nameL = name1
                            chiH = chi2
                            nameH = name2
                        else:
                            chiL = chi2
                            nameL = name2
                            chiH = chi1
                            nameH = name1
                        if nameL not in largest_chi_neighbor or chiH > largest_chi_value[nameL]:
                            largest_chi_neighbor[nameL] = nameH
                            largest_chi_value[nameL] = chiH
                            largest_chi_R2[nameL] = R2
                        if nameH not in largest_chi_neighbor or chiL > largest_chi_value[nameH]:
                            largest_chi_neighbor[nameH] = nameL
                            largest_chi_value[nameH] = chiL
                            largest_chi_R2[nameH] = R2
                        if nameH not in win_snp:
                            lead.append(nameH)
                        chi_adj = chiL*((1-math.sqrt(R2*chiH/chiL))**2)
                        if chi_adj < 29.7168:
                            not_lead.append(nameL)
                        if nameL not in lost_snp or chi_adj < lost_chi_adj[nameL]:
                                lost_snp[nameL] = nameH
                                lost_num[nameL] = R2*chiH/chiL
                                lost_chi[nameL] = chiH
                                lost_R2[nameL] = R2
                                lost_chi_adj[nameL] = chi_adj


#Compare winning and losing SNPs and add their annotations in one place
with gzip.open(args.assoc,'rb') as lines_assoc:
    for line in lines_assoc:
            if len(line) > 2:
                    data = line.split()
                    if len(data) > 1 :
                            if data[1] == args.chrom:
                                    ids  = data[0]
                                    if ids not in dup_snps:
                                            if ids in rs_to_name:
                                                    name = rs_to_name[ids]
                                            else:
                                                    name = ids
                                            if name not in not_lead:
                                                win_status = "won"
                                            else:
                                                win_status = "lost"
                                            if name in lost_snp:
                                                num_lost = lost_num[name]
                                                snp_temp = lost_snp[name]
                                                if snp_temp in name_to_rsid:
                                                    snp_lost = name_to_rsid[snp_temp]
                                                else:
                                                    snp_lost = snp_temp
                                                chi_lost = lost_chi[name]
                                                R2_lost = lost_R2[name]
                                                lose_beta = betas[snp_temp]
                                                chi_adj = lost_chi_adj[name]
                                            else:
                                                num_lost = "N/A"
                                                snp_lost = "N/A"
                                                chi_lost = "N/A"
                                                lose_beta = "N/A"
                                                chi_adj = "N/A"
                                                R2_lost = "N/A"
                                            if name in largest_chi_neighbor: #annotate variants with the name of the SNP they lost too or significant neighbors
                                                name_large_chi_temp = largest_chi_neighbor[name]
                                                value_large_chi = largest_chi_value[name]
                                                R2_large_chi = largest_chi_R2[name]
                                                if name_large_chi_temp in name_to_rsid:
                                                    name_large_chi = name_to_rsid[name_large_chi_temp]
                                                else:
                                                    name_large_chi = name_large_chi_temp
                                            else:
                                                name_large_chi = "N/A"
                                                value_large_chi = "N/A"
                                                R2_large_chi = "N/A"
                                            af = mafs[name] #Add additional information about variant
                                            beta = betas[name]
                                            info = INFO[name]
                                            A1 = Allele1[name]
                                            A0 = Allele0[name]
                                            if name in caddscores: #Add CADD score annotation
                                                    cd = caddscores[name]
                                            else:
                                                    cd = 'LOW'
                                            if name in ds_ags: #Add Splice AI annotations if present
                                                ds_ag = ds_ags[name]
                                                ds_al = ds_als[name]
                                                ds_dg = ds_dgs[name]
                                                ds_dl = ds_dls[name]
                                                splice_gene = splice_genes[name]
                                            else:
                                                ds_ag = "N/A"
                                                ds_al = "N/A"
                                                ds_dg = "N/A"
                                                ds_dl = "N/A"
                                                splice_gene = "N/A"
                                            if name in genes: #Add gene name if present in VEP file
                                                    gene_id = genes[name]
                                                    effect_id = effects[name]
                                            else:
                                                    gene_id = 'NONE'
                                                    effect_id = "N/A"
                                            pvalue = p_values[name]
                                            sitestate = chi_stats[name]
                                            a += "\n" + str(ids) + '\t' + str(data[1]) + '\t' + str(data[2]) + '\t' + A1 + "\t" + A0 + "\t" + \
                                            str(af) + "\t" + str(info) + '\t' +  str(sitestate) + '\t' + str(pvalue)  + "\t" + beta + "\t" +  str(cd) + "\t" + \
                                            str(gene_id) + "\t" + str(effect_id) + "\t" + win_status + "\t" + snp_lost + "\t" + str(num_lost) + "\t" + str(chi_lost) \
                                            + "\t" + str(R2_lost) + "\t" + str(lose_beta) + "\t" + str(chi_adj) + "\t" + str(name_large_chi) + "\t" + str(value_large_chi) + "\t" + str(R2_large_chi) \
                                            + "\t" + str(splice_gene) + "\t" + str(ds_ag) + "\t" + str(ds_al) + "\t" + str(ds_dg) + "\t"  + ds_dl

print(a)
