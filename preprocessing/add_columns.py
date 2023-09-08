import argparse
import math
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom

# write normal GT to allele count (if 2 ref, than 2,0 etc.)
def get_normal_genotype(df: pd.DataFrame) -> str:
    gt = df["sample_control"].split(":")[0]
    if gt == "1/0" or gt == "0/1":
        return "1,1"
    if gt == "0/0":
        return "2,0"
    if gt == "1/1":
        return "0,2"

# write CNA from segment file, so that it becomes clear, which copy is altered / lost / gained. E.g. genotype segment file = 2:1, you don't know which allele is carried 2 times and which 1 time
# genotype in the segment file means sth different! It means how many copies you have of each allele! 1:1 does not mean 2 times ALT, but one copy from both parents. There does not have to be any variant!
def get_tumor_genotype_with_qualities(df: pd.DataFrame, purity: float) -> tuple[str, float]:
    gt = df["genotype"]
    
    # exception to eliminate errors
    if (gt != gt):
        print("gt is NA")
        genotype = "no_gt"
        quality_score = float('nan')
    # meaning that there are multiple clones
    elif "sub" in gt:
        genotype = "sub"
        quality_score = float('nan')
    # main calculation done here
    else:
        x = int(gt.split(":")[0])
        y = int(gt.split(":")[1])
        # if both copies (one of each parent) are present, the order does not matter
        if x == y:
            genotype = f"{x},{y}"
            quality_score = 0
        else:
            try:
                # maximum likelihood with n = NR, k = NV
                n = df["sample_tumor"].split(":")[-2]
                k = df["sample_tumor"].split(":")[-1]

                # multiallelic variants are written to no_gt
                if "," in n or "," in k:
                    print("no_gt because of multiallelic variant")
                    genotype = "no_gt"
                    quality_score = float('nan')
                
                # calculation of variant allele fraction
                else:
                    n = int(n)
                    k = int(k)
                    tcn_tumor = float(df["TCN"])
                    
                    # special case for purity = 1, because log(0) cannot be calculated
                    if purity == 1:
                        purity == 0.9999
    
                    # case 1: x,y (ref, alt)
                    vaf_1 = ((y * purity) + (1 - purity)) / ((tcn_tumor * purity) + (2 * (1 - purity)))
                    prob_1 = binom.pmf(k=k, n=n, p=vaf_1)
                    
                    # case 2: y,x (ref, alt)
                    vaf_2 = ((x * purity) + (1 - purity)) / ((tcn_tumor * purity) + (2 * (1 - purity)))
                    prob_2 = binom.pmf(k=k, n=n, p=vaf_2)

                    # calculating raw phred-scaled likelihoods
                    pl_1 = -10 * math.log10(prob_1)
                    pl_2 = -10 * math.log10(prob_2)

                    # calculating quality score (capped at 99)
                    quality_score = min(99, abs(pl_1 - pl_2))

                    # print correct order of alleles, so that allele count matches to sample_tumor GT (SNP_GT) --> important for eQTLs later on
                    if prob_1 > prob_2:
                        ref = x
                        alt = y
                    else:
                        ref = y
                        alt = x
                    
                    genotype = f"{ref},{alt}"
                
            except ValueError:
                print(f"no_gt because of ValueError")
                genotype = "no_gt"
                quality_score = float('nan')
    
    return (genotype, quality_score)

def get_total_read_count_tumor(df: pd.DataFrame) -> str:
    trc = df["sample_tumor"].split(":")[-2]
    if "," in trc:
        return "NA"
    return trc

def get_total_read_count_normal(df: pd.DataFrame) -> str:
    trc = df["sample_control"].split(":")[-2]
    if "," in trc:
        return "NA"
    return trc               

def main(header_file: str, usecols_file: str, snv_file_path: str, output_path: str):
    
    # read required information
    with open(header_file) as file:
        header_names = file.readlines()

    header_names = [h.strip("\n") for h in header_names]

    with open(usecols_file) as file:
        usecols_list = file.readlines()

    usecols_list = [c.strip("\n") for c in usecols_list]

    df = pd.read_csv(snv_file_path, sep='\t', names=header_names, usecols=usecols_list, low_memory=False, compression="gzip")
    df_info = pd.read_csv("../Identitymatrix_merged_additionalInfo.tsv", delimiter="\t")
    PID = df.PID.values[0]    
    purity = df_info.set_index("PID").loc[str(PID), "purity"]
    print("PID:" + str(PID) + ", purity:" + str(purity))

    # exit if no purity was assigned
    if np.isnan(purity):
        try:
            sys.exit("No purity for this PID!\n")
        except SystemExit as message:
            print(message)
    
    # process columns
    df["normal_genotype"] = df.apply(get_normal_genotype, axis=1)
    df["tumor_genotype"], df["quality_score"] = zip(*df.apply(lambda x: get_tumor_genotype_with_qualities(x, purity), axis=1))
    df["reads_normal"] = df.apply(get_total_read_count_normal, axis=1)
    df["reads_tumor"] = df.apply(get_total_read_count_tumor, axis=1)

    # add dummy position columns in col 3
    df["position_dummy"] = df["POS"]
    df["POS"] = df["POS"]-1
    df_sorted = df[["#CHROM", "POS", "position_dummy"] + usecols_list[2:] + ["normal_genotype", "tumor_genotype", "quality_score", "reads_normal", "reads_tumor"]]

    # write df to csv
    df_sorted.to_csv(output_path, sep='\t', na_rep="NA", index=False, compression="gzip")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="path to SNV file to add columns to", dest="file_path")
    parser.add_argument("-o", "--output", help="path to file to write output to", dest="out_path")
    args = parser.parse_args()

    file_path = args.file_path
    out_path = args.out_path

    main("header_names.txt", "snv_usecols.txt", file_path, out_path)

    

