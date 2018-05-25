#!/usr/bin/python3

import sys
import ast
import copy
from vcf_file import *


def variant_evaluation(gene_loci):
    variant_caller_lst = [("Freebayes", "freebayes")]
    file_path = "./variants/"
    vcf_files_names = [["Bowtie2 best-match", "bowtie-mapping-best-match-sorted"],
                       ["Bowtie2 report-all", "bowtie-mapping-report-all-sorted"],
                       ["Bowtie2 + MMR", "bowtie-mmr-sorted"],
                       ["Bowtie2 + REMU", "bowtie-remu-sorted"]]

    evaluation_results = open("./results/variants-comparison-freebayes.txt", 'w')

    for variant_caller in variant_caller_lst:
        vcf_files = copy.deepcopy(vcf_files_names)
        for i in range(len(vcf_files)):
            vcf_files[i][1] = file_path + vcf_files[i][1] + "-variants-{}.vcf".format(variant_caller[1])

        for vcf_file in vcf_files:
            method_name = vcf_file[0]
            vcf_file_name = vcf_file[1]
            # If the extension is missing
            if vcf_file_name[-4:].lower() != ".vcf":
                vcf_file_name += ".vcf"
            # Reading variants for this mapping
            called_variants = read_vcf_file(vcf_file_name, gene_loci)
            if not called_variants:
                print("No variants called in {} !".format(vcf_file))
                evaluation_results.write("No variants called in {} !\n".format(vcf_file))
                continue
            called_variants = ["\t".join([str(e) for e in v]) for v in called_variants]
            evaluation_results.write("{}\n{}\n".format(method_name, "\n".join(called_variants)))

    evaluation_results.close()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        variant_evaluation(ast.literal_eval(sys.argv[1]))
    else:
        print("Error: list of genes locations not provided!")

