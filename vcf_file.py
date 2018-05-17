def read_vcf_file(vcf_file_name, gene_loci):
    """
    Returns a list of (position, variant) from a VCF file
    """
    variants = []
    with open(vcf_file_name) as vcf_file:
        for line in vcf_file:
            # Skip header lines
            if line[0] == "#":
                continue
            # 0.CHROM   1.POS   2.ID    3.REF   4.ALT	5.QUAL	6.FILTER	7.INFO  8.FORMAT
            fields = line.rstrip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = float(fields[5])
            # AO: Count of full observations of this alternate haplotype
            # DP: Total read depth at the locus
            info = dict(item.split("=") for item in fields[7].split(";"))
            if qual > 0:
                for locus in gene_loci:
                    # locus = int(row[1]) - 1956487 + 14
                    if locus[0] <= pos <= locus[1]:
                        start_dist = pos - (locus[0] - 1) + 14
                        end_dist = locus[1] - (pos - 1) + 14
                        variants.append((chrom, pos, ref, alt, start_dist, end_dist, info['AO'], info['DP']))
    return variants
