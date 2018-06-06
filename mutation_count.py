import sys
import pysam


def bases_at_pos(file_name, chromosome, pos):
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total_count = 0
    # BAM file coordinates are 0-based
    pos -= 1
    bam_file = pysam.AlignmentFile(file_name, 'rb')
    for pileup_column in bam_file.pileup(chromosome, pos, pos + 1, truncate=True, stepper='nofilter'):
        # print("\nThere are %s reads overlapping the coordinate at %s:%s" %
        #       (pileup_column.n, chromosome, pileup_column.pos + 1))
        coverage = pileup_column.n
        # Nucleotides at the query position for each read
        for pileup_read in pileup_column.pileups:
            if not pileup_read.is_del and not pileup_read.is_refskip:
                base_qual = pileup_read.alignment.query_alignment_qualities[pileup_read.query_position]
                base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                if base_qual > 0:
                    base_counts[base.upper()] += 1
                    total_count += 1
                # print("Base {} with mapping quality score {}".format(base, base_qual))
    return total_count, base_counts

def nearest_quarter(x):
    return round(x * 4) / 4

def find_mutation_counts(samples):
    results = open("all-samples-23S-count.txt", 'w')
    positions = [(1263703, 'A'), (1621191, 'A'), (1725090, 'A'), (1959084, 'T')]
    ref_seq = "CP001050.1"
    report_all_mapping = "bowtie-mapping-report-all-sorted.bam"
    our_method = "bowtie-remu-sorted.bam"
    other_method = ""
    for sample in samples:
        results.write("{}\n".format(sample))
        file_path = "./{}/mappings/bowtie/".format(sample)
        # The results from multi-mapping resolution
        results.write("Bowtie+REMU\n")
        for pos in positions:
            all_total_count, all_base_counts = bases_at_pos(file_path + report_all_mapping, ref_seq, pos[0])
            total_count, base_counts = bases_at_pos(file_path + our_method, ref_seq, pos[0])
            allele_count = base_counts[pos[1]]
            ratio = allele_count / all_total_count
            results.write("{}\t{}\t{}\t{}\t{}\n".format(pos[0], allele_count, all_total_count, round(ratio, 2), nearest_quarter(ratio)))

        # The results for masking three copies of 23S gene
        results.write("Bowtie-masked\n")
        other_method = "bowtie-mapping-best-match-masked-sorted.bam"
        total_count, base_counts = bases_at_pos(file_path + other_method, ref_seq, positions[3][0])
        allele_count = base_counts['T']
        ref_count = base_counts['C']
        ratio = allele_count / ref_count
        results.write("{}\t{}\t{}\t{}\t{}\n".format(positions[3][0], allele_count, ref_count, round(ratio, 2),
                                                    nearest_quarter(ratio)))
    results.close()
    return True

if __name__ == "__main__":
    if len(sys.argv) > 1:
        find_mutation_counts(sys.argv[1:])
    else:
        print("Error: list of sample names/folders not provided!")

# bases_at_pos("bowtie-remu-sorted.bam", "CP001050.1", 1959084)
