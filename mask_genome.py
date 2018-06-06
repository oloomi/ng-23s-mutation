import sys
import ast


def read_genome(genome_file):
    genome_seq = {}  # chrom_name : chrom_seq
    chrom_seq = ""
    with open(genome_file) as ref_genome:
        chrom_name = next(ref_genome).rstrip().split(' ')[0][1:]
        for line in ref_genome:
            if line[0] == ">":
                if chrom_seq:
                    genome_seq[chrom_name] = list(chrom_seq)
                    chrom_seq = ""
                    chrom_name = line.rstrip().split(' ')[0][1:]
            else:
                chrom_seq += line.rstrip()
        genome_seq[chrom_name] = list(chrom_seq)

    return genome_seq


def write_genome(genome_seq, output_file):
    num_chrom = len(genome_seq)
    with open(output_file, 'w') as genome_file:
        for chrom_name, chrom_seq in genome_seq.items():
            # Sequence name
            genome_file.write('>' + chrom_name + '\n')
            # Writing sequence to file, 70 characters per line
            line_width = 70
            length = len(chrom_seq)
            for i in range(length // line_width):
                genome_file.write(''.join(chrom_seq[i * line_width: (i + 1) * line_width]))
                genome_file.write('\n')
            # Writing the last remainder part of genome
            if length % line_width != 0:
                genome_file.write(''.join(chrom_seq[-(length % line_width):]))
            # If there are more chromosomes to be written
            if num_chrom > 1:
                genome_file.write('\n')
                num_chrom -= 1


def mask_genome(genome_file, locations):
    genome_seq = read_genome(genome_file)
    chrom_name = sorted(list(genome_seq.keys())[0])
    for loc in locations[:3]:
        for pos in range(loc[0]-1, loc[1]):
            genome_seq[chrom_name][pos] = 'N'
    write_genome(genome_seq, "{}-masked.fna".format(genome_file[:-4]))
    print("Masked genome written to file.")


if __name__ == "__main__":
    if len(sys.argv) == 3:
        mask_genome(sys.argv[1], ast.literal_eval(sys.argv[2]))
    else:
        print("Error: the path to reference genome file and/or gene locations not provided!")
