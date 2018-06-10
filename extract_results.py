
with open('samples-1-26-23S-count.txt') as infile:
    with open('samples-1-26-23S-count-summary.txt', 'w') as outfile:
        # For 26 samples
        for s in range(26):
            # Sample name
            outfile.write(infile.readline().rstrip() + '\t')
            next(infile)
            sample_result = infile.readline().rstrip().split(sep='\t')
            for i in range(3):
                curr_result = infile.readline().rstrip().split(sep='\t')
                if float(curr_result[3]) > float(sample_result[3]):
                    sample_result = curr_result
            outfile.write('\t'.join(sample_result) + '\t')
            next(infile)
            sample_result = infile.readline().rstrip().split(sep='\t')
            outfile.write('\t'.join(sample_result) + '\n')