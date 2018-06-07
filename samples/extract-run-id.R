selected.strains <- read.csv('strain-allele-count.csv', colClasses = c(rep('character', 3), 'numeric'))

# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERP001405
data1 <- read.table('SraRunTable1.txt', header = TRUE, sep = '\t')
selected.data1 <- subset(data1, strain %in% selected.strains$Strain)

# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERP000144
data2 <- read.table('SraRunTable2.txt', header = TRUE, sep = '\t')
selected.data2 <- subset(data2, strain %in% selected.strains$Strain)

union(selected.data1$strain, selected.data2$strain)

full.data1 <- merge(selected.strains, selected.data1, by.x = 'Strain', by.y = 'strain')
full.data2 <- merge(selected.strains, selected.data2, by.x = 'Strain', by.y = 'strain')

fields <- c('Strain', 'Run', 'Estimate', 'Actual', 'MIC')
full.selected.strains <- rbind(full.data1[, fields], full.data2[, fields])

final.strain.data <- full.selected.strains[order(full.selected.strains$MIC, full.selected.strains$Actual), ]
write.csv(final.strain.data, file = 'johnson-et-al-data.csv')

write.table(t(final.strain.data), file = 'strain-run.txt', 
            row.names = FALSE, col.names = FALSE, sep = ' ', fileEncoding = 'UTF-8')

rep(100, 26)
