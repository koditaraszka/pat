# All TSV files are downloaded from Neale lab, and they are the raw trait 
print('variant')
variant <- read.table('../UKBB/src/variants.tsv', sep = '\t', header = T, quote = "", stringsAsFactors = F)
variant <- subset(variant, select = c("chr", "pos", "ref", "alt", "rsid", "info", "minor_allele"))
colnames(variant) <- c("CHR", "BP", "A1", "A2", "RSID", "INFO", "minor_allele")
variant2 <- variant
variant2$A1 <- variant$A2
variant2$A2 <- variant$A1
variant <- rbind(variant, variant2)
rm(variant2)
variant <- variant[which(variant$CHR %in% 1:22 & variant$A1 %in% c('A','T','C','G') & variant$A2 %in% c('A','T','C','G') & variant$A2 == variant$minor_allele & variant$INFO > 0.8 & variant$A1 != variant$A2),]

#chr	pos	ref	alt	variant	minor_allele	minor_AF	low_confidence_variant	n_complete_samplesAC	ytx	beta	se	tstat	pval
print('bmi')
trait <- read.table('bmi.tsv', sep = '\t', header = T, quote = "", stringsAsFactors = F)
trait <- trait[which(trait$low_confidence_variant == 'false'),]
trait <- subset(trait, select = c("chr", "pos", "ref", "alt", "minor_allele", "minor_AF", "n_complete_samples", "ytx", "beta", "se", "tstat", "pval"))
colnames(trait) <- c("CHR", "BP", "A1", "A2", "minor_allele", "MAF", "N", "YTX", "BETA", "SE", "Z", "P")
trait2 <- trait
trait2$A1 <- trait$A2
trait2$A2 <- trait$A1
trait2$BETA <- -1*(trait2$BETA)
trait2$Z <- -1*(trait2$Z)
trait <- rbind(trait, trait2)
trait <- trait[which(trait$CHR %in% 1:22 & trait$A1 %in% c('A','T','C','G') & trait$A2 %in% c('A','T','C','G') & trait$A2 == trait$minor_allele & trait$MAF > 0.01),]
trait <- merge(trait, variant, by = c("CHR", "BP", "A1", "A2"))
trait <- subset(trait, select = c("RSID", "CHR", "BP", "A1", "A2", "MAF", "N", "YTX", "BETA", "SE", "Z", "P"))
write.table(trait, "bmi.txt", sep = '\t', col.names = T, row.names = F, quote = F)

print('diabp')
trait <- read.table('diabp.tsv', sep = '\t', header = T, quote = "", stringsAsFactors = F)
trait <- trait[which(trait$low_confidence_variant == 'false'),]
trait <- subset(trait, select = c("chr", "pos", "ref", "alt", "minor_allele", "minor_AF", "n_complete_samples", "ytx", "beta", "se", "tstat", "pval"))
colnames(trait) <- c("CHR", "BP", "A1", "A2", "minor_allele", "MAF", "N", "YTX", "BETA", "SE", "Z", "P")
trait2 <- trait
trait2$A1 <- trait$A2
trait2$A2 <- trait$A1
trait2$BETA <- -1*(trait2$BETA)
trait2$Z <- -1*(trait2$Z)
trait <- rbind(trait, trait2)
trait <- trait[which(trait$CHR %in% 1:22 & trait$A1 %in% c('A','T','C','G') & trait$A2 %in% c('A','T','C','G') & trait$A2 == trait$minor_allele & trait$MAF > 0.01),]
trait <- merge(trait, variant, by = c("CHR", "BP", "A1", "A2"))
trait <- subset(trait, select = c("RSID", "CHR", "BP", "A1", "A2", "MAF", "N", "YTX", "BETA", "SE", "Z", "P"))
write.table(trait, "diabp.txt", sep = '\t', col.names = T, row.names = F, quote = F)

print('height')
trait <- read.table('height.tsv', sep = '\t', header = T, quote = "", stringsAsFactors = F)
trait <- trait[which(trait$low_confidence_variant == 'false'),]
trait <- subset(trait, select = c("chr", "pos", "ref", "alt", "minor_allele", "minor_AF", "n_complete_samples", "ytx", "beta", "se", "tstat", "pval"))
colnames(trait) <- c("CHR", "BP", "A1", "A2", "minor_allele", "MAF", "N", "YTX", "BETA", "SE", "Z", "P")
trait2 <- trait
trait2$A1 <- trait$A2
trait2$A2 <- trait$A1
trait2$BETA <- -1*(trait2$BETA)
trait2$Z <- -1*(trait2$Z)
trait <- rbind(trait, trait2)
trait <- trait[which(trait$CHR %in% 1:22 & trait$A1 %in% c('A','T','C','G') & trait$A2 %in% c('A','T','C','G') & trait$A2 == trait$minor_allele & trait$MAF > 0.01),]
trait <- merge(trait, variant, by = c("CHR", "BP", "A1", "A2"))
trait <- subset(trait, select = c("RSID", "CHR", "BP", "A1", "A2", "MAF", "N", "YTX", "BETA", "SE", "Z", "P"))
write.table(trait, "height.txt", sep = '\t', col.names = T, row.names = F, quote = F)

print('sysbp')
trait <- read.table('sysbp.tsv', sep = '\t', header = T, quote = "", stringsAsFactors = F)
trait <- trait[which(trait$low_confidence_variant == 'false'),]
trait <- subset(trait, select = c("chr", "pos", "ref", "alt", "minor_allele", "minor_AF", "n_complete_samples", "ytx", "beta", "se", "tstat", "pval"))
colnames(trait) <- c("CHR", "BP", "A1", "A2", "minor_allele", "MAF", "N", "YTX", "BETA", "SE", "Z", "P")
trait2 <- trait
trait2$A1 <- trait$A2
trait2$A2 <- trait$A1
trait2$BETA <- -1*(trait2$BETA)
trait2$Z <- -1*(trait2$Z)
trait <- rbind(trait, trait2)
trait <- trait[which(trait$CHR %in% 1:22 & trait$A1 %in% c('A','T','C','G') & trait$A2 %in% c('A','T','C','G') & trait$A2 == trait$minor_allele & trait$MAF > 0.01),]
trait <- merge(trait, variant, by = c("CHR", "BP", "A1", "A2"))
trait <- subset(trait, select = c("RSID", "CHR", "BP", "A1", "A2", "MAF", "N", "YTX", "BETA", "SE", "Z", "P"))
write.table(trait, "sysbp.txt", sep = '\t', col.names = T, row.names = F, quote = F)
