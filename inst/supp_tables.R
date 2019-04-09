# get supp tables

#no Norm
vcf <- VariantAnnotation::readVcf(vcfPath)
seqs = as.character(seqnames(vcf)@values)[-length(as.character(seqnames(vcf)@values))]
profileX <- vcf2mut_cat(vcf, BSgenome.Hsapiens.1000genomes.hs37d5, seqs)
statsX <- statistics(profileX)

genomeCounts <- get_context_freq(BSgenome.Hsapiens.UCSC.hg19)
exomeCounts <- get_context_freq(BSgenome.Hsapiens.UCSC.hg19, exomePath)

# genome Norm
profileXgN <- normalize(profileX, exomeCounts, genomeCounts)
statsXgN <- statistics(profileXgN, m = colSums(profileX))

# exome Norm
exSigs <- normalize(cosmicSigs, genomeCounts, exomeCounts)
statsXeN = statistics(profileX, exSigs)

# uniNorm
uni <- rep(1,32)
uniSigs <- normalize(cosmicSigs, genomeCounts, uni)
profileXuN <- normalize(profileX[-length(profileX)], v4contexts, uni)
statsXuN <- statistics(profileXuN, P=uniSigs, m = colSums(profileX))


write.table(profile99_4uN, 
            file='~/Documents/sigPaper/signaturepackage/Paper/supp_tables/99/99_4profile_uniNorm.tsv', 
            quote=FALSE, 
            sep='\t', 
            col.names = NA)

stats99_4uN <- statistics(profile99_4uN, P=uniSigs , m=colSums(as.matrix(profile99_4[-length(profile99_4)])))
write.table(stats99_4uN, 
            file='~/Documents/sigPaper/signaturepackage/Paper/supp_tables/99/99_4stats_uniNorm.tsv', 
            quote=FALSE, 
            sep='\t', 
            col.names = NA)

write.table(uniSigs, 
            file='~/Documents/sigPaper/signaturepackage/Paper/supp_tables/cosmicSigs_uniNorm.tsv', 
            quote=FALSE, 
            sep='\t', 
            col.names = NA)


# comparing figures were done with genome norm profiles
all69 = cbind(profile69_1gN, profile69_2gN, profile69_4gN, profile69_5gN)
set.seed(123)
statistics(all69, m=colSums(as.matrix(profile99_4[-length(profile99_4)])))

set.seed(123)
ortho = summarize_exposures(create_mut_catalogues(1,1000, sig_set=c(7,13,21,24,28), c_exposure = c(0.2,0.2,0.2,0.2,0.2))[['catalogues']][,1])
write.table(ortho, 
            file='~/Documents/sigPaper/signaturepackage/Paper/supp_tables/set1.tsv', 
            quote=FALSE, 
            sep='\t', 
            col.names = NA)

set.seed(123)
set2 = summarize_exposures(create_mut_catalogues(1,1000, sig_set=c(3,5,8,16,25), c_exposure = c(0.2,0.2,0.2,0.2,0.2))[['catalogues']][,1])
write.table(set2, 
            file='~/Documents/sigPaper/signaturepackage/Paper/supp_tables/set2.tsv', 
            quote=FALSE, 
            sep='\t', 
            col.names = NA)
