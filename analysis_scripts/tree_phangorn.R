# need to run decontam on the 4 blank samples
library(decontam)
library(dada2)
library(phyloseq)
library(ggplot2)
library(DECIPHER)
library(phangorn)
t2 <- outtab$newmap

seqs <- getSequences(inputtable)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
dm <- phangorn::dist.ml(phang.align)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
fit = phangorn::pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = phangorn::pml.control(trace = 0))


ps <- phyloseq(tax_table(taxa), sample_data(t2),otu_table(outtab$newtable, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

UFout <- UniFrac(ps,TRUE)
UFoutF <- UniFrac(ps,FALSE)
outT <- cmdscale(UFout)
outF <- cmdscale(UFoutF)
all <- cbind(t2,outT,outF)
colnames(all)[21] <- "Weighted1"
colnames(all)[22] <- "Weighted2"
colnames(all)[23] <- "UnWeighted1"
colnames(all)[24] <- "UnWeighted2"
write.table(all, file = "metadata_with_unifrac.txt", quote = FALSE, sep = "\t")
