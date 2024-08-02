# DECONTAM

# need to run decontam on the 4 blank samples
library(decontam)
library(phyloseq)
library(ggplot2)
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
t2 <- outtab$newmap
t2$ContaminateCheck <- "Real"
t2$ContaminateCheck[66:71] <- "Blank"

ps <- phyloseq(otu_table(outtab$newtable, taxa_are_rows=FALSE), sample_data(t2), tax_table(taxa))
#ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 30958 taxa and 113 samples ]
#sample_data() Sample Data:       [ 113 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 30958 taxa by 7 taxonomic ranks ]

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize ),]
df$Index <- seq(nrow(df))
p1 <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=ContaminateCheck)) + geom_point()
sample_data(ps)$is.neg <- sample_data(ps)$ContaminateCheck == "Blank"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
ggsave(p1, file = "library_size_1.png", dpi = 800, height = 6, width = 7, units = "in")
table(contamdf.prev$contaminant)

# DECONTAM removed 6 ASV sequences present in blank samples > real samples
# FALSE  TRUE 
# 30952     6 
# 

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$ContaminateCheck == "Blank", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$ContaminateCheck == "Real", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
p2<- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave(p2, file = "library_size_2.png", dpi = 800, height = 6, width = 7, units = "in")

ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
#ps.noncontam
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 30952 taxa and 113 samples ]
#sample_data() Sample Data:       [ 113 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 30952 taxa by 7 taxonomic ranks ]
#brayWmeta <- diversity(outtab$newmap,outtab$newtable)

t1b <- as.data.frame(otu_table(ps.noncontam))
t2b <- as.data.frame(sample_data(ps.noncontam))
t3b <- as.data.frame(tax_table(ps.noncontam))
# remove added columns to metadata
t2b <- t2b[,1:4]
outtab <- Create_Tables(t1b,t2b,t3b)

