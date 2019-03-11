# --------------------------------------------------------------#
# Sargassum Fungi Analyses - Depends on "Process_Raw_Reads.R"
# Author: Geoffrey Zahn
# --------------------------------------------------------------#


# Packages ####

# source("https://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")


# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(ggpubr)
library(dada2); packageVersion("dada2")
library(stringr)
library(purrr)
library(deseq2)
library(vegan)
library(dplyr)
library(reshape2)
library(ade4)
library(seqinr)

source("plot_bar2.R")
source("./summarize_taxa_Joey711.R")
source("./heatmap_left.R")


# custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44","#b3bf2d",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000")
colorblindr::palette_plot(pal)

# Load data ####
ps = readRDS(file = "./Output/clean_phyloseq_object.RDS")
meta = as.data.frame(sample_data(ps))
tax_table(ps)[,2]
# Subset to Sargassum samples ####
ps <- prune_samples(ps@sam_data$Species == "Sargassum",ps)

# Remove empty ESVs and singletons and doublets ####
ps <- prune_taxa(colSums(otu_table(ps)) > 2, ps)


# Remove non-fungi ####
ps <- subset_taxa(ps, Kingdom == "k__Fungi")

# Remove any newly-empty samples ####
ps <- prune_samples(sample_sums(ps) != 0, ps)
# Sanity check
which(sample_sums(ps) == 0)
summary(taxa_sums(ps))
summary(sample_sums(ps))



# Update metadata names ####
names(ps@sam_data)[5] <- "Structure"
meta = as.data.frame(sample_data(ps))


# Normalize (relative abundance) ####
# Merge Samples by Island
psm = merge_samples(ps, "Island")


# merge by island AND Structure ####
variable1 = as.character(get_variable(ps, "Structure"))
variable2 = as.character(get_variable(ps, "Island"))

sample_data(ps)$NewPastedVar <- mapply(paste0, variable1, variable2, collapse = "_")
psm2 = merge_samples(ps, "NewPastedVar")
  # remember this trick!

# repair merged values
sample_data(psm)$Island <- levels(sample_data(ps)$Island)
sample_data(psm2)
sample_data(psm2)$Island <- levels(sample_data(ps)$Island)
sample_data(psm2)$Structure <- rep(levels(sample_data(ps)$Structure), each = 8)


#sample_data(psm)$Structure <- levels(sample_data(ps)$Structure)

# relative abundance transformation
ps_ra <- transform_sample_counts(psm, function(x) x / sum(x) )
ps_ra2 <- transform_sample_counts(psm2, function(x) x / sum(x) )

as.character((ps_ra2@tax_table[,1]))


# taxonomy table
tax = as.data.frame(ps_ra@tax_table)


# Fix taxonomy name strings ####

# Change "Kingdom" names to remove "k__"
kingdomnames = as.character((ps_ra@tax_table[,1]))
ps_ra@tax_table[,1] <- stringr::str_remove(kingdomnames, "k__")
kingdomnames2 = as.character((ps_ra2@tax_table[,1]))
ps_ra2@tax_table[,1] <- stringr::str_remove(kingdomnames2, "k__")


# Change "Class" names to remove "c__"
classnames = as.character((ps_ra@tax_table[,3]))
ps_ra@tax_table[,3] <- stringr::str_remove(classnames, "c__")
classnames2 = as.character((ps_ra2@tax_table[,3]))
ps_ra2@tax_table[,3] <- stringr::str_remove(classnames2, "c__")

# Change "Order" names to remove "o__"
ordernames = as.character((ps_ra@tax_table[,4]))
ps_ra@tax_table[,4] <- stringr::str_remove(ordernames, "o__")
ordernames2 = as.character((ps_ra2@tax_table[,4]))
ps_ra2@tax_table[,4] <- stringr::str_remove(ordernames2, "o__")

# convert "NA" to NA
unique(as.character(ps_ra@tax_table[,1]))
unique(as.character(ps_ra@tax_table[,3]))
unique(as.character(ps_ra@tax_table[,4]))

ps_ra@tax_table[,1][is.na(ps_ra@tax_table[,1])] <- "Unassigned"
ps_ra@tax_table[,3][is.na(ps_ra@tax_table[,3])] <- "Unassigned"
ps_ra@tax_table[,4][is.na(ps_ra@tax_table[,4])] <- "Unassigned"
ps_ra2@tax_table[,1][is.na(ps_ra2@tax_table[,1])] <- "Unassigned"
ps_ra2@tax_table[,3][is.na(ps_ra2@tax_table[,3])] <- "Unassigned"
ps_ra2@tax_table[,4][is.na(ps_ra2@tax_table[,4])] <- "Unassigned"


ps_ra@tax_table[,1][ps_ra@tax_table[,1] == "NA"] <- "Unassigned"
ps_ra@tax_table[,3][ps_ra@tax_table[,3] == "NA"] <- "Unassigned"
ps_ra@tax_table[,4][ps_ra@tax_table[,4] == "NA"] <- "Unassigned"
ps_ra2@tax_table[,1][ps_ra2@tax_table[,1] == "NA"] <- "Unassigned"
ps_ra2@tax_table[,3][ps_ra2@tax_table[,3] == "NA"] <- "Unassigned"
ps_ra2@tax_table[,4][ps_ra2@tax_table[,4] == "NA"] <- "Unassigned"


# Fix Phylum names
classes = as.character(unique(tax_table(ps_ra)[,3]))
phyla = c("Ascomycota","Unassigned","Ascomycota","Ascomycota","Basidiomycota",
  "Ascomycota","Glomeromycota","Ascomycota","Ascomycota","Basidiomycota")

phyla = plyr::mapvalues(tax_table(ps_ra)[,3],from=classes,to=phyla)
tax_table(ps_ra)[,2] <- phyla
tax_table(ps_ra)[,2]

# Bar plot of relative abundance, split by structure and island ####
plot_bar(ps_ra, fill = "Class",x="Island") +
   geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + theme_bw() +
  theme(axis.title = element_text(size = 24,face = "bold"), axis.text = element_text(size=18),
        legend.title = element_text(size=20), legend.text = element_text(size=16)) + scale_fill_manual(values = pal[c(1:9,11)])  
ggsave("./Output/BarPlot_Fungal_Class_by_Island.png", height = 8, width = 12, dpi=300)

plot_bar(ps_ra, fill = "Phylum",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + theme_bw() +
  theme(axis.title = element_text(size = 24,face = "bold"), axis.text = element_text(size=18),
        legend.title = element_text(size=20), legend.text = element_text(size=16)) + scale_fill_manual(values = pal[c(1:9,11)])  
ggsave("./Output/BarPlot_Fungal_Phylum_by_Island.png", height = 8, width = 12, dpi=300)




# rearrange orders for plot
orders = as.character(unique(tax_table(ps_ra)[,4]))
orders = orders[order(orders)]
orders = orders[c(1:16,18,17)]
plot_bar(ps_ra, fill = "Order",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~ps_ra@sam_data$Structure) +
  labs(x="Site",y="Relative Abundance")+ scale_fill_manual(values = pal[c(1:10,12:17,11:12)], breaks = orders) + theme_bw() 
ggsave("./Output/BarPlot_Fungal_Order_by_Island.png", height = 8, width = 12)

plot_bar(ps_ra, fill = "Phylum",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~ps_ra@sam_data$Structure) +
  labs(x="Site",y="Relative Abundance") + scale_fill_manual(values = pal) + theme_bw() 
ggsave("./Output/BarPlot_Fungal_Phylum_by_Island.png", height = 8, width = 12)

######################################## ---- pick up here
plot_bar(ps_ra2, x = "Island", fill = "Class") +
  geom_bar(stat = "identity") + coord_flip()  + facet_wrap(~Structure) +
  labs(x="Site",y="Relative Abundance") + lims(y=c(0,1)) + theme_bw() + scale_fill_manual(values = pal)
ggsave("./Output/BarPlot_Fungal_Class_by_Island_and_Structure.png", height = 8, width = 12)

plot_bar(ps_ra2, fill = "Order",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + lims(y=c(0,1)) + theme_bw() + scale_fill_manual(values = pal[c(1:10,12:17,11:12)], breaks = orders)
ggsave("./Output/BarPlot_Fungal_Order_by_Island_and_Structure.png", height = 8, width = 12)


# Heatmap Order level taxa relative abundance in all samples, colored by structure ####



# Make matrix
otus = as(otu_table(ps), "matrix")
otus.df = as.data.frame(otus)
orders = str_remove(tax_table(ps)[,4], "o__")
orders[is.na(orders)] <- "NA"
orders[orders == "NA"] <- "Unassigned"
names(otus.df) <- orders
df = as.data.frame(t(otus.df))
df$Order = row.names(df)
order.df = df %>% group_by(Order) %>% summarise_all(sum) %>% as.data.frame()
row.names(order.df) <- order.df$Order
order.df <- order.df %>% select(-Order)
order.df = order.df[which(row.names(order.df) != "Unassigned"),]
# reorder Cols by structure
meta <- as(sample_data(ps),"data.frame")
meta = meta %>% arrange(Structure)
order.df = order.df[,meta$IlluminaName]
order.df = decostand(order.df,"total", MARGIN = 2)

order.matrix = as.matrix(order.df)
hm.pal = c("#ffffff",RColorBrewer::brewer.pal(8,"OrRd"))

# make colors based on Structure
structures = meta$Structure
str.cols = as.character(plyr::mapvalues(structures, from = levels(structures), to=c("#492c24","#397c1c","#3791b2")))

# Heatmap of orders by structure (for individual samples)
png("./Output/Heatmap_of_Order_by_Structure.png", width = 16,height = 16,units = 'in',res = 300)
heatmap_left(t(order.matrix), col = hm.pal, Rowv = NA, Colv=NA, labRow = NA, RowSideColors = str.cols,
             margins = c(50,10), cexCol = 5)
dev.off()


colorblindr::palette_plot(hm.pal)
summary(order.df["Eurotiales",])
# Write sequences to file to BLAST

sink("./Output/rep_set.txt")
(taxa_names(ps))
sink(NULL)

# Summarise taxa ####

order_summary = summarize_taxa(ps,"Order","Structure")
# Remove NAs
order_summary = order_summary[order_summary$Order != "NA",]
# reorder
setorder(order_summary, -meanRA)
# plot
ggplot(order_summary, aes(x= meanRA, y=Order, facet=Structure)) +
  geom_point() + geom_errorbarh(aes(xmin=meanRA-sdRA,xmax=meanRA+sdRA)) + labs(x="Relative Abundance") + # lims(x=c(-.2,1)) +
  facet_wrap(~Structure) +
  theme_bw()
ggsave("./Output/Summarized_taxa_Order_by_Structure.png", height = 8, width = 10, dpi = 300)


class_summary = summarize_taxa(ps,"Class","Structure")
# Remove NAs
class_summary = class_summary[grep(pattern = "c__",class_summary$Class),]
# reorder
setorder(class_summary, -meanRA)

# plot
ggplot(class_summary, aes(x= meanRA, y=Class, facet=Structure)) +
  geom_point() + geom_errorbarh(aes(xmin=meanRA-sdRA,xmax=meanRA+sdRA)) + labs(x="Relative Abundance") +
  facet_wrap(~Structure) +
  theme_bw()
ggsave("./Output/Summarized_taxa_Class_by_Structure.png", height = 8, width = 10, dpi = 300)



# Taxonomy table ... frequency of each taxa (number of samples found in)
write.csv(table(tax_table(ps)),"./Output/TaxTable_no_of_Samples.csv", row.names = FALSE, quote = FALSE)
     
# Separate ESV tables for each island ####
otu = as.data.frame(ps_ra@otu_table)
levels(meta$Island)
rowSums(otu)

Hantu = colSums(otu[meta$Island == "Hantu",]) ; names(Hantu) <- paste(tax$Genus,tax$Species,sep = "_") ; Hantu <- Hantu[names(Hantu) != "NA_NA"]
Jong = colSums(otu[meta$Structure == "Jong",]) ; names(Jong) <- paste(tax$Genus,tax$Species,sep = "_") ; Jong <- Jong[names(Jong) != "NA_NA"]
Kusu = colSums(otu[meta$Structure == "Kusu",]) ; names(Kusu) <- paste(tax$Genus,tax$Species,sep = "_") ; Kusu <- Kusu[names(Kusu) != "NA_NA"]
Semaku = colSums(otu[meta$Structure == "Semaku",]) ; names(Semaku) <- paste(tax$Genus,tax$Species,sep = "_") ; Semaku <- Semaku[names(Semaku) != "NA_NA"]
Sisters = colSums(otu[meta$Structure == "Sisters",]) ; names(Sisters) <- paste(tax$Genus,tax$Species,sep = "_") ; Sisters <- Sisters[names(Sisters) != "NA_NA"]
St_John = colSums(otu[meta$Structure == "St John",]) ; names(St_John) <- paste(tax$Genus,tax$Species,sep = "_") ; St_John <- St_John[names(St_John) != "NA_NA"]
TPL = colSums(otu[meta$Structure == "TPL",]) ; names(TPL) <- paste(tax$Genus,tax$Species,sep = "_") ; TPL <- TPL[names(TPL) != "NA_NA"]
TPT = colSums(otu[meta$Structure == "TPT",]) ; names(TPT) <- paste(tax$Genus,tax$Species,sep = "_") ; TPT <- TPT[names(TPT) != "NA_NA"]

# Top 10 taxa from each island (relative abundance) ####

sink(file = "./Output/Top-Ten_Taxa_From_Each_Island.txt")
for(i in c("Hantu","Jong","Kusu","Semaku","Sisters","St_John","TPL","TPT")){
  print(i)
  print(head(unique(names(get(i))[order(get(i), decreasing = TRUE)]),10))
}
sink(NULL)


df = as.data.frame(otu_table(ps_ra))
names(df) <- c(1:length(names(df)))
max_tax <- apply(df,1,max)

max_tax_id = c()
for(i in c(1:dim(df)[1])){
  max_tax_id[i] <- (which(df[i,] == max_tax[i]))
}

# find top-ten taxa (taxa that are most commonly the most abundant in a sample)
topten_RA_tax = as.numeric(names(head(sort(table(max_tax_id),decreasing = TRUE),10)))
topten_tax_df = as.data.frame(tax[topten_RA_tax,])
topten_tax_df$Rank = (c(1:10))
topten_tax_df$No.Samples = head(sort(table(max_tax_id),decreasing = TRUE),10)
topten_tax_df <- topten_tax_df[,c(1,2,3,4,5,6,7,9,10,8)]

# Write to file
write.csv(topten_tax_df, "./Output/Top-Ten_Most_Abundant_Taxa.csv", row.names = FALSE, quote = FALSE)


# # Look at Aspergillus sydowii (It's a known coral pathogen) ####
# grep(pattern = "sydowii", tax$Species)
# unique(tax$Species)
# sydowii = which(tax$Species == "s__sydowii")
# sydowii_seqs = (tax[sydowii,"seq"])
# ps_ra2 = transform_sample_counts(psm2, function(x) x / sum(x) )
# ps_sydowii = merge_taxa(ps_ra2,sydowii)
# ps_sydowii = subset_taxa(ps_sydowii,Species == "s__sydowii")
# 
# structures = as.factor(sample_data(ps_sydowii)$Structure)
# structures = as.character(plyr::mapvalues(structures, from = "Hold Fast", to="Holdfast"))
# ps_sydowii@sam_data$Structure <- structures
# 
# plot_bar2(na.exclude(ps_sydowii),fill="Island",x="Island") + 
#     facet_wrap(~Structure) +
#   labs(x="Site",fill="Site",y="Relative Abundance") +
#   theme(axis.title = element_text(size = 16), axis.text = element_text(size=12),
#         legend.title = element_text(size=16),
#         strip.text.x = element_text(size = 12, face = "bold")) +
#   scale_fill_manual(values = pal)
# 
# 
# 
#         # panel.background = element_rect(fill = "White"),
#         # panel.grid.major = element_line(size = 0.5, linetype = 'solid',
#         #                                 colour = "Gray"), 
#         # panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
#         #                                 colour = "Gray")) 
# ggsave("./Output/A_sydowii_Relative_Abund_by_Site_and_Structure.png", dpi=300,height = 10,width = 12)
# 
# # in which 6 samples was A. sydowii the most abundant
# sydowii_sample_ids = which(sydowii_seqs == 2)
# sydowii_sample_info = sample_data(ps_ra)[sydowii_sample_ids]
# write.csv(sydowii_sample_info, "./Output/Sydowii_sample_info.csv", row.names = FALSE, quote = FALSE)
# 
# #double-check in otu table for same order
# identical(which(names(as.data.frame(ps@otu_table)) %in% sydowii_seqs), sydowii)
# sydowii_otu = as.data.frame(ps_ra@otu_table)[,sydowii] #using relabund
# 
# # Sydowii otu table and stats ####
# sydowii_otu = otu[,sydowii]
# 
# # Row (sample) mean relative abundance of sydowii
# apply(sydowii_otu,1,mean)
# 
# #mean relabund
# syd_mean = apply(sydowii_otu,1,mean)
# syd_sd = apply(sydowii_otu,1,sd)
# 
# ggplot(mapping = aes(x=meta$Structure[syd_mean>0],y=syd_mean[syd_mean>0])) + 
#   geom_boxplot()
# 
# summary(aov(syd_mean ~ meta$Structure))
# 
# #permanova
# adonis(sydowii_otu[syd_mean>0,]~meta$Structure[syd_mean>0])
# 
#     #no difference in location of A. sydowii
# 

# plot against sample info
ggplot() +
  geom_point(aes(y=ps_ra@sam_data$IlluminaName,x=sydowii_mean_ra,color=ps_ra@sam_data$Island),cex=2) +
  theme(axis.text.y = element_text(size = 4, face = "bold")) +
  labs(x=expression(paste("Mean relative abundance of ", italic("A. sydowii"))),
       y="SampleID", color = "Site")
ggsave("./Output/Relative_Abund_of_A_sydowii_by_Site.png", width = 10, height = 8, dpi=300)  

# Density plots (sqrt of relabund)
ggplot() +
  geom_density(aes(x=sqrt(sydowii_mean_ra),fill=ps_ra@sam_data$Island))
ggsave("./Output/DensityPlot_Aspergillus_sydowii_by_Island.png", dpi = 300)

ggplot() +
  geom_density(aes(x=sqrt(sydowii_mean_ra),fill=ps_ra@sam_data$Structure))
ggsave("./Output/DensityPlot_Aspergillus_sydowii_by_Structure.png", dpi = 300)


#adonis

# Remove empty rows
sydowii_samples = which(rowSums(sydowii_otu) > 0)
sydowii_otu_noempty = sydowii_otu[sydowii_samples,]

# Run permanova
adonis(sydowii_otu_noempty ~ ps_ra@sam_data$Island[sydowii_samples] * ps_ra@sam_data$Structure[sydowii_samples])
  # A. sydowii RA not predicted by structure or site!



# ordinations ####
ps_ra <- transform_sample_counts(ps, function(x) x / sum(x) )
# Phyloseq method trying several ordination methods
NMDS = ordinate(ps_ra, method = "NMDS")
DCA = ordinate(ps_ra, method = "DCA")
CCA = ordinate(ps_ra, method = "CCA")
RDA = ordinate(ps_ra, method = "RDA")
MDS = ordinate(ps_ra, method = "MDS")
PCoA = ordinate(ps_ra, method = "PCoA")

# Plot them  ####
nmdsplot=plot_ordination(ps_ra, NMDS, color = "Island",shape = "Structure") + xlim(c(-0.05,0.05)) + theme_minimal() + labs(color="Site")
ggsave(nmdsplot, "./Output/NMDS_by_Island_and_shape.png", dpi=300)
ggsave("./Output/NMDS_Plot_by_Island_and_Structure.png", dpi=300)

# Top-ten taxa from each structure for each island ####
gen_spp = paste(ps_ra@tax_table[,6], ps_ra@tax_table[,7],sep = "___")
taxa_w_name = grep("NA___NA",gen_spp, invert = TRUE)
gen_spp = gen_spp[taxa_w_name]

ps_namedtaxa = subset_taxa(ps_ra, Genus!="NA"&Species!="NA")
filter_taxa(ps_ra,function(x) x %in% head(sort(colSums(otu_table(ps_ra),decreasing = TRUE),10)))


# 
# for(i in levels(ps_ra@sam_data$Island)) {
# 
#   ps_ras <- subset_samples(ps_ra, Island == i)
#   
#   HoldFast = subset_samples(ps_ras, Structure == "Hold Fast")
#   HoldFast = subset_taxa(HoldFast, Genus != "NA")
#   Leaf = subset_samples(ps_ras, Structure == "Leaf")
#   Leaf = subset_taxa(Leaf, Genus != "NA")
#   Vesicle = subset_samples(ps_ras, Structure == "Vesicle")
#   Vesicle = subset_taxa(Vesicle, Genus != "NA")
#   
#   topten_seqs = names(sort(taxa_sums(HoldFast), TRUE)[1:10])
#   assign(paste0(i,"_holdfast"), tax_table(HoldFast)[topten_seqs])
#   assign(paste0(i,"_holdfast"), 
#          paste(get(paste0(i,"_holdfast"))[,6],get(paste0(i,"_holdfast"))[,7],sep = "___"))
#   sink(file = "./Output/Top_Ten_Taxa_From_Each_Island_and_Structure.txt", append = TRUE)
#     print(paste(i, "HoldFast", sep="___"))
#     print(get(paste0(i,"_holdfast")))
#   sink(NULL)
#   
#   topten_seqs = names(sort(taxa_sums(Leaf), TRUE)[1:10])
#   assign(paste0(i,"_leaf"), tax_table(HoldFast)[topten_seqs])
#   assign(paste0(i,"_leaf"), 
#          paste(get(paste0(i,"_leaf"))[,6],get(paste0(i,"_leaf"))[,7],sep = "___"))
#   sink(file = "./Output/Top_Ten_Taxa_From_Each_Island_and_Structure.txt", append = TRUE)
#     print(paste(i, "Leaf", sep="___"))
#     print(get(paste0(i,"_leaf")))
#   sink(NULL)
#   
#   topten_seqs = names(sort(taxa_sums(Vesicle), TRUE)[1:10])
#   assign(paste0(i,"_vesicle"), tax_table(HoldFast)[topten_seqs])
#   assign(paste0(i,"_vesicle"), 
#          paste(get(paste0(i,"_vesicle"))[,6],get(paste0(i,"_vesicle"))[,7],sep = "___"))
#   sink(file = "./Output/Top_Ten_Taxa_From_Each_Island_and_Structure.txt", append = TRUE)
#     print(paste(i, "Vesicle", sep="___"))
#     print(get(paste0(i,"_vesicle")))
#   sink(NULL)
# }
# 
# plot_bar
# 
# # PermANOVA ####
# sink(file = "./Output/Adonis_Table_Island.txt")
# adonis(otu_table(ps_ra) ~ ps_ra@sam_data$Island * ps_ra@sam_data$Structure)
# sink(NULL)


####====================== Revisions ===============================####

# Venn diagram of ESVs by island and structure

# Heatmap of ESV Class by structure

# Add structure to permanova (done)

# Look closer at structure differences - univariate level

# Check "Unknowns" against bigger database to make sure these are/aren't fungi

# Add error bars to fig. 3

# Double-check A. sydowii sequence against NCBI, etc (BLAST, NJ tree???)


#  left off here -> redoing taxonomy assignment in previous script. Pick up tuesday at beginning of this script.

# Change color palette for all barplots that will go into manuscript






