# --------------------------------------------------------------#
# Sargassum Fungi Analyses - Depends on "Process_Raw_Reads.R"
# Author: Geoffrey Zahn
# --------------------------------------------------------------#


# Packages ####
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
library(VennDiagram)

# functions
source("./R/plot_bar2.R")
source("./R/summarize_taxa_Joey711.R")
source("./R/heatmap_left.R")


# custom color palette
pal = c("#c4a113","#c1593c","#643d91","#820616","#477887","#688e52",
        "#12aa91","#705f36","#8997b2","#753c2b","#3c3e44","#b3bf2d",
        "#82b2a4","#894e7d","#a17fc1","#262a8e","#abb5b5","#000000")
colorblindr::palette_plot(pal[11])

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


# Merge Samples by Island ####
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


# relative abundance transformation
ps_ra <- transform_sample_counts(psm, function(x) x / sum(x) )
ps_ra2 <- transform_sample_counts(psm2, function(x) x / sum(x) )


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

# Change "Family" names to remove "f__"
familynames = as.character((ps_ra@tax_table[,5]))
ps_ra@tax_table[,5] <- stringr::str_remove(familynames, "f__")
familynames2 = as.character((ps_ra2@tax_table[,5]))
ps_ra2@tax_table[,5] <- stringr::str_remove(familynames2, "f__")

# Change "Genus" names to remove "g__"
genusnames = as.character((ps_ra@tax_table[,6]))
ps_ra@tax_table[,6] <- stringr::str_remove(genusnames, "g__")
genusnames2 = as.character((ps_ra2@tax_table[,6]))
ps_ra2@tax_table[,6] <- stringr::str_remove(genusnames2, "g__")

# convert "NA" to NA
unique(as.character(ps_ra@tax_table[,1]))
unique(as.character(ps_ra@tax_table[,3]))
unique(as.character(ps_ra@tax_table[,4]))
unique(as.character(ps_ra@tax_table[,5]))
unique(as.character(ps_ra@tax_table[,6]))

ps_ra@tax_table[,1][is.na(ps_ra@tax_table[,1])] <- "Unassigned"
ps_ra@tax_table[,3][is.na(ps_ra@tax_table[,3])] <- "Unassigned"
ps_ra@tax_table[,4][is.na(ps_ra@tax_table[,4])] <- "Unassigned"
ps_ra@tax_table[,5][is.na(ps_ra@tax_table[,5])] <- "Unassigned"
ps_ra@tax_table[,6][is.na(ps_ra@tax_table[,6])] <- "Unassigned"

ps_ra2@tax_table[,1][is.na(ps_ra2@tax_table[,1])] <- "Unassigned"
ps_ra2@tax_table[,3][is.na(ps_ra2@tax_table[,3])] <- "Unassigned"
ps_ra2@tax_table[,4][is.na(ps_ra2@tax_table[,4])] <- "Unassigned"
ps_ra2@tax_table[,5][is.na(ps_ra2@tax_table[,5])] <- "Unassigned"
ps_ra2@tax_table[,6][is.na(ps_ra2@tax_table[,6])] <- "Unassigned"


ps_ra@tax_table[,1][ps_ra@tax_table[,1] == "NA"] <- "Unassigned"
ps_ra@tax_table[,3][ps_ra@tax_table[,3] == "NA"] <- "Unassigned"
ps_ra@tax_table[,4][ps_ra@tax_table[,4] == "NA"] <- "Unassigned"
ps_ra@tax_table[,5][ps_ra@tax_table[,5] == "NA"] <- "Unassigned"
ps_ra@tax_table[,6][ps_ra@tax_table[,6] == "NA"] <- "Unassigned"

ps_ra2@tax_table[,1][ps_ra2@tax_table[,1] == "NA"] <- "Unassigned"
ps_ra2@tax_table[,3][ps_ra2@tax_table[,3] == "NA"] <- "Unassigned"
ps_ra2@tax_table[,4][ps_ra2@tax_table[,4] == "NA"] <- "Unassigned"
ps_ra2@tax_table[,5][ps_ra2@tax_table[,5] == "NA"] <- "Unassigned"
ps_ra2@tax_table[,6][ps_ra2@tax_table[,6] == "NA"] <- "Unassigned"


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

p.family = plot_bar(ps_ra, fill = "Family",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + theme_bw() +
  theme(axis.title = element_text(size = 24,face = "bold"), axis.text = element_text(size=18),
        legend.title = element_text(size=20), legend.text = element_text(size=16)) #+ scale_fill_manual(values = pal[c(1:9,11)])  
p = ggplot_build(p.family)
family.colors = unique(p$data[[1]]["fill"])[,1]
family.colors[26] <- "#3c3e44"

plot_bar(ps_ra, fill = "Family",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + theme_bw() +
  theme(axis.title = element_text(size = 24,face = "bold"), axis.text = element_text(size=18),
        legend.title = element_text(size=20), legend.text = element_text(size=16)) + scale_fill_manual(values = family.colors)  

ggsave("./Output/BarPlot_Fungal_Family_by_Island.png", height = 8, width = 12)

p.genus = plot_bar(ps_ra, fill = "Genus",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + theme_bw() +
  theme(axis.title = element_text(size = 24,face = "bold"), axis.text = element_text(size=18),
        legend.title = element_text(size=20), legend.text = element_text(size=16)) #+ scale_fill_manual(values = pal[c(1:9,11)])  
p = ggplot_build(p.genus)
genus.colors = unique(p$data[[1]]["fill"])[,1]
genus.colors[39] <- "#3c3e44"

plot_bar(ps_ra, fill = "Genus",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + #facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + theme_bw() +
  theme(axis.title = element_text(size = 24,face = "bold"), axis.text = element_text(size=18),
        legend.title = element_text(size=20), legend.text = element_text(size=16)) + scale_fill_manual(values = genus.colors)  
ggsave("./Output/BarPlot_Fungal_Genus_by_Island.png", height = 8, width = 12)



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

plot_bar(ps_ra2, x = "Island", fill = "Class") +
  geom_bar(stat = "identity") + coord_flip()  + facet_wrap(~Structure) +
  labs(x="Site",y="Relative Abundance") + lims(y=c(0,1)) + theme_bw() + scale_fill_manual(values = pal)
ggsave("./Output/BarPlot_Fungal_Class_by_Island_and_Structure.png", height = 8, width = 12)

plot_bar(ps_ra2, fill = "Order",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + lims(y=c(0,1)) + theme_bw() + scale_fill_manual(values = pal[c(1:10,12:17,11:12)], breaks = orders)
ggsave("./Output/BarPlot_Fungal_Order_by_Island_and_Structure.png", height = 8, width = 12)


# Family level for reviewer #1
plot_bar(ps_ra2, fill = "Family",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + lims(y=c(0,1)) + theme_bw()# + scale_fill_manual(values = pal[c(1:10,12:17,11:12)], breaks = orders)
ggsave("./Output/BarPlot_Fungal_Family_by_Island_and_Structure.png", height = 8, width = 12)

# Genus level for reviewer #1
plot_bar(ps_ra2, fill = "Genus",x="Island") +
  geom_bar(stat = "identity") + coord_flip()  + facet_wrap(~levels(sample_data(ps)$Structure)) +
  labs(x="Site",y="Relative Abundance") + lims(y=c(0,1)) + theme_bw()# + scale_fill_manual(values = pal[c(1:10,12:17,11:12)], breaks = orders)
ggsave("./Output/BarPlot_Fungal_Genus_by_Island_and_Structure.png", height = 8, width = 12)



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

rowSums(order.df)

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

# Again, for fungal Class
otus = as(otu_table(ps), "matrix")
otus.df = as.data.frame(otus)
classes = str_remove(tax_table(ps)[,3], "c__")
classes[is.na(classes)] <- "NA"
classes[classes == "NA"] <- "Unassigned"
names(otus.df) <- classes
df = as.data.frame(t(otus.df))
df$Class = row.names(df)
class.df = df %>% group_by(Class) %>% summarise_all(sum) %>% as.data.frame()
row.names(class.df) <- class.df$Class
class.df <- class.df %>% select(-Class)
class.df = class.df[which(row.names(class.df) != "Unassigned"),]
# reorder Cols by structure
meta <- as(sample_data(ps),"data.frame")
meta = meta %>% arrange(Structure)
class.df = class.df[,meta$IlluminaName]
class.df = decostand(class.df,"total", MARGIN = 2)

class.matrix = as.matrix(class.df)
hm.pal = c("#ffffff",RColorBrewer::brewer.pal(8,"OrRd"))

# make colors based on Structure
structures = meta$Structure
str.cols = as.character(plyr::mapvalues(structures, from = levels(structures), to=c("#492c24","#397c1c","#3791b2")))

# Heatmap of orders by structure (for individual samples)
png("./Output/Heatmap_of_Class_by_Structure.png", width = 16,height = 16,units = 'in',res = 300)
heatmap_left(t(class.matrix), col = hm.pal, Rowv = NA, Colv=NA, labRow = NA, RowSideColors = str.cols,
             margins = c(50,10), cexCol = 5)
dev.off()



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

# fix phyla
tax_table(ps)[,2] <- paste0("p__",phyla)

# Taxonomy table ... frequency of each taxa (number of samples found in)
write.csv(table(tax_table(ps)),"./Output/TaxTable_no_of_Samples.csv", row.names = FALSE, quote = FALSE)
     

# PermANOVA ####
psra = transform_sample_counts(ps, function(x) x/sum(x))
island = as.character(sample_data(psra)$Island)
structure = as.character(sample_data(psra)$Structure)


permanova = adonis(otu_table(psra) ~ island * structure)
sink("./Output/Adonis_Table.txt")
permanova
sink(NULL)


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



# Venn diagram ####


pspa = transform_sample_counts(ps, function(abund) 1*(abund>0))  
hf.pa = subset_samples(pspa, Structure == "Hold Fast")
lf.pa = subset_samples(pspa, Structure == "Leaf")
vs.pa = subset_samples(pspa, Structure == "Vesicle")

# subset
area.hf = length(which(taxa_sums(hf.pa) > 0))
area.lf = length(which(taxa_sums(lf.pa) > 0))
area.vs = length(which(taxa_sums(vs.pa) > 0))

# find areas
hf.pa = subset_taxa(hf.pa,taxa_sums(hf.pa) > 0)
lf.pa = subset_taxa(lf.pa,taxa_sums(lf.pa) > 0)
vs.pa = subset_taxa(vs.pa,taxa_sums(vs.pa) > 0)

# find shared sets 
nhf_lf = sum(taxa_names(hf.pa) %in% taxa_names(lf.pa))
nlf_vs = sum(taxa_names(lf.pa) %in% taxa_names(vs.pa))
nhf_vs = sum(taxa_names(hf.pa) %in% taxa_names(vs.pa))

hf.and.lf = which(taxa_names(hf.pa) %in% taxa_names(lf.pa))
hf.and.vs = which(taxa_names(hf.pa) %in% taxa_names(vs.pa))
nhf_lf_vs = sum(hf.and.lf  %in% hf.and.vs )
hf.and.lf = length(which(taxa_names(hf.pa) %in% taxa_names(lf.pa)))
hf.and.vs = length(which(taxa_names(hf.pa) %in% taxa_names(vs.pa)))


dev.off()
png("./Output/VennDiagram_Structure.png", res = 300, width = 2000,height = 1600)
draw.triple.venn(area.hf,area.lf,area.vs,nhf_lf,nlf_vs,nhf_vs,nhf_lf_vs,
                 fill=pal[c(10,6,5)],
                 category = c("Holdfast","Leaf","Vesicle"),
                 alpha = .75, cex = rep(2,7), cat.cex = rep(2,3), cat.dist = .1)
dev.off()

