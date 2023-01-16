library(ggplot2)
library(ggpubr)
library(tidyverse)
library(tibble)
library(readxl)
library(RColorBrewer)
library(vegan)
library(VennDiagram)
library(gplots)
library(purrr)
library(RVenn)
library(labdsv)
library(indicspecies)
library(lubridate)
library(wesanderson)
library(viridis)
library(tmap)
library(sf)
library(spData)
library(ggmap)
library(maps)
library(mapdata)
library(maptools)
library(ggthemes)
library(marmap)
library(datasets)
library(data.table)
library(rstatix)
library(datarium)


## Figure 1 Map ##
world <- world

ak <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group), fill = NA, colour = "black", size = 1) +
  geom_point(aes(x=-149.6, y=68.6), size = 3) +
  coord_quickmap(xlim = c(-174, -138), 
                 ylim = c(50, 72), 
                 expand = TRUE) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16))

pdf('map.pdf', height = 5, width = 6)
print(ak)
dev.off()

#Significance comparison of samples
#Only using OTUs with more than 10 reads in all sample sites
toolm.tax <- read_tsv("toolm.tax") %>% 
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string = taxonomy, pattern = ";$", replacement = "")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")
toolm.otu.single <- read_tsv("toolm.otu", col_types = cols(Group=col_character()))
metadata<-read.table("toolm.meta.txt", header = TRUE, sep = "\t")
sample <- toolm.otu.single$Group
toolm.otu.single <- toolm.otu.single[-which(names(toolm.otu.single) %in% c("label", "numOtus", "Group"))]
toolm.otu.single <- toolm.otu.single[,colSums(toolm.otu.single)>10]
toolm.otu.single.names <- cbind(toolm.otu.single, sample)

rel_abun_rowSums <- function(x, y) {
  mutate(y, n=rowSums(x))
}

total_reads_single <- toolm.otu.single.names %>% 
  group_by(sample) %>% 
  group_map(rel_abun_rowSums) %>%
  bind_rows()

toolm.otu.single.names <- toolm.otu.single.names %>% mutate(reads=total_reads_single$n)

otu.table.relabund <- toolm.otu.single.names %>%
  pivot_longer(cols = c(-sample, -reads), names_to = "otu", values_to = "count") %>%
  mutate(rel_abund=(count/reads*100))

otu.meta <- inner_join(toolm.otu.single.names, metadata, by=c("sample"="group"))

otu.meta.fe <- dplyr::filter(otu.meta, iron_abun == "High Fe")
otu.meta.fe <- arrange(otu.meta.fe, sample)
fe.otumetarownames <- otu.meta.fe$sample
toolm.otu.fe <- otu.meta.fe[-which(names(otu.meta.fe) %in% c(colnames(metadata), "sample"))]
row.names(toolm.otu.fe) <- fe.otumetarownames
meta.fe <- dplyr::filter(metadata, iron_abun == "High Fe")
meta.fe <- arrange(meta.fe, group)


otu.meta<-arrange(otu.meta, sample)
metadata<-arrange(metadata, group)
otumetarownames <- otu.meta$sample
toolm.otu <- otu.meta[-which(names(otu.meta) %in% c(colnames(metadata), "sample"))]
row.names(toolm.otu) <- otumetarownames

# @input reformatted OTU table and OTU taxonomy table
# @return % rel abund and taxonomy for every OTU from every sample
otu.tax <- inner_join(otu.table.relabund, toolm.tax)

# @input metadata for all samples and % rel abund and taxonomy for every OTU from every sample
# @return complete dataframe of all OTUs from every sample with everything we want to know about each OTU
otu.tax.meta <- inner_join(otu.tax, metadata, by=c("sample"="group"))

#toolm.otu.single <- read_tsv("toolm.otu", col_types = cols(Group=col_character())) %>%
  #select(-label, -numOtus)
#sample <- toolm.otu.single$Group
#toolm.otu.single <- select(toolm.otu.single, -Group)
#toolm.otu.single <- toolm.otu.single[,colSums(toolm.otu.single)>10]
#toolm.otu.single <- cbind(toolm.otu.single, sample)

#rel_abun_rowSums <- function(x, y) {
#  mutate(y, n=rowSums(x))
#}
#the toolm.otu.single below is the same as toolm.otu.single.names above
#total_reads_single <- toolm.otu.single %>% 
#  group_by(sample) %>% 
#  group_map(rel_abun_rowSums) %>%
#  bind_rows()

#toolm.otu.single <- toolm.otu.single %>% mutate(reads=total_reads_single$n)

#toolm.otu.single <- toolm.otu.single %>%
#  pivot_longer(cols = c(-sample, -reads), names_to = "otu", values_to = "count") %>%
#  mutate(rel_abund=(count/reads*100))

#ANOSIM Analysis
m_toolm.otu <- as.matrix(toolm.otu)
ano <- anosim(m_toolm.otu, metadata$iron_abun, distance = "bray", permutations = 9999)
ano
m_toolm.otu.fe <- as.matrix(toolm.otu.fe)
ano.feature <- anosim(m_toolm.otu.fe, meta.fe$feature_type, distance = "bray", permutations = 9999)
ano.feature

#Another way to do the ANOSIM, gives same result
#Dissimilarity matrix made here is for all sample heat map below
toolm.bc.dist <- vegdist(toolm.otu, distance="bray")
toolm.bc.env <- with(metadata, anosim(toolm.bc.dist, iron_abun))
summary(toolm.bc.env)
plot(toolm.bc.env)
anosim(toolm.bc.dist, metadata$iron_abun, permutations = 1000)
toolm.fe.bc.dist <- vegdist(toolm.otu.fe, distance="bray")
toolm.fe.bc.env <- with(meta.fe, anosim(toolm.fe.bc.dist, feature_type))
summary(toolm.fe.bc.env)
plot(toolm.fe.bc.env)

#Indicator species analysis
iron_abun <- metadata$iron_abun
ind <- multipatt(toolm.otu, iron_abun, func = "r.g", print.perm = TRUE, control = how(nperm = 1000))
ind <- as.data.frame(ind$sign, keep.rownames=TRUE)

#adjusted p-values for multiple comparisons
ind <- mutate(ind, p.value.bh=p.adjust(p.value, method="BH"))
ind <- filter(ind, p.value.bh <= 0.05)
indi.sp <- rownames_to_column(ind, "otu")
indi.sp <- inner_join(indi.sp, toolm.tax, by=c("otu"="otu"))
indi.sp <- indi.sp %>% rename_with(make.names)
indi.sp <- as.data.frame(indi.sp)
indi.sp$s.High.Fe = as.logical(indi.sp$s.High.Fe)
indi.sp$s.Low.Fe = as.logical(indi.sp$s.Low.Fe)

#Diversity statistics
H <- diversity(toolm.otu.single)
richness <- specnumber(toolm.otu.single)
eveness <- H/log(richness)
alpha.div <- cbind(sample, H, richness, eveness, metadata)

#Plot diversity
div_test_iron <- alpha.div %>% t_test(H ~ iron_abun) %>% add_significance() %>% add_xy_position(x = "iron_abun")

div_test_feature <- alpha.div %>% filter(iron_abun == "High Fe") %>% t_test(H ~ feature_type) %>% add_significance() %>% add_xy_position(x = "feature_type")

div_test_pond <- alpha.div %>% filter(feature_type == "Pond") %>% t_test(H ~ iron_abun) %>% add_significance() %>% add_xy_position(x = "iron_abun")

alpha.div.ponds <- ggplot(dplyr::filter(alpha.div, feature_type == "Pond")) +
  geom_boxplot(aes(x=iron_abun, y=H)) +
  geom_jitter(aes(x=iron_abun, y=H)) +
  stat_pvalue_manual(div_test_pond, tip.length = 0) + labs(subtitle = get_test_label(div_test_pond, detailed = FALSE)) +
  labs(#title = "High and Low Fe Ponds Only" 
    x = NULL, y = NULL) +
  ylim(2.5, 6.5) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.key=element_blank(), legend.position = c(0.8,0.2),
        axis.text.y = element_text(color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 16))

alpha.div.highfe <- ggplot(dplyr::filter(alpha.div, iron_abun == "High Fe")) +
  geom_boxplot(aes(x=feature_type, y=H)) +
  geom_jitter(aes(x=feature_type, y=H)) +
  stat_pvalue_manual(div_test_feature, tip.length = 0) +
  labs(#title = "High Fe Sites by Feature Type Only" 
    x = NULL, y = NULL) +
  ylim(2.5, 6.5) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.key=element_blank(), legend.position = c(0.8,0.2),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"))

alpha.div.all <- ggplot(alpha.div) +
  geom_boxplot(aes(x=iron_abun, y=H)) +
  geom_jitter(aes(x=iron_abun, y=H)) +
  stat_pvalue_manual(div_test_iron, tip.length = 0) + labs(subtitle = get_test_label(div_test_iron, detailed = FALSE)) +
  labs(#title = "All High and Low Fe Sites"
       x = NULL, y = "Shannon's 'H'") +
  ylim(2.5, 6.5) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.key=element_blank(), legend.position = c(0.8,0.2),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"))

pdf(file="diversity.pdf", height = 4, width = 12)
print(ggarrange(alpha.div.all, alpha.div.highfe, alpha.div.ponds,  
                labels = c("A", "B", "C"), 
                ncol = 3, nrow = 1))
dev.off()

##Fe site indicator species
fe.indi.sp <- dplyr::filter(indi.sp, 
              s.High.Fe == "TRUE" & stat > 0.5 & p.value.bh < 0.01 & grepl("*ceae", family))
nofe.indi.sp <- dplyr::filter(indi.sp, 
              s.Low.Fe == "TRUE" & stat > 0.5 & p.value.bh < 0.01 & grepl("*ceae", family))

# plot indicator species from high and low Fe sites
x <- dplyr::filter(indi.sp, 
             s.Low.Fe == "TRUE" & stat > 0.5 & p.value.bh < 0.01 & grepl("*ceae", family))
n_distinct(x$family)
nrow(x)
x_otu <- x$otu

x_relabund <- otu.table.relabund[otu.table.relabund$otu %in% x_otu, ]
x_relabund.meta <- inner_join(x_relabund, metadata, by=c("sample"="group")) %>% 
  filter(., iron_abun == "Low Fe", rel_abund > 0.1) %>% 
  inner_join(., x, by=c("otu"="otu"))
x_unique <- x_relabund.meta %>% distinct(otu, .keep_all = TRUE)
x_unique.family <- x_relabund.meta %>% distinct(family, .keep_all = TRUE)

indi.sp.plot.nofe <- ggplot(x_relabund.meta) +
  geom_point(aes(x=stat, y=reorder(family, stat), size=rel_abund)) +
  labs(x="Stat Value (OTU is more strongly associated with Low Fe site -->)", y=NULL, size="Relative Abundance") +
  theme(text = element_text(size=17),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.key=element_blank(), legend.position = c(0.8,0.3),
        axis.text.y = element_text(face="italic", color = "black"),
        axis.text.x = element_text(color = "black"))

y <- dplyr::filter(indi.sp, 
              s.High.Fe == "TRUE" & stat > 0.5 & p.value.bh < 0.01 & grepl("*ceae", family))
n_distinct(y$family)
nrow(y)
y_otu <- y$otu

y_relabund <- otu.table.relabund[otu.table.relabund$otu %in% y_otu, ]
y_relabund.meta <- inner_join(y_relabund, metadata, by=c("sample"="group")) %>% 
  filter(., iron_abun == "High Fe", rel_abund > 0.1) %>% 
  inner_join(., y, by=c("otu"="otu"))
y_unique <- y_relabund.meta %>% distinct(otu, .keep_all = TRUE)
y_unique.family <- y_relabund.meta %>% distinct(family, .keep_all = TRUE)

indi.sp.plot.fe <- ggplot(y_relabund.meta) +
  geom_point(aes(x=stat, y=reorder(family, stat), size=rel_abund)) +
  labs(x="Stat Value (OTU is more strongly associated with High Fe site -->)", y="Family", size="Relative Abundance") +
  theme(text = element_text(size=17),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(), legend.position = c(0.8,0.3),
        axis.text.y = element_text(face="italic", color = "black"),
        axis.text.x = element_text(color = "black"))

pdf(file="IndicatorSpecies.pdf", height = 8, width = 20)
print(ggarrange(indi.sp.plot.fe, indi.sp.plot.nofe,  
                labels = c("A", "B"),
                ncol = 2, nrow = 1))
dev.off()

#Plotting of bray-curtis dissimilarities using a heatmap
## #f35a00 is orange and #0d01ff is blue used in Fig. 1 map (color blind friendly)
hmp.color.rows <- c("#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", 
                    "#0d01ff", "#0d01ff", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#f35a00",
                    "#0d01ff", "#f35a00", "#0d01ff", "#f35a00", "#f35a00", "#f35a00", "#0d01ff", "#0d01ff",
                    "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#0d01ff", "#0d01ff", "#0d01ff", "#f35a00", "#f35a00",
                    "#0d01ff", "#0d01ff", "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#0d01ff", 
                    "#f35a00", "#f35a00", "#f35a00", "#f35a00", "#0d01ff", "#0d01ff", "#0d01ff", "#f35a00",
                    "#0d01ff", "#0d01ff", "#0d01ff", "#0d01ff", "#f35a00", "#f35a00", "#f35a00", "#0d01ff", "#0d01ff", "#0d01ff",
                    "#f35a00", "#f35a00", "#f35a00", "#0d01ff", "#f35a00", "#f35a00", "#f35a00")
## #F5C75D is gold, #34DDEB is cyan, #5DDB00 is green (Color blind friendly)   
hmp.color.cols <- c("#F5C75D", "#34DDEB", "#34DDEB", "#34DDEB", "#34DDEB", "#5DDB00", "#34DDEB", "#34DDEB", 
                    "#F5C75D", "#34DDEB", "#34DDEB", "#34DDEB", "#34DDEB", "#5DDB00", "#F5C75D", "#34DDEB", "#34DDEB",
                    "#5DDB00", "#5DDB00", "#34DDEB", "#34DDEB", "#34DDEB", "#34DDEB",
                    "#F5C75D", "#F5C75D", "#F5C75D", "#34DDEB", "#34DDEB", "#5DDB00", "#5DDB00", "#5DDB00", "#F5C75D", 
                    "#34DDEB", "#34DDEB", "#34DDEB", "#F5C75D", "#34DDEB", "#34DDEB", "#34DDEB", "#5DDB00", "#F5C75D", "#F5C75D",
                    "#34DDEB", "#34DDEB", "#5DDB00", "#5DDB00", "#5DDB00", "#5DDB00", "#34DDEB", "#34DDEB", "#34DDEB", "#34DDEB",
                    "#34DDEB", "#34DDEB", "#34DDEB", "#34DDEB", "#5DDB00", "#F5C75D", "#F5C75D", 
                    "#34DDEB", "#34DDEB", "#34DDEB", "#5DDB00", "#5DDB00", "#5DDB00", "#34DDEB",
                    "#5DDB00", "#5DDB00", "#34DDEB")
hmp.col <- colorRampPalette(brewer.pal(11, "RdYlBu"))

pdf("toolm_heatmap.pdf", width = 12, height = 12)
heatmap.2(as.matrix(toolm.bc.dist), col = hmp.col, 
          labRow = NULL, labCol = NULL, vline = NULL, hline = NULL, tracecol = NULL,
          RowSideColors = hmp.color.rows, ColSideColors = hmp.color.cols, cexRow = 0.01, cexCol = 0.01,
          key = TRUE,
          keysize = 1.5,
          density.info=c("histogram","density","none"),
          densadj = 0.25,
          key.title = NULL,
          key.xlab = "Dissimilarity value",
          key.ylab = "Count")
dev.off()

#Heat map of High Fe sites only
#Uses Bray-Curtis dissimilarity matrix calculated above

##ADD COLOR BARS TO THE TOP FOR POND WET SEDGE AND SEEP

pdf("toolm_fe_heatmap.pdf", width = 12, height = 12)
heatmap.2(as.matrix(toolm.fe.bc.dist), col = hmp.col, 
          labRow = NULL, labCol = NULL, vline = NULL, hline = NULL, tracecol = NULL,
          cexRow = 0.01, cexCol = 0.01,
          key = TRUE,
          keysize = 1.5,
          density.info=c("histogram","density","none"),
          densadj = 0.25,
          key.title = NULL,
          key.xlab = "Dissimilarity value",
          key.ylab = "Count")
dev.off()

#Site water quality parameters stats
summary(filter(metadata, iron_abun == "Low Fe"))
summary(filter(metadata, iron_abun == "Low Fe"))
summary(filter(metadata, feature_type == "Seep"))
summary(filter(metadata, feature_type == "Pond"))
summary(filter(metadata, feature_type == "Wet Sedge"))
summary(metadata)
ggboxplot(metadata, x = "date", y = "oxy_uM", add = "jitter")

# NMDS of High Fe sites
toolm.fe.bc.metanmds <- metaMDS(m_toolm.otu.fe, distance = "bray")
nmds.fe.data.scores <- as.data.frame(scores(toolm.fe.bc.metanmds)$sites)
nmds.fe.meta <- cbind(nmds.fe.data.scores, meta.fe)
toolm.fe.bc.metanmds$stress

toolm.fe.nmds <- ggplot(nmds.fe.meta, aes(x=NMDS1, y=NMDS2, color=feature_type)) + 
  geom_point(size=4) + 
  scale_color_manual(name=NULL, 
                     values=brewer.pal(n=5, name="Dark2"), 
                     breaks=c("Pond", "Seep", "Wet Sedge"), 
                     labels=c("Pond", "Seep", "Wet Sedge")) +
  labs(x="PCo Axis 1", y="PCo Axis 2") +
  theme_classic()

pdf('toolm_Fe_PCoA_BC.pdf', height = 5, width = 6)
print(toolm.fe.pcoa)
dev.off()

# NMDS of All sites
#Pair this plot with the ANOSIM data that shows significant difference between high and low fe sites
#Uses m_toolm.otu matrix made above
toolm.bc.metanmds <- metaMDS(m_toolm.otu, distance = "bray")
nmds.data.scores <- as.data.frame(scores(toolm.bc.metanmds)$sites)
nmds.meta <- cbind(nmds.data.scores, metadata)
toolm.bc.metanmds$stress

toolm.nmds <- ggplot(nmds.meta, aes(x=NMDS1, y=NMDS2, shape=iron_abun, color=feature_type)) + 
  geom_point(size=4) + 
  scale_color_manual(name=NULL, 
                     values=brewer.pal(n=5, name="Dark2"), 
                     breaks=c("Pond", "Seep", "Wet Sedge"), 
                     labels=c("Pond", "Seep", "Wet Sedge")) +
  scale_shape_manual(name=NULL,
                     values=c(19, 17, 15),
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(x="NMDS Axis 1", y="NMDS Axis 2") +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(), legend.position = c(0.9,0.2),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"))

pdf('toolm_nmds_BC.pdf', height = 5, width = 6)
print(toolm.fe.pcoa)
dev.off()

vector.meta <- metadata[-which(names(metadata) %in% c("date","dna_ex", "lat", "long"))]

en = envfit(toolm.bc.metanmds, vector.meta, permutations = 999, na.rm = TRUE)
plot(en)

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en, fill = 0.35)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en, fill = 0.35)

en_coord_cat_rownames = row.names(en_coord_cat)
en_coord_cat_labels = cbind(en_coord_cat, en_coord_cat_rownames)
en_coord_cont_rownames <- c("Day of Year", "Methane (nM)", "Overlying Water Fe2+ (uM)", "Temperature (C)", 
                                      "Dissolved O2 (mg L)", "Dissolved O2 (uM)", "pH", "Total Extract Fe(III)",
                                      "0.5 M HCl Extract Fe(III)")
en_coord_cont = cbind(en_coord_cont, en_coord_cont_rownames)
en_coord_cont <- filter(en_coord_cont, en_coord_cont_rownames %in% c("Methane (nM)","Dissolved O2 (mg L)", "pH", 
                                                                     "Total Extract Fe(III)",
                                                                     "0.5 M HCl Extract Fe(III)"))
summary(nmds.meta)

env.nmds <- ggplot(nmds.meta) +
  geom_point(aes(x=NMDS1, y=NMDS2, color=feature_type, shape=iron_abun), size = 3) +
  scale_color_manual(name=NULL, 
                     values=brewer.pal(n=5, name="Dark2"), 
                     breaks=c("Pond", "Seep", "Wet Sedge"), 
                     labels=c("Pond", "Seep", "Wet Sedge")) +
  scale_shape_manual(name=NULL,
                     values=c(19, 17, 15),
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  #geom_point(data = dplyr::filter(en_coord_cat_labels, grepl("feature_type*", en_coord_cat_rownames)), 
                                #  aes(x=NMDS1, y=NMDS2), color = "pink") +
  #geom_text(data = dplyr::filter(en_coord_cat_labels, grepl("feature_type*", en_coord_cat_rownames)), 
                                 #aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                                 #fontface = "bold", label = row.names(en_coord_cat_labels)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = en_coord_cont$en_coord_cont_rownames,
            vjust = 1, hjust = -0.05) +
  labs(x="NMDS Axis 1", y="NMDS Axis 2") +
  #xlim(-1.25,2) +
  theme(text = element_text(size=24),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(), legend.position = c(0.9,0.2),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"))

pdf('toolm_env_nmds.pdf', height = 10, width = 12)
print(env.nmds)
dev.off()

x=dplyr::filter(en_coord_cat_labels, grepl("feature_type*", en_coord_cat_rownames))


## Finding the "Core Community" of Toolm samples
toolm.tax <- read_tsv("toolm.tax") %>% 
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string = taxonomy, pattern = ";$", replacement = "")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")
taxrownames <- toolm.tax$otu
toolm.tax <- toolm.tax[-which(names(toolm.tax) %in% c("size", "otu"))]
row.names(toolm.tax) <- taxrownames

# Fe sites 
toolm.otu <- read_tsv("toolm.otu", col_types = cols(Group=col_character()))
sample <- toolm.otu$Group
toolm.otu <- toolm.otu[,-which(names(toolm.otu) %in% c("label", "Group", "numOtus"))]
toolm.otu <- toolm.otu/rowSums(toolm.otu)
toolm.otu <- toolm.otu*100
metadata<-read.table("toolm.meta.txt", header = TRUE, sep = "\t")
metadata.fe <- metadata[,-c(2,3,4,6:11)]
metadata.fe.names <- metadata.fe$group
metadata.fe.names <- as.vector(metadata.fe.names)
View(metadata.fe.names)
toolm.otu <- cbind(toolm.otu, metadata.fe.names)
toolm.otu.meta <- inner_join(toolm.otu, metadata.fe, by=c("metadata.fe.names" = "group"))
row.names(toolm.otu) <- sample
toolm.otu[1:5,1:5]
has_rownames(toolm.otu.fe)
toolm.otu.fe <- filter(toolm.otu.meta, iron_abun == "High Fe")
toolm.otu.fe[1:5,1:5]
otu.fe.names <- toolm.otu.fe$metadata.fe.names

toolm.otu.fe <- toolm.otu.fe[,-c(50087:50088)]

row.names(toolm.otu.fe) <- otu.fe.names
toolm.otu.fe[1:5,1:5]

dim(toolm.otu.fe)[1] #always reports the dimensions
fe.core <- colSums(toolm.otu.fe > 0) == dim(toolm.otu.fe)[1]
fe.core <- as.data.frame(fe.core)
fe.core.names <- row.names(fe.core)
I <- which(fe.core[,1])
fe.core.names[I]
toolm.otu.fe.means <- colMeans(toolm.otu.fe[,I])
toolm.otu.fe.means <- as.data.frame(toolm.otu.fe.means)
has_rownames(toolm.otu.fe.means)
toolm.fe.core <- merge(toolm.tax, toolm.otu.fe.means, by = 0)
write.csv(toolm.fe.core, "toolm.fe.core.csv")

# No-Fe sites core community
toolm.otu <- read_tsv("toolm.otu", col_types = cols(Group=col_character()))
sample <- toolm.otu$Group
toolm.otu <- toolm.otu[,-c(1:3)]
toolm.otu <- toolm.otu/rowSums(toolm.otu)
toolm.otu <- toolm.otu*100
metadata.nofe<-read.table("toolm.meta.txt", header = TRUE, sep = "\t")


metadata.nofe <- metadata[,-c(2,3,4,6:11)]
metadata.nofe.names <- metadata.nofe$group
metadata.nofe.names <- as.vector(metadata.nofe.names)
View(metadata.nofe.names)
toolm.otu <- cbind(toolm.otu, metadata.nofe.names)
toolm.otu.meta <- inner_join(toolm.otu, metadata.nofe, by=c("metadata.nofe.names" = "group"))
row.names(toolm.otu) <- sample
toolm.otu[1:5,1:5]
has_rownames(toolm.otu.nofe)
toolm.otu.nofe <- filter(toolm.otu.meta, iron_abun == "Low Fe")
toolm.otu.nofe[1:5,1:5]
otu.nofe.names <- toolm.otu.nofe$metadata.nofe.names

toolm.otu.nofe <- toolm.otu.nofe[,-c(50087:50088)]

row.names(toolm.otu.nofe) <- otu.nofe.names
toolm.otu.nofe[1:5,1:5]

dim(toolm.otu.nofe)[1] #always reports the dimensions
nofe.core <- colSums(toolm.otu.nofe > 0) == dim(toolm.otu.nofe)[1]
nofe.core <- as.data.frame(nofe.core)
nofe.core.names <- row.names(nofe.core)
I <- which(nofe.core[,1])
nofe.core.names[I]
toolm.otu.nofe.means <- colMeans(toolm.otu.nofe[,I])
toolm.otu.nofe.means <- as.data.frame(toolm.otu.nofe.means)
has_rownames(toolm.otu.nofe.means)
toolm.nofe.core <- merge(toolm.tax, toolm.otu.nofe.means, by = 0)
write.csv(toolm.nofe.core, "toolm.nofe.core.csv")

# All toolm sites core community
toolm.otu <- read_tsv("toolm.otu", col_types = cols(Group=col_character()))
sample <- toolm.otu$Group
toolm.otu <- toolm.otu[,-c(1:3)]
toolm.otu <- toolm.otu/rowSums(toolm.otu)
toolm.otu <- toolm.otu*100
row.names(toolm.otu) <- sample
toolm.otu[1:5,1:5]

dim(toolm.otu)[1] #always reports the dimensions
all.core <- colSums(toolm.otu > 0) == dim(toolm.otu)[1]
all.core <- as.data.frame(all.core)
all.core.names <- row.names(all.core)
I <- which(all.core[,1])
all.core.names[I]
toolm.otu.means <- colMeans(toolm.otu[,I])
toolm.otu.means <- as.data.frame(toolm.otu.means)
has_rownames(toolm.otu.means)
toolm.core <- merge(toolm.tax, toolm.otu.means, by = 0)
write.csv(toolm.core, "toolm.core.csv")

filter(I) #filter my OTU table by I, especially after subsetting the OTU table by high Fe sites

##This worked from remote sensing Fe project R code
#High Fe sites core community
tfs.otu <- read_tsv("toolm.otu", col_types = cols(Group=col_character()))
sample <- tfs.otu$Group
tfs.otu.fe <- tfs.otu[,-which(names(tfs.otu) %in% c("label", "Group", "numOtus"))]
tfs.otu.fe <- tfs.otu.fe/rowSums(tfs.otu.fe)
tfs.otu.fe <- tfs.otu.fe*100
metadata<-read.table("toolm.meta.txt", header = TRUE, sep = "\t")

metadata.fe.names <- metadata$group
metadata.fe.names <- as.vector(metadata.fe.names)
tfs.otu.fe <- cbind(tfs.otu.fe, metadata.fe.names)
tfs.otu.meta <- inner_join(tfs.otu.fe, metadata, by=c("metadata.fe.names" = "group"))
row.names(tfs.otu.fe) <- sample
tfs.otu[1:5,1:5]
tfs.otu.fe <- filter(tfs.otu.meta, iron_abun == "High Fe")
otu.fe.names <- tfs.otu.fe$metadata.fe.names
tfs.otu.fe <- tfs.otu.fe[,-which(names(tfs.otu.fe) %in% c("metadata.fe.names", 
                                                                "wshd", "material", "site", "feature_type",
                                                                "iron_abun", "doy", "date", "methane_nM",
                                                                "water_fe2_uM", "temp", "oxy_mgL", "oxy_uM", "pH", 
                                                                "totfe3_ext", "X05fe3_ext", "dna_ex"))]

row.names(tfs.otu.fe) <- otu.fe.names
tfs.otu.fe[1:5,1:5]
dim(tfs.otu.fe)[1] #always reports the dimensions
tfs.fe.core <- colSums(tfs.otu.fe > 0) == dim(tfs.otu.fe)[1]
tfs.fe.core <- as.data.frame(tfs.fe.core)
tfs.fe.core.names <- row.names(tfs.fe.core)
I <- which(tfs.fe.core[,1])
tfs.fe.core.names[I]
tfs.otu.fe.means <- colMeans(tfs.otu.fe[,I])
tfs.otu.fe.means <- as.data.frame(tfs.otu.fe.means)
tfs.fe.names <- rownames(tfs.otu.fe.means)
rownames(tfs.otu.fe.means) <- c()
tfs.otu.fe.means <- cbind(tfs.otu.fe.means, tfs.fe.names)
tfs.fe.core <- inner_join(toolm.tax, tfs.otu.fe.means, by = c("otu" = "tfs.fe.names"))
tfs.fe.core$genus <- factor(tfs.fe.core$genus, levels = tfs.fe.core$genus[order(tfs.fe.core$tfs.otu.fe.means)])

#Plotting the Fe core community
fe.core <- ggplot(tfs.fe.core, aes(x=tfs.otu.fe.means, y=genus)) + 
  geom_point(size=4) +
  scale_x_log10() +
  ylab("Taxa") +
  xlab("Mean Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#Low Fe core community
tfs.otu <- read_tsv("toolm.otu", col_types = cols(Group=col_character()))
sample <- tfs.otu$Group
tfs.otu.nofe <- tfs.otu[,-which(names(tfs.otu) %in% c("label", "Group", "numOtus"))]
tfs.otu.nofe <- tfs.otu.nofe/rowSums(tfs.otu.nofe)
tfs.otu.nofe <- tfs.otu.nofe*100
metadata<-read.table("toolm.meta.txt", header = TRUE, sep = "\t")

metadata.nofe.names <- metadata$group
metadata.nofe.names <- as.vector(metadata.nofe.names)
tfs.otu.nofe <- cbind(tfs.otu.nofe, metadata.nofe.names)
tfs.otu.meta.nofe <- inner_join(tfs.otu.nofe, metadata, by=c("metadata.nofe.names" = "group"))
row.names(tfs.otu.nofe) <- sample
tfs.otu.nofe <- filter(tfs.otu.meta.nofe, iron_abun == "Low Fe")
otu.nofe.names <- tfs.otu.nofe$metadata.nofe.names
tfs.otu.nofe <- tfs.otu.nofe[,-which(names(tfs.otu.nofe) %in% c("metadata.nofe.names", 
                                                                "wshd", "material", "site", "feature_type",
                                                                "iron_abun", "doy", "date", "methane_nM",
                                                                "water_fe2_uM", "temp", "oxy_mgL", "oxy_uM", "pH", 
                                                                "totfe3_ext", "X05fe3_ext", "dna_ex"))]

row.names(tfs.otu.nofe) <- otu.nofe.names

dim(tfs.otu.nofe)[1] #always reports the dimensions
tfs.nofe.core <- colSums(tfs.otu.nofe > 0) == dim(tfs.otu.nofe)[1]
tfs.nofe.core <- as.data.frame(tfs.nofe.core)
tfs.nofe.core.names <- row.names(tfs.nofe.core)
I <- which(tfs.nofe.core[,1])
tfs.nofe.core.names[I]
tfs.otu.nofe.means <- colMeans(tfs.otu.nofe[,I])
tfs.otu.nofe.means <- as.data.frame(tfs.otu.nofe.means)
tfs.nofe.names <- rownames(tfs.otu.nofe.means)
rownames(tfs.otu.nofe.means) <- c()
tfs.otu.nofe.means <- cbind(tfs.otu.nofe.means, tfs.nofe.names)
tfs.nofe.core <- inner_join(toolm.tax, tfs.otu.nofe.means, by = c("otu" = "tfs.nofe.names"))


#Plotting the nofe core community
nofe.core <- ggplot(tfs.nofe.core, aes(x=tfs.otu.nofe.means, y=genus)) + 
  geom_point(size=4) +
  scale_x_log10() +
  ylab("Taxa") +
  xlab("Mean Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

pdf("toolm_core.pdf", width = 6, height = 12)
print(ggarrange(fe.core, nofe.core,  
                labels = c("A", "B"), 
                ncol = 1, nrow = 2))
dev.off()

## Supplemental Figure 1 -- Family-level abundances
agg_fam<- otu.tax.meta %>%
  group_by(family, sample) %>%
  summarize(rel_abund = sum(rel_abund))%>%
  filter(mean(rel_abund) > 0.75) %>%
  inner_join(., metadata, by=c("sample" = "group"))

rhodo <- otu.tax.meta %>%
  filter(genus == "Rhodoferax")

wes_colors <- wes_palette("FantasticFox1", n = 28, type = c("continuous"))
agg_fam_box <- ggplot(agg_fam) +
  geom_boxplot(aes(x=rel_abund, y=family, color=family), outlier.shape=NA, lwd=1) +
  geom_jitter(aes(x=rel_abund, y=reorder(family, rel_abund)), shape=19, size=1) +
  scale_color_manual(values = wes_colors) +
  scale_x_log10() +
  labs(x="Relative Abundance (%)", y="Family") +
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  legend.position = "none", legend.key=element_blank(), 
  axis.text.y = element_text(face="italic", color = "black"),
  axis.text.x = element_text(color = "black"))

pdf(file="supplemental_fig1.pdf", height = 8, width = 6)
print(agg_fam_box)
dev.off()

# High Fe sites only
agg_fam_fe<- otu.tax.meta %>%
  filter(., iron_abun == "High Fe") %>%
  group_by(family, sample) %>%
  summarize(rel_abund = sum(rel_abund))%>%
  filter(mean(rel_abund) > 0.75) %>%
  inner_join(., metadata, by=c("sample" = "group"))

rhodo <- otu.tax.meta %>%
  filter(genus == "Rhodoferax")

wes_colors <- wes_palette("FantasticFox1", n = 28, type = c("continuous"))
agg_fam_box_fe <- ggplot(agg_fam_fe) +
  geom_boxplot(aes(x=rel_abund, y=family, color=family), outlier.shape=NA, lwd=1) +
  geom_jitter(aes(x=rel_abund, y=reorder(family, rel_abund)), shape=19, size=1) +
  scale_color_manual(values = wes_colors) +
  scale_x_log10() +
  labs(x="Relative Abundance (%)", y="Family") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none", legend.key=element_blank(), 
        axis.text.y = element_text(face="italic", color = "black"),
        axis.text.x = element_text(color = "black"))

pdf(file="supplemental_fig1_feonly.pdf", height = 8, width = 6)
print(agg_fam_box_fe)
dev.off()

# Significance testing of families
agg_phylum_meta <- otu.tax.meta %>% filter(iron_abun != "Low Fe") %>%
  group_by(sample, phylum) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  filter(rel_abund > 0.5) %>%
  inner_join(., metadata, by=c("sample" = "group"))


toolm_agg_phylum <- ggplot(agg_phylum_meta) +
  geom_boxplot(aes(x=reorder(phylum, -rel_abund), y=rel_abund)) +
  geom_jitter(aes(x=reorder(phylum, -rel_abund), y=rel_abund, color=feature_type), shape=19, size=2, 
              position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Pond", "Seep", "Wet Sedge"),
                     labels=c("Pond", "Seep", "Wet Sedge")) +
  scale_y_log10() +
  labs(x="Phylum", y="Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=12, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.key=element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(color="black", size=15))

pdf(file="toolm_agg_phylum.pdf", height = 4, width = 7)
print(toolm_agg_phylum)
dev.off()

#Using abundance
agg_phylum_meta_abun <- otu.abun.tax %>% filter(iron_abun != "Low Fe") %>%
  group_by(sample, phylum) %>%
  summarize(abun = sum(abun)) %>%
  filter(abun > 5) %>%
  inner_join(., metadata, by=c("sample" = "group"))


toolm_agg_phylum_abun <- ggplot(agg_phylum_meta_abun) +
  geom_boxplot(aes(x=reorder(phylum, -abun), y=abun)) +
  geom_jitter(aes(x=reorder(phylum, abun), y=abun, color=feature_type), shape=19, size=2, 
              position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Pond", "Seep", "Wet Sedge"),
                     labels=c("Pond", "Seep", "Wet Sedge")) +
  scale_y_log10() +
  labs(x="Phylum", y="Abundance (Reads ng DNA-1 gww-1)") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=12, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.key=element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(color="black", size=15))

pdf(file="toolm_agg_phylum.pdf", height = 4, width = 7)
print(toolm_agg_phylum)
dev.off()

##Comparison of phyla and classes that are abundant in Arctic soils
proteo<-filter(otu.tax.meta, phylum%in%"Proteobacteria") %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))
summary(proteo)

acido<-filter(otu.tax.meta, phylum%in%"Acidobacteriota") %>%
  group_by(sample) %>%
  summarise(acido_rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group")) %>%
  filter(iron_abun == "High Fe") %>%
  rename(acido_sample = sample)

actino<-filter(otu.tax.meta, phylum%in%"Actinobacteriota") %>%
  group_by(sample) %>%
  summarise(actino_rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group")) %>%
  filter(iron_abun == "High Fe") %>%
  rename(actino_sample = sample)

actino_acido_mean_Fe <- data_frame(actino$actino_sample, actino$actino_rel_abund, acido$acido_sample, acido$acido_rel_abund) %>%
  mutate(sum=(actino$actino_rel_abund + acido$acido_rel_abund)) %>%
  summarize(mean(sum))

actino_acido_sum <- data_frame(actino$actino_sample, actino$actino_rel_abund, acido$acido_sample, acido$acido_rel_abund) %>%
  mutate(sum=(actino$actino_rel_abund + acido$acido_rel_abund)) 

cyano<-filter(otu.tax.meta, phylum%in%c("Cyanobacteria")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

alpha<-filter(otu.tax.meta, class%in%"Alphaproteobacteria") %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

verruco<-filter(otu.tax.meta, phylum%in%c("Verrucomicrobiota")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

comma<-filter(otu.tax.meta, family%in%c("Comamonadaceae")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

proteo_feabund<-ggplot(proteo) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=iron_abun, y=rel_abund, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  labs(title="Proteobacteria", x="Iron Abundance",
       y="Relative Abundance (%)") +
  theme_classic()

acido_feabund<-ggplot(acido) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=iron_abun, y=rel_abund, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  labs(title="Acidobacteriota", x="Iron Abundance",
       y="Relative Abundance (%)") +
  theme_classic()

cyano_feabund<-ggplot(cyano) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=iron_abun, y=rel_abund, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  labs(title="Cyanobacteria", x="Iron Abundance",
       y="Relative Abundance (%)") +
  theme_classic()

actino_feabund<-ggplot(actino) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=iron_abun, y=rel_abund, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  labs(title="Actinobacteriota", x="Iron Abundance",
       y="Relative Abundance (%)") +
  theme_classic()

verruco_feabund<-ggplot(verruco) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=iron_abun, y=rel_abund, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  labs(title="Verrucomicrobiota", x="Iron Abundance",
       y="Relative Abundance (%)") +
  theme_classic()

alpha_feabund<-ggplot(alpha) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=iron_abun, y=rel_abund, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  labs(title="Alphaproteobacteria", x="Iron Abundance",
       y="Relative Abundance (%)") +
  theme_classic()

#Taxon of interest in each feature
otu.tax.meta.fe <- filter(otu.tax.meta, iron_abun != "Low Fe")
otu.tax.meta.nofe <- filter(otu.tax.meta, iron_abun != "High Fe")

feox<-filter(otu.tax.meta, genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella")) %>%
  filter(., rel_abund>0) %>%
  group_by(sample) %>%
  summarize(sum=sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

feox.fe <- filter(otu.tax.meta.fe, genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella")) %>%
  filter(., rel_abund>0) %>%
  group_by(sample, genus) %>%
  summarize(sum=sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

feox_otu <- filter(otu.tax.meta, genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella") &
                     rel_abund>0) %>%
  group_by(sample, otu) %>%
  summarize(sum = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group")) %>%
  filter(sum > 0.05)

feox_OKS_otu_time <- ggplot(feox) +
  geom_boxplot(aes(x=factor(date), y=rel_abund), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=factor(date), y=rel_abund, color=genus)) +
  scale_y_log10() + 
  labs(x="Sample Date", y="Relative Abundance (%)", color="OTU") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.key=element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(color="black", size=15))

feox_OKS<-filter(otu.tax.meta.fe, genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella") & wshd == "OKS") %>%
  filter(., rel_abund>0) %>%
  group_by(sample, genus) %>%
  summarize(sum=sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

gal_OKS<-filter(otu.tax.meta.fe, family%in%"Gallionellaceae" & wshd == "OKS") %>%
  filter(., rel_abund>0) %>%
  group_by(sample, family)

lepto_OKS<-filter(otu.tax.meta.fe, genus%in%"Leptothrix"& wshd == "OKS") %>%
  filter(., rel_abund>0) %>%
  group_by(sample, genus)

feox_OKS<- rbind(gal_OKS, lepto_OKS)

feox_OKS_time <- ggplot(feox_OKS) +
  #geom_boxplot(aes(x=factor(rdate), y=sum), outlier.shape = NA) +
  geom_point(aes(x=rdate, y=rel_abund, color=family), se=T) +
  geom_smooth(aes(x=rdate, y=rel_abund, color=family), se=T) +
  #geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 4, jitter.width = 0.4), 
              #aes(x=rdate, y=rel_abund, color=family)) +
  scale_y_log10() + 
  #ylim(0.01, 10) +
  labs(x="Sample Date", y="Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.key=element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(color="black", size=15))

feox_TFS<-filter(otu.tax.meta.fe, wshd == "TFS" & 
                   genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella") &
                   date != 20190622 &
                   date != 20190713) %>%
  group_by(sample, genus) %>%
  summarize(sum=sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

gal_TFS<-filter(otu.tax.meta.fe, family%in%"Gallionellaceae" & wshd == "TFS") %>%
  filter(., rel_abund>0) %>%
  group_by(sample, family)

lepto_TFS<-filter(otu.tax.meta.fe, genus%in%"Leptothrix"& wshd == "TFS") %>%
  filter(., rel_abund>0) %>%
  group_by(sample, genus)

feox_TFS<- rbind(gal_TFS, lepto_TFS)

feox_TFS_time <- ggplot(feox_TFS) +
  #geom_boxplot(aes(x=factor(rdate), y=sum), outlier.shape = NA) +
  geom_smooth(aes(x=rdate, y=rel_abund, color=family), se=T) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 4, jitter.width = 0.4), 
              aes(x=rdate, y=rel_abund, color=family)) +
  scale_y_log10() + 
  labs(x="Sample Date", y="Relative Abundance (%)", color = "Genera") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.key=element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(color="black", size=15))

pdf(file="Figure6.pdf", height = 5, width = 9)
print(ggarrange(feox_OKS_time, feox_TFS_time, 
                labels = c("A", "B"), common.legend = TRUE, legend = "bottom",
                ncol = 2, nrow = 1))
dev.off()

feox_otu_TFS <- filter(otu.tax.meta.fe, wshd == "TFS" &
                     genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella") &
                       date != 20190622 &
                       date != 20190713) %>%
  group_by(sample, otu) %>%
  summarize(sum = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group")) %>%
  filter(sum > 0.01)

feox_TFS_otu_time <- ggplot(feox_otu_TFS) +
  geom_boxplot(aes(x=factor(rdate), y=sum), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=factor(rdate), y=sum, color=otu)) +
  scale_y_log10() + 
  labs(x="Sample Date", y="Relative Abundance (%)", color="OTU") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.key=element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(color="black", size=15))

pdf(file="feox_success_boxplot_otu.pdf", height = 5, width = 9)
print(ggarrange(feox_OKS_otu_time, feox_TFS_otu_time, 
                labels = c("A", "B"), common.legend = TRUE, legend = "bottom",
                ncol = 2, nrow = 1))
dev.off()


# Figure 4
otu.tax.meta.fe <- filter(otu.tax.meta, iron_abun != "Low Fe")

feox.fe <- filter(otu.tax.meta.fe, genus%in%c("Sideroxydans", "Leptothrix", "Gallionellaceae_unclassified", "Gallionella")) %>%
  filter(., rel_abund>0) %>%
  group_by(sample, genus) %>%
  summarize(sum=sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

feob_groups <- ggplot(feox.fe) +
  geom_boxplot(aes(x=genus, y=sum, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=genus, y=sum, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("red", "blue", "black"),
                     breaks=c("Pond","Seep", "Wet Sedge"),
                     labels=c("Pond","Seep", "Wet Sedge")) +
  scale_x_discrete(limits=c("Sideroxydans", "Leptothrix", "Gallionella", "Gallionellaceae_unclassified"),
                   labels=c("Sideroxydans", "Leptothrix", "Gallionella", "Gallionellaceae")) +
  scale_y_log10() + 
  labs(x=NULL, y="Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 30, hjust = 1, color="black", size=13, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black", size=13),
        axis.title.y = element_text(color = "black", size = 13),
        legend.key=element_blank())

ferb.fe <- filter(otu.tax.meta.fe, genus%in%c("Geobacter", "Geobacteraceae_unclassified", "Geotalea", "Citrifermentans", "Shewanella", 
                                      "Anaeromyxobacter", "Pelobacter", "Desulfuromonas", "Carboxydothermus", "Geothrix")) %>%
  filter(., rel_abund>0) %>%
  group_by(sample, genus) %>%
  summarize(sum=sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

ferb_groups <- ggplot(ferb.fe) +
  geom_boxplot(aes(x=genus, y=sum, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=genus, y=sum, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("Geobacter", "Geobacteraceae_unclassified", "Geothrix", "Anaeromyxobacter", 
                            "Citrifermentans", "Desulfuromonas"),
                   labels=c("Geobacter", "Geobacteraceae_unclassified", "Geothrix", "Anaeromyxobacter", 
                            "Citrifermentans", "Desulfuromonas")) +
  scale_y_log10() + 
  labs(x=NULL, y="Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 30, hjust = 1, color="black", size=13, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black", size=13),
        axis.title.y = element_text(color = "black", size = 13),
        legend.key=element_blank())

methox_groups <- ggplot(methox_genus %>% filter(iron_abun == "High Fe")) +
  geom_boxplot(aes(x=family, y=rel_abund, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=family, y=rel_abund, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("Methylomonadaceae", "Methylophilaceae", "Beijerinckiaceae", "Methylococcaceae"),
                   labels=c("Methylomonadaceae", "Methylophilaceae", "Beijerinckiaceae", "Methylococcaceae")) +
  scale_y_log10() + 
  labs(x=NULL, y="Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 30, hjust = 1, color="black", size=13, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black", size=13),
        axis.title.y = element_text(color = "black", size = 13),
        legend.key=element_blank())

ggplot(methox_genus %>% filter(iron_abun == "Low Fe")) +
  geom_boxplot(aes(x=family, y=rel_abund, color=feature_type), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=family, y=rel_abund, color=feature_type)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  scale_x_discrete(limits=c("Methylomonadaceae", "Methylophilaceae", "Beijerinckiaceae", "Methylococcaceae"),
                   labels=c("Methylomonadaceae", "Methylophilaceae", "Beijerinckiaceae", "Methylococcaceae")) +
  scale_y_log10() + 
  labs(x=NULL, y="Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 30, hjust = 1, color="black", size=13, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black", size=13),
        axis.title.y = element_text(color = "black", size = 13),
        legend.key=element_blank())

## Archeal communities in iron mats
arch_class_fe <- otu.tax.meta %>%
  filter(kingdom == "Archaea" & iron_abun == "High Fe") %>%
  group_by(sample, class) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  filter(rel_abund > 0.01) %>%
  inner_join(., metadata, by=c("sample" = "group"))

arch_class_boxplot <- ggplot(arch_class_fe) +
  geom_boxplot(aes(x=reorder(class, -rel_abund), y=rel_abund)) +
  geom_jitter(aes(x=reorder(class, -rel_abund), y=rel_abund, color=feature_type), shape=19, size=2, 
              position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  scale_color_manual(name=NULL,
                     values=c("red", "blue", "black"),
                     breaks=c("Pond", "Seep", "Wet Sedge"),
                     labels=c("Pond", "Seep", "Wet Sedge")) +
  scale_y_log10() +
  labs(x=NULL, y="Relative Abundance (%)") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 30, hjust = 1, color="black", size=13, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black", size=13),
        axis.title.y = element_text(color = "black", size = 13),
        legend.key=element_blank())

pdf(file="Figure4.pdf", height = 16, width = 6)
print(ggarrange(feob_groups, ferb_groups, methox_groups, arch_class_boxplot, 
                labels = c("A", "B", "C", "D"), common.legend = TRUE, legend = "bottom",
                ncol = 1, nrow = 4))
dev.off()

# Organize otu.tax.meta to be by family of interest, then sum all those otus with taxonomic strings containing
# the family of interest and join with the metadata file
gal<-filter(otu.tax.meta, family%in%"Gallionellaceae") %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

geo<-filter(otu.tax.meta, family%in%"Geobacteraceae") %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

ferb<-filter(otu.tax.meta, genus%in%c("Geobacter", "Geobacteraceae_uncultured", "Geotalea", "Citrifermentans", "Shewanella", 
                                      "Anaeromyxobacter", "Pelobacter", "Desulfuromonas", "Carboxydothermus", "Geothrix")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

lepto<-filter(otu.tax.meta, genus%in%"Leptothrix") %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

summary(lepto.nofe)

# used Dedysh, S. N., & Knief, C. (2018). Diversity and Phylogeny of Described Aerobic Methanotrophs. 
# Methane Biocatalysis: Paving the Way to Sustainability, 17–42. and
# ISME Journal (2019) 13:1209–1225 that shows Methylomonadaceae is a methanotroph
methyl<-filter(otu.tax.meta, family%in%c("Methylobacteriaceae", "Methylophilaceae", "Methylomonadaceae", "Methylococcaceae", "Methylothermaceae", 
                                         "Methylocystaceae", "Beijerinckiaceae", "Methylacidophilaceae")) %>% 
  group_by(sample) %>%
  filter(rel_abund>0) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))


methyl_class<-filter(otu.tax.meta, family%in%c("Methylobacteriaceae", "Methylophilaceae", "Methylomonadaceae", "Methylococcaceae", "Methylothermaceae", 
                                         "Methylocystaceae", "Beijerinckiaceae", "Methylacidophilaceae")) %>% 
  group_by(sample, class) %>%
  filter(rel_abund>0) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

methox_class<-filter(otu.tax.meta, genus%in%c("Methylomonadaceae_unclassified", 
                                              "Methylophilaceae_unclassified","Crenothrix", 
                                              "Methylococcus", "Methylomonas", "Methylobacter", 
                                              "Methylomicrobium", "Methylosarcina", "Methylocaldum",
                                              "Methylogaea", "Methylogaea", "Methylosoma", "Methyloparacoccus",
                                               "Methylocaldum", "Methyloglobulus",
                                              "Methyloprofundus", "Methyloprofundus", "Methylomarinum",
                                              "Methylovulum", "Methylomagnum", "Methylosphaera",
                                              "Methylothermus", "Methylohalobius", "Methylomarinovum",
                                               "Methylosinus", "Methylocystis", "Methylocella", 
                                              "Methylocapsa", "Methyloferula", "Methylobacterium-Methylorubrum")) %>% 
  group_by(sample, class) %>%
  filter(rel_abund>0) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

methox_genus<-filter(otu.tax.meta, genus%in%c("Methylomonadaceae_unclassified", 
                                              "Methylophilaceae_unclassified","Crenothrix", 
                                              "Methylococcus", "Methylomonas", "Methylobacter", 
                                              "Methylomicrobium", "Methylosarcina", "Methylocaldum",
                                              "Methylogaea", "Methylogaea", "Methylosoma", "Methyloparacoccus",
                                              "Methylocaldum", "Methyloglobulus",
                                              "Methyloprofundus", "Methyloprofundus", "Methylomarinum",
                                              "Methylovulum", "Methylomagnum", "Methylosphaera",
                                              "Methylothermus", "Methylohalobius", "Methylomarinovum",
                                              "Methylosinus", "Methylocystis", "Methylocella", 
                                              "Methylocapsa", "Methyloferula", "Methylobacterium-Methylorubrum")) %>% 
  group_by(sample, family) %>%
  filter(rel_abund>0) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

methox_class_plot <- ggplot(methox_class) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund, shape=class), outlier.shape = NA) +
  geom_jitter(aes(x=iron_abun, y=rel_abund, shape=class, color=methane_nM), size=2, 
              position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  scale_colour_gradientn(colours = c("blue", "yellow", "red")) +
  xlab("Iron Abundance") +
  ylab("Relative Abundance (%)") +
  labs(shape="Class", color="Methane (nM)") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(color="black", size=13),
        axis.title.x = element_text(color = "black", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black", size=13),
        axis.title.y = element_text(color = "black", size = 13),
        legend.key=element_blank())

pdf("Figure7.pdf")
print(methox_class_plot)
dev.off()


# using Holmes and Smith 2016 Advances in Applied Microbiology section 2.2 Phylogeny of Methanogens to determine methanogens families
methano<-filter(otu.tax.meta, order%in%c("Methanobacteriales", "Methanococcales", "Methanomicrobiales", "Methanosarcinales", "Methanopyrales")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

summary(methano.nofe)

methano<-filter(otu.abun.tax, order%in%c("Methanobacteriales", "Methanococcales", "Methanomicrobiales", "Methanosarcinales", "Methanopyrales")) %>% 
  group_by(sample) %>%
  summarise(abun = sum(abun)) %>%
  filter(abun>1) %>%
  inner_join(., metadata, by=c("sample" = "group"))

nophotocyano<-filter(otu.tax.meta, class%in%c("Melainabacter", "Sericytochromatia")) %>% 
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))



## sum of Fe(II)-oxidizing taxa for high and low Fe sites for comparison with
## Remi's remote sensing analysis
gallepto <- mutate(gal, sum=gal$rel_abund+lepto$rel_abund)
write.table(gallepto, "~/Documents/Toolm/gallepto.txt", sep = "\t")

# make box plots with dots behind to show where there is unequal sampling
# combine into a multipanel figure
## Figure 6
## box plots by Fe abundance
gal_fe<-ggplot(gal) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund), outlier.shape = NA) +
  geom_jitter(aes(x=iron_abun, y=rel_abund), shape=19, size=2)  +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  scale_y_log10() + 
  labs(x=NULL, y="Relative Abundance (%)") +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"))

lepto_fe<-ggplot(lepto) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund), outlier.shape = NA) +
  geom_jitter(aes(x=iron_abun, y=rel_abund), shape=19, size=2)  +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  scale_y_log10() + 
  labs(x=NULL, y=NULL) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"))

ferb_fe<-ggplot(ferb) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund), outlier.shape = NA) +
  geom_jitter(aes(x=iron_abun, y=rel_abund), shape=19, size=2)  +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  scale_y_log10() + 
  labs(x=NULL, y=NULL) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"))

methox_sample<-filter(otu.tax.meta, genus%in%c("Methylomonadaceae_unclassified", 
                                              "Methylophilaceae_unclassified","Crenothrix", 
                                              "Methylococcus", "Methylomonas", "Methylobacter", 
                                              "Methylomicrobium", "Methylosarcina", "Methylocaldum",
                                              "Methylogaea", "Methylogaea", "Methylosoma", "Methyloparacoccus",
                                              "Methylocaldum", "Methyloglobulus",
                                              "Methyloprofundus", "Methyloprofundus", "Methylomarinum",
                                              "Methylovulum", "Methylomagnum", "Methylosphaera",
                                              "Methylothermus", "Methylohalobius", "Methylomarinovum",
                                              "Methylosinus", "Methylocystis", "Methylocella", 
                                              "Methylocapsa", "Methyloferula", "Methylobacterium-Methylorubrum")) %>% 
  group_by(sample) %>%
  filter(rel_abund>0) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "group"))

methyl_fe<-ggplot(methox_sample) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund), outlier.shape = NA) +
  geom_jitter(aes(x=iron_abun, y=rel_abund), shape=19, size=2)  +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  scale_y_log10() + 
  labs(x=NULL, y=NULL) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"))

methano_fe<-ggplot(methano) +
  geom_boxplot(aes(x=iron_abun, y=rel_abund), outlier.shape = NA) +
  geom_jitter(aes(x=iron_abun, y=rel_abund), shape=19, size=2)  +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  scale_y_log10() + 
  labs(title="Methanogens", x="Fe Presence", y=NULL) +
  theme_classic()

pdf(file="Figure6.pdf", height = 5, width = 15)
print(ggarrange(gal_fe, lepto_fe, ferb_fe, methyl_fe, 
                labels = c("A", "B", "C", "D"),
                ncol = 4, nrow = 1))
dev.off()

ggplot(metadata) +
  geom_boxplot(aes(x=iron_abun, y=totfe3_ext), outlier.shape = NA) +
  geom_jitter(aes(x=iron_abun, y=totfe3_ext), shape=19, size=2)  +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  labs(x="Fe Presence",
       y="Total Extractable Fe(III) (umol/gdw)") +
  theme_classic()

ggplot(metadata) +
  geom_boxplot(aes(x=iron_abun, y=methane_nM), outlier.shape = NA) +
  geom_jitter(aes(x=iron_abun, y=methane_nM), shape=19, size=2)  +
  scale_x_discrete(limits=c("High Fe", "Low Fe"),
                   labels=c("High Fe", "Low Fe")) +
  labs(x="Fe Presence",
       y="Methane (uM)") +
  theme_classic()



# Extractable Fe in high and low Fe sites
meta_outlier <- metadata %>% group_by(iron_abun) %>% identify_outliers("totfe3_ext")
metadata %>% group_by(iron_abun) %>% shapiro_test(totfe3_ext)
ggqqplot(metadata, x = "totfe3_ext", facet.by = "iron_abun")
metadata %>% levene_test(totfe3_ext ~ factor(iron_abun))

fe_test <- metadata %>% t_test(totfe3_ext ~ iron_abun) %>% add_significance() %>% add_xy_position(x = "iron_abun")

fe_bxp <- ggboxplot(metadata, x = "iron_abun", y = "totfe3_ext", 
                         ylab = "Total Extractable Fe(III)", add = "jitter") +
  stat_pvalue_manual(fe_test, tip.length = 0) + labs(subtitle = get_test_label(fe_test, detailed = FALSE)) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.title.x = element_blank())

#Independent groups from distinct populations and the sample don't affect each other.
#Does not need to be normally distributed data
wilcox.test(totfe3_ext ~ iron_abun, data = metadata)
#Null hyp is that methane concs of high and low iron sites are identical
#p-value is 0.6, and is greater than the 0.05 significance level, we accept the null hyp.
#same can be said about the t test analysis from above, but the data defy the normality and outlier test.

# methane in high and low Fe sites
meta_outlier <- metadata %>% group_by(iron_abun) %>% identify_outliers(methane_nM)
metadata %>% group_by(iron_abun) %>% shapiro_test(methane_nM)
ggqqplot(metadata, x = "methane_nM", facet.by = "iron_abun")
metadata %>% levene_test(methane_nM ~ factor(iron_abun))

methane_test <- metadata %>% t_test(methane_nM ~ iron_abun) %>% add_significance() %>% add_xy_position(x = "iron_abun")

methane_bxp <- ggboxplot(metadata, x = "iron_abun", y = "methane_nM", 
                         ylab = "Methane (nM)", add = "jitter") +
  stat_pvalue_manual(methane_test, tip.length = 0) + labs(subtitle = get_test_label(methane_test, detailed = FALSE)) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.title.x = element_blank())

#Independent groups from distinct populations and the sample don't affect each other.
#Does not need to be normally distributed data
wilcox.test(methane_nM ~ iron_abun, data = metadata)
#Null hyp is that methane concs of high and low iron sites are identical
#p-value is 0.6, and is greater than the 0.05 significance level, we accept the null hyp.
#same can be said about the t test analysis from above, but the data defy the normality and outlier test.

# pH in high and low Fe sites
meta_outlier <- metadata %>% group_by(iron_abun) %>% identify_outliers(pH)
metadata %>% group_by(iron_abun) %>% shapiro_test(pH)
ggqqplot(metadata, x = "pH", facet.by = "iron_abun")
metadata %>% levene_test(pH ~ factor(iron_abun))

pH_test <- metadata %>% t_test(pH ~ iron_abun) %>% add_significance() %>% add_xy_position(x = "iron_abun")

pH_bxp <- ggboxplot(metadata, x = "iron_abun", y = "pH", 
                    ylab = "pH", add = "jitter") +
  stat_pvalue_manual(pH_test, tip.length = 0) + labs(subtitle = get_test_label(pH_test, detailed = FALSE)) +
  theme(text = element_text(size=18),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.title.x = element_blank())

#Independent groups from distinct populations and the sample don't affect each other.
#Does not need to be normally distributed data
wilcox.test(pH ~ iron_abun, data = metadata)
#Null hyp is that pH values of high and low iron sites are identical
#p-value is <0.00001, and is less than the 0.05 significance level, we reject the null hyp.
#same can be said about the t test analysis from above, but the data defy the normality and outlier test.

pdf(file="geochem_stats.pdf", height = 4, width = 12)
print(ggarrange(fe_bxp, methane_bxp, pH_bxp, 
                labels = c("A", "B", "C"),
                ncol = 3, nrow = 1))
dev.off()

# sig taxa rel abund with geochem
gal_methane<-ggscatter(feox, x="sum", y="methane_nM", add = "reg.line", conf.int = TRUE, color = "iron_abun", palette = "Dark2") +
  scale_alpha_manual(name=NULL,
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(title="Fe-Oxidizing Bacteria",
       x="Relative Abundance (%)",
       y="Methane (nM)") +
  scale_y_log10() +
  stat_cor(aes(color=iron_abun), label.x = 3, label.y = c(1.6, 1.9)) +
  theme_classic()

lepto_methane<-ggscatter(lepto, x="rel_abund", y="methane_nM", add = "reg.line", conf.int = TRUE, color = "iron_abun", palette = "Dark2") +
  scale_alpha_manual(name=NULL,
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(title="Leptothrix",
       x="Relative Abundance (%)",
       y=NULL) +
  scale_y_log10() +
  stat_cor(aes(color=iron_abun), label.x = 0.5, label.y = c(1.6, 1.9)) +
  theme_classic()

ferb_methane<-ggscatter(ferb, x="rel_abund", y="methane_nM", add = "reg.line", conf.int = TRUE, color = "iron_abun", palette = "Dark2") +
  scale_alpha_manual(name=NULL,
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(title="Fe-Reducing Bacteria",
       x="Relative Abundance (%)",
       y=NULL) +
  scale_y_log10() +
  stat_cor(aes(color=iron_abun), label.x = 4, label.y = c(1.6, 1.9)) +
  theme_classic()

methyl_methane<-ggscatter(methox, x="sum", y="methane_nM", add = "reg.line", conf.int = TRUE, color = "iron_abun", palette = "Dark2") +
  scale_alpha_manual(name=NULL,
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(title="Methylotrophs",
       x="Relative Abundance (%)",
       y="Methane (nM)") +
  scale_y_log10() +
  stat_cor(aes(color=iron_abun), label.x = 0.8, label.y = c(1.6, 1.9)) +
  theme_classic()

methano_methane<-ggscatter(methano, x="rel_abund", y="methane_nM", add = "reg.line", conf.int = TRUE, color = "iron_abun", palette = "Dark2") +
  scale_alpha_manual(name=NULL,
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(title="Methanogens",
       x="Relative Abundance (%)",
       y=NULL) +
  scale_y_log10() +
  stat_cor(aes(color=iron_abun), label.x = 0.02, label.y = c(1.6, 1.9)) +
  theme_classic()

pdf(file="sigtaxa_methane.pdf", height = 3, width = 15)
print(ggarrange(gal_methane, lepto_methane, ferb_methane, methyl_methane, methano_methane, 
                labels = c("A", "B", "C", "D", "E"), common.legend = TRUE, legend = "bottom",
                ncol = 5, nrow = 1))
dev.off()

# assemble all of the summed rel abunds from each sig taxa into single dataframe
gal_x<-gal[,1:2] %>%
  rename(rel_abund_gal = rel_abund)
ferb_x<-ferb[,1:2] %>%
  rename(rel_abund_ferb = rel_abund)
lepto_x<-lepto[,1:2] %>%
  rename(rel_abund_lepto = rel_abund)
methyl_x<-methyl[,1:2] %>%
  rename(rel_abund_methyl = rel_abund)
methano_x<-methano[,1:2] %>%
  rename(rel_abund_methano = rel_abund)
cyano_x<-cyano[,1:2] %>%
  rename(rel_abund_cyano = rel_abund)
#thrix_x<-geothrix[,1:2] %>%
  #rename(rel_abund_geothrix = rel_abund)
gal.ferb<-inner_join(gal_x, ferb_x, by="sample")
gal.ferb.lepto<- inner_join(gal.ferb, lepto_x, by="sample")
#gal.ferb.lepto.thrix <- inner_join(gal.ferb.lepto, thrix_x, by="sample")
gal.ferb.lepto.methyl<-inner_join(gal.ferb.lepto, methyl_x, by="sample")
gal.ferb.lepto.methyl.methano<-inner_join(gal.ferb.lepto.methyl, methano_x, by="sample")
gal.ferb.lepto.methyl.methano.cyano<-inner_join(gal.ferb.lepto.methyl.methano, cyano_x, by="sample")
sig_taxa.meta<-inner_join(gal.ferb.lepto.methyl.methano.cyano, metadata, by=c("sample" = "group"))
sig_taxa.meta.feonly <- filter(sig_taxa.meta, iron_abun%in%"High Fe")

FeO <- mutate(gal, feosum=gal$rel_abund+lepto$rel_abund)
FeOR <- inner_join(FeO, ferb, by="sample")


# Gallionella and methanotroph correlations
gal.methyl.fe<-ggscatter(sig_taxa.meta, x="rel_abund_gal", y="rel_abund_methyl", add="reg.line", conf.int=TRUE, color = "iron_abun") +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black"),
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(x="Gallionellaceae Relative Abundance (%)",
       y="Methanotrophs Relative Abundance (%)") +
  theme_classic() +
  stat_cor(aes(color=iron_abun), label.x = 1)

gal.methyl.feature<-ggscatter(sig_taxa.meta.feonly, x="rel_abund_gal", y="rel_abund_methyl", add="reg.line", conf.int=TRUE, color = "feature_type") +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black", "blue"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  labs(x="Gallionellaceae Relative Abundance (%)",
       y=NULL) +
  theme_classic() +
  stat_cor(aes(color=feature_type), label.x = 1)

pdf("gal_methyl_relationships.pdf", height = 4, width = 12)
ggarrange(gal.methyl.fe, gal.methyl.feature, labels=c("A", "B"), nrow = 1, ncol = 2)
dev.off()

ggscatter(sig_taxa.meta.feonly, x="rel_abund_cyano", y="rel_abund_ferb", color = "feature_type", add="reg.line", conf.int=TRUE) +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black", "blue"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  labs(x="Relative Abundance of Cyanobacteria (%)", y="Relative Abundance of Fe-Reducing Bacteria (%)") +
  theme_classic() +
  stat_cor(aes(color=feature_type), label.x = 1)

feo.fer.feature<-ggscatter(sig_taxa.meta, x="rel_abund_gal", y="rel_abund_ferb", color = "feature_type", add="reg.line", conf.int=TRUE) +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black", "blue"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  labs(x="Relative Abundance of Fe-Oxidizing Bacteria (%)", y="Relative Abundance of Fe-Reducing Bacteria (%)") +
  theme_classic() +
  stat_cor(aes(color=feature_type), label.x = 1)

ferb.methano.feature<-ggscatter(sig_taxa.meta, x="rel_abund_ferb", y="rel_abund_methano", color = "feature_type", add="reg.line", conf.int=TRUE) +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black", "blue"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  labs(x="Relative Abundance of Fe-Reducing Bacteria (%)", y="Relative Abundance of Methanogens (%)") +
  theme_classic() +
  stat_cor(aes(color=feature_type), label.x = 1)

feo.fer.feature.feonly<-ggscatter(sig_taxa.meta.feonly, x="rel_abund_ferb", y="rel_abund_methano", color = "feature_type", add="reg.line", conf.int=TRUE) +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black", "blue"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  labs(x="Relative Abundance of Fe-Oxidizing Bacteria (%)", y="Relative Abundance of Fe-Reducing Bacteria (%)") +
  theme_classic() +
  stat_cor(aes(color=feature_type), label.x = 1)

ggscatter(sig_taxa.meta.feonly, x="rel_abund_gal", y="rel_abund_ferb", color = "feature_type", add="reg.line", conf.int=TRUE) +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black", "blue"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  labs(x="Relative Abundance of Fe-Oxidizing Bacteria (%)", y="Relative Abundance of Fe-Reducing Bacteria (%)") +
  theme_classic() +
  stat_cor(aes(color=feature_type), label.x = 1)

# Leptothrix and methanotroph correlations
lepto.methyl.fe<-ggscatter(sig_taxa.meta, x="rel_abund_lepto", y="rel_abund_methyl", add="reg.line", conf.int=TRUE, color = "iron_abun") +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black"),
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(x="Leptothrix Relative Abundance (%)",
       y="Methanotrophs Relative Abundance (%)") +
  theme_classic() +
  stat_cor(aes(color=iron_abun), label.x = 0.5)

lepto.methyl.feature<-ggscatter(sig_taxa.meta.feonly, x="rel_abund_lepto", y="rel_abund_methyl", add="reg.line", conf.int=TRUE, color = "feature_type") +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black", "blue"),
                     breaks=c("Seep", "Pond", "Wet Sedge"),
                     labels=c("Seep", "Pond", "Wet Sedge")) +
  labs(x="Leptothrix Relative Abundance (%)",
       y=NULL) +
  theme_classic() +
  stat_cor(aes(color=feature_type), label.x = 0.5)

pdf("lepto_methyl_relationships.pdf", height = 4, width = 12)
ggarrange(lepto.methyl.fe, lepto.methyl.feature, labels=c("A", "B"), nrow = 1, ncol = 2)
dev.off()

## no real relationship betwen O2 and sig taxa rel abund, but they are below.
ggplot(lepto, aes(x=rel_abund, y=oxy_mgL, color=iron_abun)) +
  geom_point(shape=19, size=2) +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black"),
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(title="Leptothrix",
       x="Relative Abundance (%)",
       y="Oxygen (mg L-1)") +
  theme_classic()

ggplot(methano, aes(x=rel_abund, y=oxy_mgL, color=iron_abun)) +
  geom_point(shape=19, size=2) +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black"),
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(title="Methanogens",
       x="Relative Abundance (%)",
       y="Oxygen (mg L-1)") +
  theme_classic()

ggplot(methyl, aes(x=rel_abund, y=oxy_mgL, color=iron_abun)) +
  geom_point(shape=19, size=2) +
  scale_alpha_manual(name=NULL,
                     values=c("red", "black"),
                     breaks=c("High Fe", "Low Fe"),
                     labels=c("High Fe", "Low Fe")) +
  labs(title="Methanotrophs",
       x="Relative Abundance (%)",
       y="Oxygen (mg L-1)") +
  theme_classic()

toolm.seqs<-read_excel(path="toolm.neqs.xlsx")
toolm.adiv<-read_excel(path="toolm.adiv.xlsx")
env.data<-read_excel(path = "toolm_all_geochem.xlsx")

toolm.adiv2<-toolm.adiv[-c(69:163),]
toolm.adiv.meta<-inner_join(metadata, toolm.adiv2, by=c('group'='group'))

## PCA of environmental variables
# Toolik data
?color_palette
env.data<-read_excel(path = "toolm_all_geochem.xlsx")
oks.cols<-brewer.pal(6, "Set2")
tfs.cols<-brewer.pal(8, "Set1")
pca.cols2<-c(oks.cols, tfs.cols)
scaleredblue <- colorRampPalette(c("red", "blue"), space = "rgb")(6)
scaleyelgreen <- colorRampPalette(c("yellow", "green"), space = "rgb")(8)
pca.cols<-c(scaleredblue, scaleyelgreen)
wshd.site <- env.data$wshd_site
site.data <- env.data[,2:7]
ch4<-((site.data$ch4-mean(site.data$ch4, na.rm = T))/sd(site.data$ch4, na.rm = T))
fe2<-((site.data$water_fe2-mean(site.data$water_fe2, na.rm = T))/sd(site.data$water_fe2, na.rm = T))
temp<-((site.data$temp-mean(site.data$temp, na.rm = T))/sd(site.data$temp, na.rm = T))
o2<-((site.data$o2-mean(site.data$o2, na.rm = T))/sd(site.data$o2, na.rm = T))
s_cond<-((site.data$s_cond-mean(site.data$s_cond, na.rm = T))/sd(site.data$s_cond, na.rm = T))
pH<-((site.data$pH-mean(site.data$pH, na.rm = T))/sd(site.data$pH, na.rm = T))
site.data.norm <- data.frame(ch4, fe2, temp, o2, s_cond, pH)

site.data.norm2<-scale(site.data)

site.data.colnames <- c("ch4", "water_fe2", "temp", "o2", "s_cond", "pH")
names(site.data.norm) <- site.data.colnames
?plot()

toolm.env.pca2<-princomp(na.omit(site.data.norm2), cor = FALSE)
pdf(file = "Env_pca.pdf", height = 8, width = 11)
p <- plot(toolm.env.pca2$scores, xlim = c(-3,7), pch = c(16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17), cex=2, col = pca.cols2)
legend(5.1, 3, legend = c("OKS_BSeep2", "OKS_Pond2", "OKS_Pond4", "OKS_Pond6", "OKS_Pond9", "OKS_WSedge1", "TFS_Pond1", "TFS_Pond10", "TFS_Pond11", "TFS_Pond2", "TFS_Pond3", "TFS_Pond4", "TFS_WSedge1", "TFS_WSedge10"), 
       col = pca.cols2, pch = c(16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17), bty='n', cex=1)
print(p)
dev.off()
dim(env.data)
dim(na.omit(env.data))

ggplot(metadata) +
  geom_point(aes(x=rdate, y=water_fe2_uM, color=feature_type)) +
  scale_y_log10()

rain <- read_tsv("tfs_rain.txt")
ggplot(rain) +
  geom_line(aes(x=date, y=cum_precip))
