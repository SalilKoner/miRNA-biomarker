# Need imputation or not
do.impute <- FALSE

source("Source/analysis.R")

#****************************************************************************************************
#                          Table 1: Significantly deregulated miRNA's                               *
#****************************************************************************************************

print(discoveries)

# Converting to latex output
my_table <- discoveries[,-c(4, 7, 9)]
names(my_table) <- c("MiRNA", "chr.", "strand", "p-value", "missingness", "Significant tests")
print(xtable(my_table), include.rownames=FALSE)


#****************************************************************************************************
#                         Figure 1: Heat map of the double delta CT values                          *
#****************************************************************************************************


deltaCT_dat <- mydata
names(deltaCT_dat) <- c(names(mydata)[1:5], paste0("Subj", 1:18, "_C"), paste0("Subj", 1:18, "_N"))

deltaCT_dat = deltaCT_dat %>% pivot_longer(
  cols = starts_with("Subj"),
  names_prefix = "Subj",
  names_sep = "_",
  names_to = c("Subject", "Gene_type")
) %>%
  mutate(Subject = as.numeric(Subject)) %>%
  group_by(MiR, Chr, strand, start_coordinate, centro, Subject) %>%
  summarize(deltaCT = value[Gene_type == "C"] - value[Gene_type == "N"], .groups = "drop") %>%
  pivot_wider(names_from = "Subject", values_from = "deltaCT", names_prefix = "Subj")


plotdat <- deltaCT_dat %>% inner_join(discoveries) %>% arrange(desc(pval))

delCT.matrix <- plotdat %>% dplyr::select(starts_with("Subj")) %>% as.matrix
rownames(delCT.matrix) <- plotdat$MiR
colnames(delCT.matrix) <- paste0('ID ', 1:18)

# Plot presented in Figure 1
superheat::superheat(delCT.matrix,
                     # set heatmap color map
                     heat.pal = brewer.pal(5, "RdBu"),
                     heat.na.col = "grey50",
                     yr = plotdat$pval,
                     yr.plot.type = "bar",
                     yr.axis.name = "Unadjusted p-values",
                     yr.plot.size = 0.85,
                     yr.num.ticks = 5,
                     yr.axis.size = 20,
                     yr.axis.name.size = 17,
                     yr.bar.col = "white",
                     bottom.label.text.size = 7,
                     bottom.label.size = 0.27,
                     bottom.label.col = "white",
                     bottom.label.text.angle = 90,
                     bottom.label.text.alignment = "right",
                     padding = 0.7,
                     legend.height = 0.09,
                     legend.width = 3,
                     legend.vspace = 0.01)



#****************************************************************************************************
#     Table S1: size, frequency, chromosome number, and the corresponding strand for each group.    *
#****************************************************************************************************

group_hist <- mydata %>% inner_join(group.info, by = c("Chr", "strand", "centro")) %>%
  select(MiR, group_id, Chr, strand, centro, start_coordinate)

groups_dist <- group_hist %>% group_by(group_id) %>%
  dplyr::summarise(distance = round({max(start_coordinate) - min(start_coordinate)}/
                                      1e6, 4), count = n())

print(groups_dist)

print(xtable(groups_dist), include.rownames=FALSE)


#****************************************************************************************************
#               Table S2: Significantly deregulated miRNA's using imputed data                      *
#****************************************************************************************************

rm(list=ls())

do.impute <- TRUE

source("Source/analysis.R")

print(discoveries)

# Converting to latex output
my_table <- discoveries[,-c(4, 7, 9)]
names(my_table) <- c("MiRNA", "chr.", "strand", "p-value", "missingness", "Significant tests")
print(xtable(my_table), include.rownames=FALSE)



