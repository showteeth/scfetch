# process the output

# library
library(tidyverse)
library(Seurat)

# read the output
setwd('/home/syb/data/projects/scfetch/benchtime/runtimeOut')
all_runtime_file = list.files(path = "./", pattern = "runtime.rds$", recursive = T)
all_runtime_df = data.frame()
for (rf in all_runtime_file){
  rf.name = basename(rf)
  rf.split = strsplit(x = rf.name, split = "_")[[1]]
  rf.conversion = rf.split[1]
  rf.soft = rf.split[2]
  rf.cellnumber = rev(rf.split)[2]
  rf.dataset = gsub(pattern = paste(rf.conversion, rf.soft, "(.*?)", rf.cellnumber, "runtime.rds", sep = "_"),
                    replacement = "\\1", x = rf.name)
  rf.dataset = gsub(pattern = "_+$", replacement = "", x = rf.dataset)
  rf.time = readRDS(rf)
  rf.time.user = rf.time[1]
  rf.time.system = rf.time[2]
  rf.time.elapsed = rf.time[3]
  rf.df = data.frame(rf.conversion, rf.soft, rf.dataset, rf.cellnumber, rf.time.user, rf.time.system, rf.time.elapsed)
  all_runtime_df = rbind(all_runtime_df, rf.df)
}
rownames(all_runtime_df) = NULL
all_runtime_df = all_runtime_df %>%
  dplyr::arrange(rf.conversion, rf.cellnumber, rf.dataset, rf.time.elapsed)

# add abbreviation
all_runtime_df = all_runtime_df %>%
  dplyr::mutate(dataset.abbr = ifelse(rf.dataset == "AIDA_Phase_1_Data_Freeze_v1_Chinese_Indian_Japanese_Korean_and_Malay_donors_in_Japan_Singapore_and_South_Korea", "AIDA_Phase1",
                                      ifelse(rf.dataset == "A_single_cell_multi_omic_atlas_spanning_the_adult_rhesus_macaque_brain", "Mmu_brain",
                                             ifelse(rf.dataset == "Major_cell_cluster_Definitive_erythroid", "Definitive_erythroid",
                                                    ifelse(rf.dataset == "Survey_of_human_embryonic_development_1_million_cells_subset", "Embryo_dev (1M)",
                                                           ifelse(rf.dataset == "The_single_cell_lung_cancer_atlas_LuCA_extended_atlas", "LuCA_extended", "empty"))))))
all_runtime_df$rf.cellnumber = gsub(pattern = "^C", replacement = "", x = all_runtime_df$rf.cellnumber)
all_runtime_df$rf.cellnumber = gsub(pattern = "w$", replacement = "0000", x = all_runtime_df$rf.cellnumber)
all_runtime_df$rf.cellnumber = as.numeric(all_runtime_df$rf.cellnumber)
all_runtime_df$rf.cellnumber = scales::label_number(scale_cut = scales::cut_short_scale(),accuracy = 0.1)(all_runtime_df$rf.cellnumber)
all_runtime_df$rf.cellnumber = gsub(pattern = "\\.0", replacement = "", x = all_runtime_df$rf.cellnumber)
all_runtime_df$rf.cellnumber = factor(all_runtime_df$rf.cellnumber,
                                      levels = c("10K", "20K", "30K", "50K", "60K", "80K", "100K", "200K",
                                                 "300K", "500K", "600K", "800K", "1M", "1.2M"))
all_runtime_df_plot = ggplot(all_runtime_df, aes(x=rf.cellnumber, y=rf.time.elapsed, group = rf.soft)) +
  geom_line(aes(color = rf.soft), size=0.5)+
  geom_point(aes(shape=rf.soft, color = rf.soft), size=1.2) +
  # facet_grid(vars(dataset.abbr), vars(rf.conversion), scales = "free") +
  facet_wrap(vars(rf.conversion, dataset.abbr), scales = "free_y", ncol=4, dir="v") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Cell number", y = "Time (s)")
# save the figure
setwd('/home/syb/projects/04_scfetch/benchtime')
# ggsave(filename = "all_runtime_df_plot.pdf", plot = all_runtime_df_plot, height = 9, width = 10)
ggsave(filename = "all_runtime_df_plot.pdf", plot = all_runtime_df_plot, height = 11, width = 11)



