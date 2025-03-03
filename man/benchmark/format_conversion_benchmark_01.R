# prepare the test data (download processed objects from CELLxGENE)
library(tidyverse)
library(GEfetch2R)
# all available datasets
all.cellxgene.datasets = readRDS("all.cellxgene.datasets.rds")
# one million
cellxgene.oneM.datasets = all.cellxgene.datasets[all.cellxgene.datasets$cell_count>1000000, ] %>%
  dplyr::arrange(cell_count)
# select datasets (Supplementary Table S1)
cellxgene.oneM.datasets.used = cellxgene.oneM.datasets[c(1, 2, 3, 9, 26),]
# openxlsx::write.xlsx(x = cellxgene.oneM.datasets.used, file = "supp_s1.xlsx")
# download data
cellxgene.oneM.down = ParseCELLxGENE(meta = cellxgene.oneM.datasets.used, file.ext = c("rds","h5ad"),
                                     out.folder = "/Volumes/soyabean/oneM", parallel = FALSE,
                                     timeout = 3600000000000)
