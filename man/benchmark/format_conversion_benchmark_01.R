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


#### --------------------------------------
# set up
base.folder = '/home/syb/projects/04_scfetch/benchtime'
data.folder = file.path(base.folder, "oneM")
out.folder = file.path(base.folder, "runtimeOut")
script.folder =  file.path(base.folder, "runScripts")
sce.folder = file.path(base.folder, "oneMsce")
# subsample the downloaded data (on a Linux Server)
all.rds = list.files(path = data.folder, pattern = "rds$", full.names = T)
sub.vec = c(1200000, 1000000, 800000, 600000, 500000, 300000, 200000, 100000, 80000, 60000, 50000, 30000, 20000, 10000)
for (i in 1:length(all.rds)){
  seu = readRDS(all.rds[i])
  base.name = gsub(pattern = ".rds$", replacement = "", x = basename(all.rds[i]))
  base.name = gsub(pattern = "\\.+", replacement = "_", x = base.name)
  for (cell in sub.vec){
    out.folder.name = paste0("C", as.character(cell/10000), 'w')
    out.file.name = paste0(base.name, "_", out.folder.name, ".rds")
    dir.create(file.path(data.folder, out.folder.name), showWarnings = FALSE)
    if(ncol(seu) >= cell){
      seu.subsampled <- seu[, sample(colnames(seu), size =cell, replace=F)]
      saveRDS(object = seu.subsampled, file = file.path(data.folder, out.folder.name, out.file.name))
    }
  }
}

