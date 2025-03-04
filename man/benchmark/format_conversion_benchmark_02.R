# generate scripts
# data: format_conversion_benchmark_01.R

# library -----
library(tidyverse)
library(Seurat)

# SeuratObject to AnnData (oneM) --------
cell.sub.vec = c('C120w', 'C100w', 'C80w', 'C60w', 'C50w', 'C30w', 'C20w', 'C10w', 'C8w', 'C6w', 'C5w', 'C3w', 'C2w', 'C1w')
# * 1.SeuratDisk -------
dir.create(file.path(out.folder, 'Seu2AD_SeuratDisk'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'Seu2AD_SeuratDisk'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'Seu2AD_SeuratDisk', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'Seu2AD_SeuratDisk', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.rds = list.files(path = cs.data.folder, pattern = "rds$", full.names = T)
  for (rds in sub.rds){
    filename = basename(rds)
    out.h5seurat = file.path(cs.out.folder, paste0("Seu2AD_SeuratDisk_", gsub(pattern = ".rds$", replacement = ".h5Seurat", x = filename)))
    # second: remove output h5ad file
    out.h5ad = file.path(cs.out.folder, paste0("Seu2AD_SeuratDisk_", gsub(pattern = ".rds$", replacement = ".h5ad", x = filename)))
    cat(paste0("seu=readRDS('",rds,"')"),
        "runtime = system.time({",
        paste0("SeuratDisk::SaveH5Seurat(seu, filename = '",out.h5seurat,"', overwrite = TRUE)"),
        paste0("SeuratDisk::Convert('", out.h5seurat, "', dest = 'h5ad', assay = 'RNA', overwrite = TRUE)"),
        "})",
        paste0("file.remove('",out.h5seurat,"')"),
        # second: remove output h5ad file
        paste0("file.remove('",out.h5ad,"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("Seu2AD_SeuratDisk_", gsub(pattern = ".rds$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("Seu2AD_SeuratDisk_",paste0(gsub(pattern = ".rds$", replacement = ".R", x = filename)))),sep="\n")
  }
}

# example
# system.time({
#   SeuratDisk::SaveH5Seurat(pbmc3k.final, filename = file.path(out.folder, "pbmc3k.h5Seurat"), overwrite = TRUE)
#   SeuratDisk::Convert(file.path(out.folder, "pbmc3k.h5Seurat"), dest = "h5ad", assay = "RNA", overwrite = TRUE)
# })

# * 2.sceasy --------
dir.create(file.path(out.folder, 'Seu2AD_sceasy'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'Seu2AD_sceasy'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'Seu2AD_sceasy', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'Seu2AD_sceasy', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.rds = list.files(path = cs.data.folder, pattern = "rds$", full.names = T)
  for (rds in sub.rds){
    filename = basename(rds)
    out.h5ad = file.path(cs.out.folder, paste0("Seu2AD_sceasy_", gsub(pattern = ".rds$", replacement = ".h5ad", x = filename)))
    cat(paste0("seu=readRDS('",rds,"')"),
        "runtime = system.time({",
        "reticulate::use_condaenv('testsc', required = TRUE)",
        paste0("sceasy::convertFormat(seu, from='seurat', to='anndata', drop_single_values = F, outFile = '",out.h5ad,"', main_layer = 'counts', assay='RNA')"),
        "})",
        # second: remove output h5ad file
        paste0("file.remove('",out.h5ad,"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("Seu2AD_sceasy_", gsub(pattern = ".rds$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("Seu2AD_sceasy_",paste0(gsub(pattern = ".rds$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# system.time({
#   # reticulate::use_condaenv("testsc", required = TRUE)
#   # or set RETICULATE_PYTHON = "/Applications/anaconda3/bin/python" in Renvion
#   sceasy::convertFormat(pbmc3k.final, from="seurat", to="anndata", drop_single_values = F,
#                         outFile=file.path(out.folder, "pbmc3k.h5ad"), main_layer = "counts", assay="RNA")
# })

# * 3.scDIOR ---------
dir.create(file.path(out.folder, 'Seu2AD_scdior'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'Seu2AD_scdior'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'Seu2AD_scdior', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'Seu2AD_scdior', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.rds = list.files(path = cs.data.folder, pattern = "rds$", full.names = T)
  for (rds in sub.rds){
    filename = basename(rds)
    out.h5 = file.path(cs.out.folder, paste0("Seu2AD_scdior_", gsub(pattern = ".rds$", replacement = ".h5", x = filename)))
    cat(paste0("seu=readRDS('",rds,"')"),
        "runtime = system.time({",
        paste0("dior::write_h5(data = seu, object.type='seurat', file = '",out.h5,"', assay.name='RNA', save.scale = F)"),
        "})",
        "# adata = diopy.input.read_h5(file = 'pbmc3k.h5') # require diopy to load h5 to AnnData",
        # second: remove output h5 file
        paste0("file.remove('",out.h5,"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("Seu2AD_scdior_", gsub(pattern = ".rds$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("Seu2AD_scdior_",paste0(gsub(pattern = ".rds$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# system.time({
#   dior::write_h5(data = pbmc3k.final, object.type = "seurat", file = file.path(out.folder, "pbmc3k.h5"),
#                  assay.name = "RNA", save.scale = F)
#   # adata = diopy.input.read_h5(file = 'pbmc3k.h5') # require diopy to load h5 to AnnData
# })

# AnnData to SeuratObject (oneMzell) --------
# rsync -R -r /data/user/syb/projects/scfetch/benchtime/runtimeOut/SCE2AD_zellkonverter/./*/*.h5ad ./
# for i in */*.h5ad;do new_name=$(echo ${i}|sed 's/SCE2AD_zellkonverter_//');mv ${i}  ${new_name};done
# all.h5ad = list.files(path = data.folder, pattern = "h5ad$", full.names = T)

# Error information
# scdior error: <simpleError in h5attr(h5mat, "datatype"): Attribute does not exist>
# SeuratDisk error: 'dims' must contain all (i,j) pairs

# * 1.SeuratDisk --------
cell.sub.vec = c('C120w', 'C100w', 'C80w', 'C60w', 'C50w', 'C30w', 'C20w', 'C10w', 'C8w', 'C6w', 'C5w', 'C3w', 'C2w', 'C1w')
dir.create(file.path(out.folder, 'AD2Seu_SeuratDisk'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'AD2Seu_SeuratDisk'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'AD2Seu_SeuratDisk', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'AD2Seu_SeuratDisk', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.h5ad = list.files(path = cs.data.folder, pattern = "h5ad$", full.names = T)
  for (h5ad in sub.h5ad){
    filename = basename(h5ad)
    out.h5seurat = gsub(pattern = ".h5ad$", replacement = ".h5seurat", x = h5ad)
    cat("runtime = system.time({",
        paste0("SeuratDisk::Convert('", h5ad, "', dest = 'h5seurat', assay = 'RNA', overwrite = TRUE)"),
        "# https://github.com/mojaveazure/seurat-disk/issues/109",
        paste0("f <- hdf5r::H5File$new('",out.h5seurat,"', 'r+')"),
        readLines(file.path(script.folder, 'AD2Seu_SeuratDisk_bug.R')),
        paste0("seu=SeuratDisk::LoadH5Seurat('",out.h5seurat,"', assays = 'RNA')"),
        "})",
        paste0("file.remove('",out.h5seurat,"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("AD2Seu_SeuratDisk_", gsub(pattern = ".h5ad$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("AD2Seu_SeuratDisk_",paste0(gsub(pattern = ".h5ad$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# system.time({
#   SeuratDisk::Convert("pbmc3k_final.h5ad", dest = "h5seurat", assay = "RNA", overwrite = TRUE)
#   fix bugs: https://github.com/mojaveazure/seurat-disk/issues/109
#   pbmc3k.seurat <- SeuratDisk::LoadH5Seurat("pbmc3k_final.h5seurat")
# })
# * 2.sceasy --------
dir.create(file.path(out.folder, 'AD2Seu_sceasy'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'AD2Seu_sceasy'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'AD2Seu_sceasy', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'AD2Seu_sceasy', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.h5ad = list.files(path = cs.data.folder, pattern = "h5ad$", full.names = T)
  for (h5ad in sub.h5ad){
    filename = basename(h5ad)
    cat("runtime = system.time({",
        "reticulate::use_condaenv('testsc', required = TRUE)",
        paste0("seu=sceasy::convertFormat('",h5ad,"', from='anndata', to='seurat', main_layer = 'scale.data', assay = 'RNA')"),
        "})",
        # second: don't run
        # paste0("saveRDS(seu, '", file.path(cs.out.folder, paste0("AD2Seu_sceasy_", gsub(pattern = ".h5ad$", replacement = "_seu.rds", x = filename))),"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("AD2Seu_sceasy_", gsub(pattern = ".h5ad$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("AD2Seu_sceasy_",paste0(gsub(pattern = ".h5ad$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# system.time({
#   pbmc3k.sceasy =
#     sceasy::convertFormat("pbmc3k_final.h5ad",
#                           from="anndata", to="seurat",
#                           main_layer = "scale.data", assay = "RNA")
# })

# * 3.scDIOR --------
dir.create(file.path(out.folder, 'AD2Seu_scdior'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'AD2Seu_scdior'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'AD2Seu_scdior', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'AD2Seu_scdior', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.h5ad = list.files(path = cs.data.folder, pattern = "h5ad$", full.names = T)
  for (h5ad in sub.h5ad){
    filename = basename(h5ad)
    cat("runtime = system.time({",
        "reticulate::use_condaenv('testsc', required = TRUE)",
        paste0("seu=dior::read_h5ad(file = '",h5ad,"', assay_name = 'RNA', target.object = 'seurat')"),
        "})",
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("AD2Seu_scdior_", gsub(pattern = ".h5ad$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("AD2Seu_scdior_",paste0(gsub(pattern = ".h5ad$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# system.time({
#   pbmc3k.scdior =
#     dior::read_h5ad(file = "pbmc3k_final.h5ad",
#                     assay_name="RNA", target.object = "seurat")
# })

# * 4.schard --------
dir.create(file.path(out.folder, 'AD2Seu_schard'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'AD2Seu_schard'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'AD2Seu_schard', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'AD2Seu_schard', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.h5ad = list.files(path = cs.data.folder, pattern = "h5ad$", full.names = T)
  for (h5ad in sub.h5ad){
    filename = basename(h5ad)
    cat("runtime = system.time({",
        paste0("seu=schard::h5ad2seurat(file = '",h5ad,"', use.raw = FALSE, assay = 'RNA')"),
        "})",
        # second: don't run
        # paste0("saveRDS(seu, '", file.path(cs.out.folder, paste0("AD2Seu_schard_", gsub(pattern = ".h5ad$", replacement = "_seu.rds", x = filename))),"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("AD2Seu_schard_", gsub(pattern = ".h5ad$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("AD2Seu_schard_",paste0(gsub(pattern = ".h5ad$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# system.time({
#   pbmc3k.schard =
#     schard::h5ad2seurat(file = "pbmc3k_final.h5ad",
#                         use.raw = TRUE, assay = "RNA")
# })

# * 5.SeuratDisk + scDIOR ---------
dir.create(file.path(out.folder, 'AD2Seu_SeuratDisk_add_scdior'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'AD2Seu_SeuratDisk_add_scdior'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'AD2Seu_SeuratDisk_add_scdior', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'AD2Seu_SeuratDisk_add_scdior', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.h5ad = list.files(path = cs.data.folder, pattern = "h5ad$", full.names = T)
  for (h5ad in sub.h5ad){
    filename = basename(h5ad)
    out.h5seurat = gsub(pattern = ".h5ad$", replacement = ".h5seurat", x = h5ad)
    cat("runtime = system.time({",
        "reticulate::use_condaenv('testsc', required = TRUE)",
        paste0("SeuratDisk::Convert('", h5ad, "', dest = 'h5seurat', assay = 'RNA', overwrite = TRUE)"),
        "# https://github.com/mojaveazure/seurat-disk/issues/109",
        paste0("f <- hdf5r::H5File$new('",out.h5seurat,"', 'r+')"),
        readLines(file.path(script.folder, 'AD2Seu_SeuratDisk_bug.R')),
        paste0("seu=SeuratDisk::LoadH5Seurat('",out.h5seurat,"', assays = 'RNA')"),
        paste0("seu.scdior=dior::read_h5ad(file = '",h5ad,"', assay_name = 'RNA', target.object = 'seurat')"),
        readLines(file.path(script.folder, 'AD2Seu_SeuratDisk_add_scDIOR.R')),
        "})",
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("AD2Seu_SeuratDisk_add_scdior_", gsub(pattern = ".h5ad$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("AD2Seu_SeuratDisk_add_scdior_",paste0(gsub(pattern = ".h5ad$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# AnnData to SingleCellExperiemnt (oneMh5ad) ---------
cell.sub.vec = c('C120w', 'C100w', 'C80w', 'C60w', 'C50w', 'C30w', 'C20w', 'C10w', 'C8w', 'C6w', 'C5w', 'C3w', 'C2w', 'C1w')
# * 1.scDIOR ---------
dir.create(file.path(out.folder, 'AD2SCE_scdior'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'AD2SCE_scdior'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'AD2SCE_scdior', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'AD2SCE_scdior', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.h5ad = list.files(path = cs.data.folder, pattern = "h5ad$", full.names = T)
  for (h5ad in sub.h5ad){
    filename = basename(h5ad)
    tmp.h5 = gsub(pattern = ".h5ad$", replacement = "_tmp.h5", x = h5ad)
    cat("runtime = system.time({",
        "reticulate::use_condaenv('testsc', required = TRUE)",
        "anndata <- reticulate::import('anndata')",
        paste0("adata = anndata$read_h5ad('",h5ad,"')"),
        "diopy <- reticulate::import('diopy')",
        paste0("diopy$output$write_h5(adata = adata$raw$to_adata(), file='",tmp.h5,"', assay_name = 'RNA', save_X = TRUE)"),
        paste0("sce = dior::read_h5(file = '",tmp.h5,"', target.object = 'singlecellexperiment')"),
        "})",
        paste0("file.remove('",tmp.h5,"')"),
        # second: don't run
        # paste0("saveRDS(sce, '", file.path(cs.out.folder, paste0("AD2SCE_scdior_", gsub(pattern = ".h5ad$", replacement = "_sce.rds", x = filename))),"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("AD2SCE_scdior_", gsub(pattern = ".h5ad$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("AD2SCE_scdior_",paste0(gsub(pattern = ".h5ad$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# anndata <- reticulate::import("anndata")
# adata = anndata$read_h5ad(anndata.file)
# diopy = reticulate::import("diopy")
# h5.file = gsub(pattern = ".h5ad$", replacement = "_tmp.h5", x = anndata.file)
# diopy$output$write_h5(adata = adata$raw$to_adata(), file=h5.file, assay_name=assay, save_X = TRUE)
# dior::read_h5(file = h5.file, target.object = "singlecellexperiment")

# * 2.zellkonverter ---------
dir.create(file.path(out.folder, 'AD2SCE_zellkonverter'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'AD2SCE_zellkonverter'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'AD2SCE_zellkonverter', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'AD2SCE_zellkonverter', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.h5ad = list.files(path = cs.data.folder, pattern = "h5ad$", full.names = T)
  for (h5ad in sub.h5ad){
    filename = basename(h5ad)
    cat("runtime = system.time({",
        "reticulate::use_condaenv('testsc', required = TRUE)",
        "anndata <- reticulate::import('anndata')",
        paste0("adata = anndata$read_h5ad('",h5ad,"')"),
        "sce = zellkonverter::AnnData2SCE(adata, X_name = 'scale.data', raw = TRUE)",
        "})",
        # second: don't run
        # paste0("saveRDS(sce, '", file.path(cs.out.folder, paste0("AD2SCE_zellkonverter_", gsub(pattern = ".h5ad$", replacement = "_sce.rds", x = filename))),"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("AD2SCE_zellkonverter_", gsub(pattern = ".h5ad$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("AD2SCE_zellkonverter_",paste0(gsub(pattern = ".h5ad$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# anndata <- reticulate::import("anndata")
# adata <- anndata$read_h5ad(anndata.file)
# zellkonverter::AnnData2SCE(adata, X_name = 'scale.data', raw = TRUE)

# * 3.schard ---------
dir.create(file.path(out.folder, 'AD2SCE_schard'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'AD2SCE_schard'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'AD2SCE_schard', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'AD2SCE_schard', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.h5ad = list.files(path = cs.data.folder, pattern = "h5ad$", full.names = T)
  for (h5ad in sub.h5ad){
    filename = basename(h5ad)
    cat("runtime = system.time({",
        paste0("sce = schard::h5ad2sce('",h5ad,"', use.raw = TRUE)"),
        "})",
        # second: don't run
        # paste0("saveRDS(sce, '", file.path(cs.out.folder, paste0("AD2SCE_schard_", gsub(pattern = ".h5ad$", replacement = "_sce.rds", x = filename))),"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("AD2SCE_schard_", gsub(pattern = ".h5ad$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("AD2SCE_schard_",paste0(gsub(pattern = ".h5ad$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# schard::h5ad2sce(anndata.file, use.raw = TRUE)

# SingleCellExperiemnt to AnnData (oneMsce) ----------
# rsync -R -r /data/user/syb/projects/scfetch/benchtime/runtimeOut/AD2SCE_zellkonverter/./*/*sce.rds ./
# for i in */*.rds;do new_name=$(echo ${i}|sed 's/AD2SCE_zellkonverter_//');mv ${i}  ${new_name};done
cell.sub.vec = c('C120w', 'C100w', 'C80w', 'C60w', 'C50w', 'C30w', 'C20w', 'C10w', 'C8w', 'C6w', 'C5w', 'C3w', 'C2w', 'C1w')
# * 1.sceasy -------
dir.create(file.path(out.folder, 'SCE2AD_sceasy'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'SCE2AD_sceasy'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'SCE2AD_sceasy', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'SCE2AD_sceasy', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.rds = list.files(path = cs.data.folder, pattern = "rds$", full.names = T)
  for (rds in sub.rds){
    filename = basename(rds)
    out.h5ad = file.path(cs.out.folder, paste0("SCE2AD_sceasy_", gsub(pattern = ".rds$", replacement = ".h5ad", x = filename)))
    cat(paste0("sce=readRDS('",rds,"')"),
        "runtime = system.time({",
        "reticulate::use_condaenv('testsc', required = TRUE)",
        paste0("sceasy::convertFormat(sce, from='sce', to='anndata', drop_single_values = F, outFile = '",out.h5ad,"', main_layer = 'scale.data')"),
        "})",
        # second: remove output h5ad file
        paste0("file.remove('",out.h5ad,"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("SCE2AD_sceasy_", gsub(pattern = ".rds$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("SCE2AD_sceasy_",paste0(gsub(pattern = ".rds$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# sceasy::convertFormat(sce, from="sce", to="anndata", drop_single_values = FALSE,
#                       outFile="sce.h5ad", main_layer = "scale.data")

# * 2.zellkonverter ---------
dir.create(file.path(out.folder, 'SCE2AD_zellkonverter'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'SCE2AD_zellkonverter'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'SCE2AD_zellkonverter', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'SCE2AD_zellkonverter', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.rds = list.files(path = cs.data.folder, pattern = "rds$", full.names = T)
  for (rds in sub.rds){
    filename = basename(rds)
    out.h5ad = file.path(cs.out.folder, paste0("SCE2AD_zellkonverter_", gsub(pattern = ".rds$", replacement = ".h5ad", x = filename)))
    cat(paste0("sce=readRDS('",rds,"')"),
        "runtime = system.time({",
        "reticulate::use_condaenv('testsc', required = TRUE)",
        "anndata <- reticulate::import('anndata')",
        "adata = zellkonverter::SCE2AnnData(sce, X_name = 'scale.data')",
        paste0("adata$write_h5ad('",out.h5ad,"')"),
        "})",
        # second: remove output h5ad file
        paste0("file.remove('",out.h5ad,"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("SCE2AD_zellkonverter_", gsub(pattern = ".rds$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("SCE2AD_zellkonverter_",paste0(gsub(pattern = ".rds$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# anndata <- reticulate::import("anndata")
# adata <- zellkonverter::SCE2AnnData(sce, X_name = 'scale.data')
# adata$write_h5ad("sce.h5ad")

# * 3.scDIOR ---------
dir.create(file.path(out.folder, 'SCE2AD_scdior'), showWarnings = FALSE)
dir.create(file.path(script.folder, 'SCE2AD_scdior'), showWarnings = FALSE)
for (cs in cell.sub.vec){
  # prepare folder
  cs.data.folder = file.path(data.folder, cs)
  cs.out.folder = file.path(out.folder, 'SCE2AD_scdior', cs)
  dir.create(cs.out.folder, showWarnings = FALSE)
  cs.script.folder =  file.path(script.folder, 'SCE2AD_scdior', cs)
  dir.create(cs.script.folder, showWarnings = FALSE)
  # generate scripts
  sub.rds = list.files(path = cs.data.folder, pattern = "rds$", full.names = T)
  for (rds in sub.rds){
    filename = basename(rds)
    out.h5 = file.path(cs.out.folder, paste0("SCE2AD_scdior_", gsub(pattern = ".rds$", replacement = ".h5", x = filename)))
    cat(paste0("sce=readRDS('",rds,"')"),
        "runtime = system.time({",
        paste0("dior::write_h5(data = sce, object.type='singlecellexperiment', file = '",out.h5,"')"),
        "})",
        "# adata = diopy.input.read_h5(file = 'pbmc3k.h5') # require diopy to load h5 to AnnData",
        # second: remove output h5 file
        paste0("file.remove('",out.h5,"')"),
        paste0("saveRDS(runtime, '", file.path(cs.out.folder, paste0("SCE2AD_scdior_", gsub(pattern = ".rds$", replacement = "_runtime.rds", x = filename))),"')"),
        file=file.path(cs.script.folder, paste0("SCE2AD_scdior_",paste0(gsub(pattern = ".rds$", replacement = ".R", x = filename)))),sep="\n")
  }
}
# example
# dior::write_h5(data = sce, object.type = "singlecellexperiment",
#                file = "sce.h5")

# run scripts -----------
# Linux, SCE2AD_sceasy as an example
# cd SCE2AD_sceasy
# for i in C*/SCE2AD*.R
# do
#     ~/miniconda3/envs/testsc/bin/Rscript ${i} >${i}.log 2>&1
#     sleep 1m
# done

