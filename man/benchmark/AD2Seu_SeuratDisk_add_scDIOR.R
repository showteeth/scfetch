# add additional assays
all.assays = Seurat::Assays(seu.scdior)
unused.assays = setdiff(all.assays, 'RNA')
if(length(unused.assays) > 0){
  for (ay in unused.assays){
    # https://github.com/JiekaiLab/dior/blob/2b1ea47b6661c8a10d9455f3baeeccb8f12be2f0/R/seuratIO.R#L65
    # https://github.com/satijalab/seurat-object/blob/58bf437fe058dd78913d9ef7b48008a3e24a306a/R/assay.R#L157
    assay.data <- Seurat::GetAssayData(object =  seu.scdior[[ay]], slot = 'counts')
    seu[[ay]] <- Seurat::CreateAssayObject(counts = assay.data )
  }
}
# add graphs
seu@graphs = seu.scdior@graphs
