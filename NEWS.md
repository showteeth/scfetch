# scfetch 0.4.0

## New features
* Added `ExportSeurat` to convert SeuratObject to other scRNA-seq formats.
* Added `ImportSeurat` to convert other scRNA-seq formats to SeuratObject.
* Added `SCEaAnnData` to perform data format conversion between SingleCellExperiment and AnnData.
* Added `SCELoom` to perform data format conversion between SingleCellExperiment and loom.

-------------------

## Minor changes
* Added `conda.path` in `ExportSeurat` to specify conda enviroment used.
* Fixed bug in `SCELoom` to keep gene names and cell ID.

-------------------

# scfetch 0.3.0

## New features
* Added `ExtractRun` to extract run number from GEO.
* Added `DownloadSRA` to download SRA according to run number.
* Added `SplitSRA` to split SRA to fastqs and format the fastqs to 10x format.
* Added `DownloadBam` to download 10x bam files from GEO.
* Added `Bam2Fastq` to convert bam to 10x formatted fastqs.

-------------------

## Minor changes
* Fixed bug in `DownloadSRA` when downloading bam files.

-------------------

# scfetch 0.2.0

## New features
* Provided APIs for GEO.
* Added `Read10XOnline` to load cellranger output.
* Added `ExtractCBComposition` to extract cell type composition of UCSC Cell Browser datasets.
* `ParseGEO` support creating Seurat object.
* Added `ExtractGEOMeta` to support extract metadata from GEO.

-------------------

## Minor changes
* Supported lazy mode in `ShowCBDatasets` to save time for multi-runs.
* Optimized the paras for `ExtractCBDatasets` and `ParseCBDatasets`.
* Added `cell.num` filter to `ParseCBDatasets`, `ExtractCBDatasets` and `ShowPanglaoDBMeta`.
* Extracted matrix, barcode and feature info from dataset.json (comprehensive) instead of desc.json.
* Fixed bugs in `ParseCBDatasets` when dealing with cellranger output.
* Added `timeout` in `ParseCBDatasets` to avoid possible timeout error.
* Changed `time.out` in GEO related functions to `timeout`.
* Simplified the output of `ParseGEO`. 
* Fixed bug in `ExtractGEOMeta`.

-------------------

# scfetch 0.1.0

* Provided APIs for Zenodo, PanglaoDB and UCSC Cell Browser.
