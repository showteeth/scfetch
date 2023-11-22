# scfetch 0.5.0

## New features
* Added `ExtractCELLxGENEMeta` to extract metadata of CELLxGENE datasets.
* Added `ParseCELLxGENE` to download objects from CELLxGENE.
* Added `ExtractHCAMeta` to extract metadata of Human Cell Atlas projects.
* Added `ShowHCAProjects` to show all available Human Cell Atlas projects.
* Added `ShowCELLxGENEDatasets` to show all available CELLxGENE datasets.
* Added `StatDBAttribute`to stat database attributes.
* Added `ParseHCA` to download objects from Human Cell Atlas.

-------------------

## Minor changes
* Added filters to `ExtractCELLxGENEMeta`.
* Added `data.type` in `ParseGEO` to return count matrix for bulk RNA-seq.
* Added filters to `ExtractHCAMeta`.
* Added `catalog` in `ShowHCAProjects`.
* Fixed bug in `ParseCELLxGENE` and returned failed dataframe when downloading error.
* Fixed bug in `ImportSeurat`.
* Fixed bug in `ExportSeurat`.
* Resolved installation.
* Fixed bug in `ExtractZenodoMeta` (API changed).
* Fixed bug in `ShowCBDatasets` (added `--no-check-certificate` when downloading json files).
* Fixed bug in `ShowCELLxGENEDatasets`.

-------------------

# scfetch 0.4.0

## New features
* Added `ExportSeurat` to convert SeuratObject to other scRNA-seq formats.
* Added `ImportSeurat` to convert other scRNA-seq formats to SeuratObject.
* Added `SCEaAnnData` to perform data format conversion between SingleCellExperiment and AnnData.
* Added `SCELoom` to perform data format conversion between SingleCellExperiment and loom.
* Added `ExtractPanglaoDBComposition` to extract cell type composition of PanglaoDB datasets.

-------------------

## Minor changes
* Added `conda.path` in `ExportSeurat` to specify conda enviroment used.
* Fixed bug in `SCELoom` to keep gene names and cell ID.
* `Bam2Fastq` supported normal bam files (non-10x bam files).
* Changed function names: `ShowPanglaoDBMeta` to `ExtractPanglaoDBMeta`, `PrepareZenodo` to `ExtractZenodoMeta`.
* Supported downloadding normal (non-10x) bam file in `DownloadBam`.
* `ExtractZenodoMeta` supported a vector of Zenodo dois.
* Added `fastq.type` in `SplitSRA` to deal with fastq files from 10x, other scRNA-seq protocols and bulk RNA-seq.
* Added check for `split.cmd.paras` in `SplitSRA`.
* Added `local.data` in `ExtractPanglaoDBMeta` and `ExtractPanglaoDBComposition` to use cached sample metadata and composition.

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
