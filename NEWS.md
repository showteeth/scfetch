# GEfetch2R 0.7.1

## New features
* Added `Seu2AD`/`SCE2AD`/`AD2Seu`/`AD2SCE` for benchmark (data IO).

-------------------

## Minor changes
* Added `timeout` in `ShowCBDatasets`.
* Added `Gunzip` to deal with the removal of `gunzip` in package `GEOquery`.

-------------------

# GEfetch2R 0.7.0

## New features
* Added `DownloadFastq` to download FASTQ files directly from ENA.
* Updated `DownloadSRA` to download sra files directly from ENA.
* Updated `DownloadBam` to download bam files directly from ENA.

-------------------

## Minor changes
* Updated `DownloadSRA`, `DownloadFastq`, `DownloadBam` to keep the same file structure.
* Updated `DownloadFastq` to format the fastqs to 10x format.

-------------------


# GEfetch2R 0.6.1

* Fix bugs in `LoadRDS2Seurat`.
* Give a message when there is no file with `file.ext` in `ParseHCA`, `ParseHCA`, `ExtractZenodoMeta`, `ParseZenodo`.
* Return SeuratObject in `ParseHCA`.
* Fix bug in `Load2Seurat`.
* Support `pigz` in `Bam2Fastq`.
* Return full dataframe in `ExtractRun`.
* Fix bug in `ParseCELLxGENE`, `ParseHCA`.
* Add `use.cores` in `ParseCELLxGENE`, `ParseHCA`, `ParseZenodo`.
* Fix bug in `ExtractGEOMeta`.

-------------------

# GEfetch2R 0.6.0

## New features
* Added `RunCellRanger` to run CellRanger on downloaded FASTQ files.
* Added `RunSTAR` to run STAR on downloaded FASTQ files.
* Added `Fastq2R` to pipe FASTQ files to `SeuratObject` and `DESeqDataSet`.
* Updated `StatDBAttribute` to support `cellxgene.census`.
* Updated `ParseCELLxGENE` to support `cellxgene.census` and return SeuratObject.
* Updated `ParseCBDatasets` to support subset.

-------------------

## Minor changes
* Return `DESeqDataSet` when `data.type` is bulk in `ParseGEO`.
* Support summarising multiple attributes in `StatDBAttribute`.
* Fix bugs in `mergeExperiments`.
* Fix bugs in `ParseHCA` (no file with extension specified by `file.ext`).
* Optimize `ParseHCA` to return metadata of downloaded files.

-------------------

# GEfetch2R 0.5.1

* Create output directory automatically.
* CELLxGENE new [API](https://api.cellxgene.cziscience.com/curation/ui/).

-------------------

# GEfetch2R 0.5.0

## New features
* Added `ExtractCELLxGENEMeta` to extract metadata of CELLxGENE datasets.
* Added `ParseCELLxGENE` to download objects from CELLxGENE.
* Added `ExtractHCAMeta` to extract metadata of Human Cell Atlas projects.
* Added `ShowHCAProjects` to show all available Human Cell Atlas projects.
* Added `ShowCELLxGENEDatasets` to show all available CELLxGENE datasets.
* Added `StatDBAttribute`to stat database attributes.
* Added `ParseHCA` to download objects from Human Cell Atlas.
* docker image added. 

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
* Supported GEO of 10x (separate files).

-------------------

# GEfetch2R 0.4.0

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

# GEfetch2R 0.3.0

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

# GEfetch2R 0.2.0

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

# GEfetch2R 0.1.0

* Provided APIs for Zenodo, PanglaoDB and UCSC Cell Browser.
