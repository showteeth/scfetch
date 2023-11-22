
# scfetch - Access and Format Single-cell RNA-seq Datasets from Public Resources 

<img src = "man/figures/scfetch.png" align = "right" width = "200"/>

## Introduction

`scfetch` is designed to accelerate users download and prepare single-cell datasets from public resources. It can be used to:

* **Download fastq files** from `GEO/SRA`, **foramt fastq files** to standard style that can be identified by 10x softwares (e.g. CellRanger).
* **Download bam files** from `GEO/SRA`, support **downloading original 10x generated bam files (with custom tags) and normal bam files**, and **convert bam files to fastq files**.
* Download scRNA-seq **matrix** and **annotation (e.g. cell type)** information from `GEO`, `PanglanDB` and `UCSC Cell Browser`, **load the downnloaded matrix to `Seurat`**.
* Download processed objects from `Zeenodo`, `CELLxGENE` and `Human Cell Atlas`.
* **Formats conversion between widely used single cell objects** (`SeuratObject`, `AnnData`, `SingleCellExperiment`, `CellDataSet/cell_data_set` and `loom`).

<hr />

## Framework

<div align="center">
<img src="man/figures/scfetch_workflow.png"  title="scfetch_framework"  alt="scfetch_framework" />
</div>

<hr />

## Installation

`scfetch` is an R package distributed as part of the [CRAN](https://cran.r-project.org/web/packages/scfetch/index.html). To install the package, start R and enter:
```R
# install via CRAN
install.packages("scfetch")

# you can also install the development version from GitHub
# install.packages("devtools")
devtools::install_github("showteeth/scfetch")
```

There are some conditionally used packages:

```R
# install.packages("devtools") #In case you have not installed it.
devtools::install_github("alexvpickering/GEOfastq") # download fastq
devtools::install_github("cellgeni/sceasy") # format conversion
devtools::install_github("mojaveazure/seurat-disk") # format conversion
devtools::install_github("satijalab/seurat-wrappers") # format conversion
```

**For issues about installation, please refer `INSTALL.md`.**

For data structures conversion, `scfetch` requires several python packages, you can install with:

``` bash
# install python packages
conda install -c bioconda loompy anndata
# or
pip install anndata loompy
```

<hr />

## Usage
### Vignette
Detailed usage is available in [website](https://showteeth.github.io/scfetch/).

<hr />

### Function list

<table>
<thead>
  <tr>
    <th>Type</th>
    <th>Function</th>
    <th>Usage</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="3">Download and format fastq</td>
    <td>ExtractRun</td>
    <td>Extract runs with GEO accession number or GSM number</td>
  </tr>
  <tr>
    <td>DownloadSRA</td>
    <td>Download sra files</td>
  </tr>
  <tr>
    <td>SplitSRA</td>
    <td>Split sra files to fastq files and format to 10x standard style</td>
  </tr>
  <tr>
    <td rowspan="2">Download and convert bam</td>
    <td>DownloadBam</td>
    <td>Download bam (support 10x original bam)</td>
  </tr>
  <tr>
    <td>Bam2Fastq</td>
    <td>Convert bam files to fastq files</td>
  </tr>
  <tr>
    <td rowspan="9">Download matrix and load to Seurat </td>
    <td>ExtractGEOMeta</td>
    <td>Extract sample metadata from GEO</td>
  </tr>
  <tr>
    <td>ParseGEO</td>
    <td>Download matrix from GEO and load to Seurat</td>
  </tr>
  <tr>
    <td>ExtractPanglaoDBMeta</td>
    <td>Extract sample metadata from PandlaoDB</td>
  </tr>
  <tr>
    <td>ExtractPanglaoDBComposition</td>
    <td>Extract cell type composition of PanglaoDB datasets</td>
  </tr>
  <tr>
    <td>ParsePanglaoDB</td>
    <td>Download matrix from PandlaoDB and load to Seurat</td>
  </tr>
  <tr>
    <td>ShowCBDatasets</td>
    <td>Show all available datasets in UCSC Cell Browser</td>
  </tr>
  <tr>
    <td>ExtractCBDatasets</td>
    <td>Extract UCSC Cell Browser datasets with attributes</td>
  </tr>
  <tr>
    <td>ExtractCBComposition</td>
    <td>Extract cell type composition of UCSC Cell Browser datasets</td>
  </tr>
  <tr>
    <td>ParseCBDatasets</td>
    <td>Download UCSC Cell Browser datasets and load to Seurat</td>
  </tr>
  <tr>
    <td rowspan="8">Download objects</td>
    <td>ExtractZenodoMeta</td>
    <td>Extract sample metadata from Zenodo with DOIs</td>
  </tr>
  <tr>
    <td>ParseZenodo</td>
    <td>Download rds/rdata/h5ad/loom from Zenodo with DOIs</td>
  </tr>
  <tr>
    <td>ShowCELLxGENEDatasets</td>
    <td>Show all available datasets in CELLxGENE</td>
  </tr>
  <tr>
    <td>ExtractCELLxGENEMeta</td>
    <td>Extract metadata of CELLxGENE datasets with attributes</td>
  </tr>
  <tr>
    <td>ParseCELLxGENE</td>
    <td>Download rds/h5ad from CELLxGENE</td>
  </tr>
  <tr>
    <td>ShowHCAProjects</td>
    <td>Show all available projects in Human Cell Atlas</td>
  </tr>
  <tr>
    <td>ExtractHCAMeta</td>
    <td>Extract metadata of Human Cell Atlas projects with attributes</td>
  </tr>
  <tr>
    <td>ParseHCA</td>
    <td>Download rds/rdata/h5/h5ad/loom from Human Cell Atlas</td>
  </tr>
  <tr>
    <td rowspan="4">Convert between different single-cell objects</td>
    <td>ExportSeurat</td>
    <td>Convert SeuratObject to AnnData, SingleCellExperiment, CellDataSet/cell_data_set and loom</td>
  </tr>
  <tr>
    <td>ImportSeurat</td>
    <td>Convert AnnData, SingleCellExperiment, CellDataSet/cell_data_set and loom to SeuratObject</td>
  </tr>
  <tr>
    <td>SCEAnnData</td>
    <td>Convert between SingleCellExperiment and AnnData</td>
  </tr>
  <tr>
    <td>SCELoom</td>
    <td>Convert between SingleCellExperiment and loom</td>
  </tr>
  <tr>
    <td rowspan="1">Summarize datasets based on attributes</td>
    <td>StatDBAttribute</td>
    <td>Summarize datasets in PandlaoDB, UCSC Cell Browser and CELLxGENE based on attributes</td>
  </tr>
</tbody>
</table>

<hr />

### Downloas fastq and bam

Since the downloading process is time-consuming, we provide the commands used to illustrate the usage.

#### Downloas fastq

##### Prepare run number

For fastq files stored in SRA, `scfetch` can extract sample information and run number with GEO accession number or users can also provide a dataframe contains the run number of interested samples.

Extract all samples under `GSE130636` and the platform is `GPL20301` (use `platform = NULL` for all platforms):
```r
GSE130636.runs <- ExtractRun(acce = "GSE130636", platform = "GPL20301")
```

<hr />

##### Download sra

With the dataframe contains gsm and run number, `scfetch` will download sra files using `prefetch`. The returned result is a dataframe contains failed runs. If not `NULL`, users can re-run `DownloadSRA` by setting `gsm.df` to the returned result.

```r
# a small test
GSE130636.runs <- GSE130636.runs[GSE130636.runs$run %in% c("SRR9004346", "SRR9004351"), ]
# download, you may need to set prefetch.path
out.folder <- tempdir()
GSE130636.down <- DownloadSRA(
  gsm.df = GSE130636.runs,
  out.folder = out.folder
)
# GSE130636.down is null or dataframe contains failed runs
```

The `out.folder` structure will be: `gsm_number/run_number`.

<hr />

##### Split fastq

After obtaining the sra files, `scfetch` provides function `SplitSRA` to split sra files to fastq files using `parallel-fastq-dump` (parallel, fastest and gzip output), `fasterq-dump` (parallel, fast but unzipped output) and `fastq-dump` (slowest and gzip output).

For fastqs generated with 10x Genomics, `SplitSRA` can identify read1, read2 and index files and format the read1 and read2 to 10x required format (`sample1_S1_L001_R1_001.fastq.gz` and `sample1_S1_L001_R2_001.fastq.gz`). In detail, the file with read length 26 or 28 is considered as read1, the files with read length 8 or 10 are considered as index files and the remain file is considered as read2. The read length rules is from [Sequencing Requirements for Single Cell 3'](https://www.10xgenomics.com/cn/support/single-cell-gene-expression/documentation/steps/sequencing/sequencing-requirements-for-single-cell-3) and [Sequencing Requirements for Single Cell V(D)J](https://www.10xgenomics.com/cn/support/single-cell-immune-profiling/documentation/steps/sequencing/sequencing-requirements-for-single-cell-v-d-j).

The returned result is a vector of failed sra files. If not `NULL`, users can re-run `SplitSRA` by setting `sra.path` to the returned result.

```r
# parallel-fastq-dump requires sratools.path
# you may need to set split.cmd.path and sratools.path
sra.folder <- tempdir()
GSE130636.split <- SplitSRA(
  sra.folder = sra.folder,
  fastq.type = "10x", split.cmd.threads = 4
)
```

<hr />

#### Download bam

##### Prepare run number

`scfetch` can extract sample information and run number with GEO accession number or users can also provide a dataframe contains the run number of interested samples.

```r
GSE138266.runs <- ExtractRun(acce = "GSE138266", platform = "GPL18573")
```

<hr />

##### Download bam

With the dataframe contains gsm and run number, `scfetch` provides `DownloadBam` to download bam files using `prefetch`. It suooorts 10x generated bam files and normal bam files.

* 10x generated bam: While bam files generated from 10x softwares (e.g. CellRanger) contain custom tags which are not kept when using default parameters of `prefetch`, `scfetch` adds `--type TenX` to make sure the downloaded bam files contain these tags. 
* normal bam: For normal bam files, `DownloadBam` will download sra files first and then convert sra files to bam files with `sam-dump`. After testing the efficiency of `prefetch` + `sam-dump` and `sam-dump`, the former is much faster than the latter (52G sra and 72G bam files):
```bash
# # use prefetch to download sra file
# prefetch -X 60G SRR1976036
# # real	117m26.334s
# # user	16m42.062s
# # sys	3m28.295s

# # use sam-dump to convert sra to bam
# time (sam-dump SRR1976036.sra | samtools view -bS - -o SRR1976036.bam)
# # real	536m2.721s
# # user	749m41.421s
# # sys	20m49.069s


# use sam-dump to download bam directly
# time (sam-dump SRR1976036 | samtools view -bS - -o SRR1976036.bam)
# # more than 36hrs only get ~3G bam files, too slow
```

The returned result is a dataframe containing failed runs (either failed to download sra files or failed to convert to bam files for normal bam; failed to download bam files for 10x generated bam). If not `NULL`, users can re-run `DownloadBam` by setting `gsm.df` to the returned result. The following is an example to download 10x generated bam file:

```r
# a small test
GSE138266.runs <- GSE138266.runs[GSE138266.runs$run %in% c("SRR10211566"), ]
# download, you may need to set prefetch.path
out.folder <- tempdir()
GSE138266.down <- DownloadBam(
  gsm.df = GSE138266.runs,
  out.folder = out.folder
)
# GSE138266.down is null or dataframe contains failed runs
```

The `out.folder` structure will be: `gsm_number/run_number`.

<hr />

##### Convert bam to fastq

With downloaded bam files, `scfetch` provides function `Bam2Fastq` to convert bam files to fastq files. For bam files generated from 10x softwares, `Bam2Fastq` utilizes `bamtofastq` tool developed by 10x Genomics.

The returned result is a vector of bam files failed to convert to fastq files. If not `NULL`, users can re-run `Bam2Fastq` by setting `bam.path` to the returned result.

```r
bam.folder <- tempdir()
# you may need to set bamtofastq.path and bamtofastq.paras
GSE138266.convert <- Bam2Fastq(
  bam.folder = bam.folder
)
```

<hr />

### Download count matrix

`scfetch` provides functions for users to download **count matrices** and **annotations** (e.g. cell type annotation and composition) from GEO and some single-cell databases (e.g. [PanglaoDB](https://panglaodb.se/index.html) and [UCSC Cell Browser](https://cells.ucsc.edu/?#)). `scfetch` also supports loading the downloaded data to `Seurat`.

Until now, the public resources supported and the returned results:

| Resources         | URL                               | Download Type | Returned results                              |
|-------------------|-----------------------------------|---------------|-----------------------------------------------|
| GEO               | https://www.ncbi.nlm.nih.gov/geo/ | count matrix  | SeuratObject or count matrix for bulk RNA-seq |
| PanglaoDB         | https://panglaodb.se/index.html   | count matrix  | SeuratObject                                  |
| UCSC Cell Browser | https://cells.ucsc.edu/           | count matrix  | SeuratObject                                  |

<hr />

#### GEO

[GEO is an international public repository that archives and freely distributes microarray, next-generation sequencing, and other forms of high-throughput functional genomics data submitted by the research community.](https://www.ncbi.nlm.nih.gov/geo/info/overview.html) It provides a very convenient way for users to explore and select interested scRNA-seq datasets. 

##### Extract metadata

`scfetch` provides `ExtractGEOMeta` to extract sample metadata, including sample title, source name/tissue, description, cell type, treatment, paper title, paper abstract, organism, protocol, data processing methods, et al.

```r
# extract metadata of specified platform
GSE200257.meta <- ExtractGEOMeta(acce = "GSE200257", platform = "GPL24676")
# set VROOM_CONNECTION_SIZE to avoid error: Error: The size of the connection buffer (786432) was not large enough
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 60)
# extract metadata of all platforms
GSE94820.meta <- ExtractGEOMeta(acce = "GSE94820", platform = NULL)
```

<hr />

##### Download matrix and load to Seurat

After manually check the extracted metadata, users can **download count matrix** and **load the count matrix** to Seurat with `ParseGEO`. 

For count matrix, `ParseGEO` supports downloading the matrix from **supplementary files** and extracting from `ExpressionSet`, users can control the source by specifying `down.supp` or detecting automatically (`ParseGEO` will extract the count matrix from `ExpressionSet` first, if the count matrix is NULL or contains non-integer values, `ParseGEO` will download supplementary files). While the supplementary files have two main types: single count matrix file containing all cells and CellRanger-style outputs (barcode, matrix, feature/gene), users are required to choose the type of supplementary files with `supp.type`.

With the count matrix, `ParseGEO` will load the matrix to Seurat automatically. If multiple samples available, users can choose to merge the SeuratObject with `merge`.

```r
# for cellranger output
out.folder <- tempdir()
GSE200257.seu <- ParseGEO(
  acce = "GSE200257", platform = NULL, supp.idx = 1, down.supp = TRUE, supp.type = "10x",
  out.folder = out.folder
)
# for count matrix, no need to specify out.folder, download count matrix to tmp folder
GSE94820.seu <- ParseGEO(acce = "GSE94820", platform = NULL, supp.idx = 1, down.supp = TRUE, supp.type = "count")
```

**For bulk RNA-seq**, set `data.type = "bulk"` in `ParseGEO`, this will return count matrix.

<hr />

#### PanglaoDB

[PanglaoDB](https://panglaodb.se/index.html) is a database which contains scRNA-seq datasets from mouse and human. Up to now, it contains **5,586,348 cells** from **1368 datasets (1063 from Mus musculus and 305 from	Homo sapiens)**. It has well organized metadata for every dataset, including tissue, protocol, species, number of cells and cell type annotation (computationally identified). Daniel Osorio has developed [rPanglaoDB](https://github.com/dosorio/rPanglaoDB/) to access [PanglaoDB](https://panglaodb.se/index.html) data, the functions of `scfetch` here are based on [rPanglaoDB](https://github.com/dosorio/rPanglaoDB/).

Since [PanglaoDB](https://panglaodb.se/about.html) is no longer maintained, `scfetch` has cached all metadata and cell type composition and use these cached data by default to accelerate, users can access the cached data with `PanglaoDBMeta` (all metadata) and `PanglaoDBComposition` (all cell type composition).

##### Summarise attributes

`scfetch` provides `StatDBAttribute` to summary attributes of [PanglaoDB](https://panglaodb.se/index.html):

```r
# use cached metadata
StatDBAttribute(df = PanglaoDBMeta, filter = c("species", "protocol"), database = "PanglaoDB")
```

<hr />

##### Extract metadata

`scfetch` provides `ExtractPanglaoDBMeta` to select interested datasets with specified **species**, **protocol**, **tissue** and **cell number** (The available values of these attributes can be obtained with `StatDBAttribute`). User can also choose to whether to add cell type annotation to every dataset with `show.cell.type`.

`scfetch` uses cached metadata and cell type composition by default, users can change this by setting `local.data = FALSE`.

```r
hsa.meta <- ExtractPanglaoDBMeta(
  species = "Homo sapiens", protocol = c("Smart-seq2", "10x chromium"),
  show.cell.type = TRUE, cell.num = c(1000, 2000)
)
```

<hr />

##### Extract cell type composition

`scfetch` provides `ExtractPanglaoDBComposition` to extract cell type annotation and composition (use cached data by default to accelerate, users can change this by setting `local.data = FALSE`).

```r
hsa.composition <- ExtractPanglaoDBComposition(species = "Homo sapiens", protocol = c("Smart-seq2", "10x chromium"))
```

<hr />

##### Download matrix and load to Seurat

After manually check the extracted metadata, `scfetch` provides `ParsePanglaoDB` to **download count matrix** and **load the count matrix** to Seurat. With available cell type annotation, uses can filter datasets without specified cell type with `cell.type`. Users can also include/exclude cells expressing specified genes with `include.gene`/`exclude.gene`. 

With the count matrix, `ParsePanglaoDB` will load the matrix to Seurat automatically. If multiple datasets available, users can choose to merge the SeuratObject with `merge`.

```r
# small test
hsa.seu <- ParsePanglaoDB(hsa.meta[1:3, ], merge = TRUE)
```

<hr />

#### UCSC Cell Browser

The [UCSC Cell Browser](https://cells.ucsc.edu/?#) is a web-based tool that allows scientists to interactively visualize scRNA-seq datasets. It contains **1040 single cell datasets** from **17 different species**. And, it is **organized with the hierarchical structure**, which can help users quickly locate the datasets they are interested in.

##### Show available datasets

`scfetch` provides `ShowCBDatasets` to show all available datasets. Due to the large number of datasets, `ShowCBDatasets` enables users to perform *lazy load* of dataset json files instead of downloading the json files online (time-consuming!!!). This *lazy load* requires users to provide `json.folder` to save json files and set `lazy = TRUE` (for the first time of run, `ShowCBDatasets` will download current json files to `json.folder`, for next time of run, with `lazy = TRUE`, `ShowCBDatasets` will load the downloaded json files from `json.folder`.). And, `ShowCBDatasets` supports updating the local datasets with `update = TRUE`.

```r
json.folder <- tempdir()
# first time run, the json files are stored under json.folder
# ucsc.cb.samples = ShowCBDatasets(lazy = TRUE, json.folder = json.folder, update = TRUE)

# second time run, load the downloaded json files
ucsc.cb.samples <- ShowCBDatasets(lazy = TRUE, json.folder = json.folder, update = FALSE)

# always read online
# ucsc.cb.samples = ShowCBDatasets(lazy = FALSE)
```

The number of datasets and all available species:

```r
# the number of datasets
nrow(ucsc.cb.samples)

# available species
unique(unlist(sapply(unique(gsub(pattern = "\\|parent", replacement = "", x = ucsc.cb.samples$organisms)), function(x) {
  unlist(strsplit(x = x, split = ", "))
})))
```

<hr />

##### Summarise attributes

`scfetch` provides `StatDBAttribute` to summary attributes of [UCSC Cell Browser](https://cells.ucsc.edu/?#):

```r
StatDBAttribute(df = ucsc.cb.samples, filter = c("organism", "organ"), database = "UCSC")
```

<hr />

##### Extract metadata

`scfetch` provides `ExtractCBDatasets` to filter metadata with **collection**, **sub-collection**, **organ**, **disease status**, **organism**, **project** and **cell number** (The available values of these attributes can be obtained with `StatDBAttribute` except **cell number**). All attributes except cell number support fuzzy match with `fuzzy.match`, this is useful when selecting datasets.

```{r cb_extract, eval=FALSE}
hbb.sample.df <- ExtractCBDatasets(all.samples.df = ucsc.cb.samples, organ = c("brain", "blood"), organism = "Human (H. sapiens)", cell.num = c(1000, 2000))
```

<hr />

##### Extract cell type composition

`scfetch` provides `ExtractCBComposition` to extract cell type annotation and composition.

```r
hbb.sample.ct <- ExtractCBComposition(json.folder = json.folder, sample.df = hbb.sample.df)
```

<hr />

##### Load the online datasets to Seurat

After manually check the extracted metadata, `scfetch` provides `ParseCBDatasets` to **load the online count matrix** to Seurat. All the attributes available in `ExtractCBDatasets` are also same here. Please note that the loading process provided by `ParseCBDatasets` will load the online count matrix instead of downloading it to local. If multiple datasets available, users can choose to merge the SeuratObject with `merge`.

```r
hbb.sample.seu <- ParseCBDatasets(sample.df = hbb.sample.df)
```

<hr />

### Download object

`scfetch` provides functions for users to download processed single-cell RNA-seq data from [Zenodo](https://zenodo.org/), [CELLxGENE](https://cellxgene.cziscience.com/) and [Human Cell Atlas](https://www.humancellatlas.org/), including `RDS`, `RData`, `h5ad`, `h5`, `loom` objects. 

Until now, the public resources supported and the returned results:

| Resources        | URL                               | Download Type                          | Returned results        |
|------------------|-----------------------------------|----------------------------------------|-------------------------|
| Zenodo           | https://zenodo.org/               | count matrix, rds, rdata, h5ad, et al. | NULL or failed datasets |
| CELLxGENE        | https://cellxgene.cziscience.com/ | rds, h5ad                              | NULL or failed datasets |
| Human Cell Atlas | https://www.humancellatlas.org/   | rds, rdata, h5, h5ad, loom             | NULL or failed projects |

<hr />

#### Zenodo

[Zenodo](https://zenodo.org/) contains various types of processed objects, such as `SeuratObject` which has been clustered and annotated, `AnnData` which contains processed results generated by `scanpy`.

##### Extract metadata

`scfetch` provides `ExtractZenodoMeta` to extract dataset metadata, including dataset title, description, available files and corresponding md5. Please note that when the dataset is restricted access, the returned dataframe will be empty.

```r
# single doi
zebrafish.df <- ExtractZenodoMeta(doi = "10.5281/zenodo.7243603")

# vector dois
multi.dois <- ExtractZenodoMeta(doi = c("1111", "10.5281/zenodo.7243603", "10.5281/zenodo.7244441"))
```

<hr />

##### Download object

After manually check the extracted metadata, users can **download the specified objects** with `ParseZenodo`. The downloaded objects are controlled by `file.ext` and **the provided object formats should be in lower case (e.g. rds/rdata/h5ad).**

The returned result is a dataframe containing failed objects. If not `NULL`, users can re-run `ParseZenodo` by setting `doi.df` to the returned result.

```r
out.folder <- tempdir()
multi.dois.parse <- ParseZenodo(
  doi = c("1111", "10.5281/zenodo.7243603", "10.5281/zenodo.7244441"),
  file.ext = c("rdata", "rds"), out.folder = out.folder
)
```

<hr />

#### CELLxGENE

The [CELLxGENE](https://cellxgene.cziscience.com/) is a web server contains **910** single-cell datasets, users can explore, download and upload own datasets. The downloaded datasets provided by [CELLxGENE](https://cellxgene.cziscience.com/) have two formats: `h5ad (AnnData v0.8)` and `rds (Seurat v4)`.

##### Show available datasets

`scfetch` provides `ShowCELLxGENEDatasets` to extract dataset metadata, including dataset title, description, contact, organism, ethnicity, sex, tissue, disease, assay, suspension type, cell type, et al.

```r
# all available datasets
all.cellxgene.datasets <- ShowCELLxGENEDatasets()
```

<hr />

##### Summarise attributes

`scfetch` provides `StatDBAttribute` to summary attributes of [CELLxGENE](https://cellxgene.cziscience.com/):

```r
StatDBAttribute(df = all.cellxgene.datasets, filter = c("organism", "sex"), database = "CELLxGENE")
```

<hr />

##### Extract metadata

`scfetch` provides `ExtractCELLxGENEMeta` to filter dataset metadata, the available values of attributes can be obtained with `StatDBAttribute` except **cell number**:

```r
# human 10x v2 and v3 datasets
human.10x.cellxgene.meta <- ExtractCELLxGENEMeta(
  all.samples.df = all.cellxgene.datasets,
  assay = c("10x 3' v2", "10x 3' v3"), organism = "Homo sapiens"
)
```

<hr />

##### Download object

After manually check the extracted metadata, users can **download the specified objects** with `ParseCELLxGENE`. The downloaded objects are controlled by `file.ext` (choose from `"rds"` and `"h5ad"`).

The returned result is a dataframe containing failed datasets. If not `NULL`, users can re-run `ParseCELLxGENE` by setting `meta` to the returned result.

```r
out.folder <- tempdir()
ParseCELLxGENE(
  meta = human.10x.cellxgene.meta[1:5, ], file.ext = "rds",
  out.folder = out.folder
)
```

<hr />

### Format conversion

There are many tools have been developed to process scRNA-seq data, such as [Scanpy](https://scanpy.readthedocs.io/en/stable/), [Seurat](https://satijalab.org/seurat/), [scran](https://bioconductor.org/packages/release/bioc/html/scran.html) and [Monocle](http://cole-trapnell-lab.github.io/monocle-release/). These tools have their own objects, such as `Anndata` of `Scanpy`, `SeuratObject` of `Seurat`, `SingleCellExperiment` of `scran` and `CellDataSet`/`cell_data_set` of `Monocle2`/`Monocle3`. There are also some file format designed for large omics datasets, such as [loom](http://loompy.org/). To perform a comprehensive scRNA-seq data analysis, we usually need to combine multiple tools, which means we need to perform object conversion frequently. To facilitate user analysis of scRNA-seq data, `scfetch` provides multiple functions to perform object conversion between widely used tools and formats. The object conversion implemented in `scfetch` has two main advantages: 

* **one-step conversion between different objects**. There will be no conversion to intermediate objects, thus preventing unnecessary information loss.
* **tools used for object conversion are developed by the team of the source/destination object as far as possible**. For example, we use `SeuratDisk` to convert SeuratObject to loom, use `zellkonverter` to perform conversion between `SingleCellExperiment` and `Anndata`. When there is no such tools, we use `sceasy` to perform conversion.

<hr />

#### Test data

```r
# library
library(Seurat) # pbmc_small
library(scRNAseq) # seger
```

`SeuratObject`:

```r
# object
pbmc_small
```

`SingleCellExperiment`:
```r
seger <- scRNAseq::SegerstolpePancreasData()
```

<hr />

#### Convert SeuratObject to other objects

Here, we will convert SeuratObject to `SingleCellExperiment`, `CellDataSet`/`cell_data_set`, `Anndata`, `loom`.

##### SeuratObject to SingleCellExperiment

The conversion is performed with functions implemented in `Seurat`:
```r
sce.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", to = "SCE")
```

<hr />

##### SeuratObject to CellDataSet/cell_data_set

To `CellDataSet` (The conversion is performed with functions implemented in `Seurat`):

```r
# BiocManager::install("monocle") # reuqire monocle
cds.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", reduction = "tsne", to = "CellDataSet")
```

To `cell_data_set` (The conversion is performed with functions implemented in `SeuratWrappers`):

```r
# remotes::install_github('cole-trapnell-lab/monocle3') # reuqire monocle3
cds3.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", to = "cell_data_set")
```

<hr />

##### SeuratObject to AnnData

`AnnData` is a Python object, `reticulate` is used to communicate between Python and R. User should create a Python environment which contains `anndata` package and specify the environment path with `conda.path` to ensure the exact usage of this environment.

The conversion is performed with functions implemented in `sceasy`:
```r
# remove pbmc_small.h5ad first
anndata.file <- tempfile(pattern = "pbmc_small_", fileext = ".h5ad")
# you may need to set conda.path
ExportSeurat(
  seu.obj = pbmc_small, assay = "RNA", to = "AnnData",
  anndata.file = anndata.file
)
```

<hr />

##### SeuratObject to loom

The conversion is performed with functions implemented in `SeuratDisk`:
```r
loom.file <- tempfile(pattern = "pbmc_small_", fileext = ".loom")
ExportSeurat(
  seu.obj = pbmc_small, assay = "RNA", to = "loom",
  loom.file = loom.file
)
```

<hr />

#### Convert other objects to SeuratObject

##### SingleCellExperiment to SeuratObject

The conversion is performed with functions implemented in `Seurat`:
```r
seu.obj.sce <- ImportSeurat(obj = sce.obj, from = "SCE", count.assay = "counts", data.assay = "logcounts", assay = "RNA")
```

<hr />

##### CellDataSet/cell_data_set to SeuratObject

`CellDataSet` to `SeuratObject` (The conversion is performed with functions implemented in `Seurat`):
```r
seu.obj.cds <- ImportSeurat(obj = cds.obj, from = "CellDataSet", count.assay = "counts", assay = "RNA")
```

`cell_data_set` to `SeuratObject` (The conversion is performed with functions implemented in `Seurat`):
```{r cds2seu2, eval=FALSE}
seu.obj.cds3 <- ImportSeurat(obj = cds3.obj, from = "cell_data_set", count.assay = "counts", data.assay = "logcounts", assay = "RNA")
```

<hr />

##### AnnData to SeuratObject

`AnnData` is a Python object, `reticulate` is used to communicate between Python and R. User should create a Python environment which contains `anndata` package and specify the environment path with `conda.path` to ensure the exact usage of this environment.

The conversion is performed with functions implemented in `sceasy`:
```r
# you may need to set conda.path
seu.obj.h5ad <- ImportSeurat(
  anndata.file = anndata.file, from = "AnnData", assay = "RNA"
)
```

<hr />

##### loom to SeuratObject

The conversion is performed with functions implemented in `SeuratDisk` and `Seurat`:

```r
# loom will lose reduction
seu.obj.loom <- ImportSeurat(loom.file = loom.file, from = "loom")
```

<hr />

#### Conversion between SingleCellExperiment and AnnData

The conversion is performed with functions implemented in `zellkonverter`.

##### SingleCellExperiment to AnnData

```r
# remove seger.h5ad first
seger.anndata.file <- tempfile(pattern = "seger_", fileext = ".h5ad")
SCEAnnData(
  from = "SingleCellExperiment", to = "AnnData", sce = seger, X_name = "counts",
  anndata.file = seger.anndata.file
)
```

<hr />

##### AnnData to SingleCellExperiment

```r
seger.anndata <- SCEAnnData(
  from = "AnnData", to = "SingleCellExperiment",
  anndata.file = seger.anndata.file
)
```

<hr />

#### Conversion between SingleCellExperiment and loom

The conversion is performed with functions implemented in `LoomExperiment`.

##### SingleCellExperiment to loom

```r
# remove seger.loom first
seger.loom.file <- tempfile(pattern = "seger_", fileext = ".loom")
SCELoom(
  from = "SingleCellExperiment", to = "loom", sce = seger,
  loom.file = seger.loom.file
)
```

<hr />

##### loom to SingleCellExperiment

```r
seger.loom <- SCELoom(
  from = "loom", to = "SingleCellExperiment",
  loom.file = seger.loom.file
)
```

<hr />


## Contact
For any question, feature request or bug report please write an email to `songyb0519@gmail.com`.

<hr />

## Code of Conduct
  
Please note that the `scfetch` project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

<br />


