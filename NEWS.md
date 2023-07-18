# scfetch 0.2.0

## New features
* Provided APIs for GEO.
* Added `Read10XOnline` to load cellranger output.
* Added `ExtractCBComposition` to extract cell type composition of UCSC Cell Browser datasets.

-------------------

## Minor changes
* Supported lazy mode in `ShowCBDatasets` to save time for multi-runs.
* Optimized the paras for `ExtractCBDatasets` and `ParseCBDatasets`.
* Added `cell.num` filter to `ParseCBDatasets`, `ExtractCBDatasets` and `ShowPanglaoDBMeta`.
* Extracted matrix, barcode and feature info from dataset.json (comprehensive) instead of desc.json.
* Fixed bugs in `ParseCBDatasets` when dealing with cellranger output.
* Added `timeout` in `ParseCBDatasets` to avoid possible timeout error.

-------------------

# scfetch 0.1.0

* Provided APIs for Zenodo, PanglaoDB and UCSC Cell Browser.
