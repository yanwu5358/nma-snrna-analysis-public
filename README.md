# NMA single-nucleus RNA-seq analysis: draft

This repository contains preprocessing, integration and differential-expression analysis scripts. Sample identifiers in the integration script are illustrative placeholders.

## Contents

- `R/CellBender_B4_AF_Compare.r`: exploratory comparison of raw and CellBender-corrected expression and marker plots. It contains interactive fragments that require existing objects and settings; it is not a standalone script.
- `R/01_integrate_harmony.R`: loads per-sample cleaned Seurat objects, applies SCTransform, merges samples, runs PCA and Harmony, and produces clustering and UMAP diagnostics. Edit its example sample manifest and input paths for your data.
- `R/motor_neuron_gate.R`: motor-neuron selection rule applied to the original SCT corrected counts.
- `R/compartment_MAST_GO_adapted.R`: implementation of the stated compartment MAST/GO settings for each cell compartment.
- `R/mn_SCT_MAST_GO_GSEA.R`: marker-defined MN, SCT-based mixed-model MAST, GO BP and signed-z GSEA.&#x20;

The intended workflow includes quality control, per-sample SCTransform, Harmony, clustering and annotation, cell-compartment MAST/GO, and marker-defined motor-neuron MAST/GSEA. Raw expression data, private sample IDs and local paths are not included.&#x20;

The integration script defaults to Harmony grouping by sample, 2,000 features and 30 PCs; Its sample-specific filenames and manifest are examples and must be configured for a new dataset.

Dependencies for the marker function: Seurat, SeuratObject and Matrix. The code expects an already annotated object with `celltype_res05`, an SCT assay and the SCT corrected counts layer, before any operation that changes these counts.

The compartment script calls `PrepSCTFindMarkers` on the selected cells before MAST and uses the SCT data layer. This preparation may also require the original UMI assay referenced by the fitted SCT models; expression testing still uses SCT. Save motor-neuron labels before preparation because SCT counts can change.&#x20;

## Motor-neuron analysis

Run `Rscript R/mn_SCT_MAST_GO_GSEA.R INPUT.rds OUTPUT_DIR [SAMPLE_COLUMN]`; the optional sample column defaults to `sample`. The object must contain `condition` (CTRL, FUS, TDP), `celltype_res05`, sample metadata, and the original SCT counts and data layers.&#x20;

Each comparison filters genes detected in at least 10% of cells in either group. Cellular detection rate is the standardized number of detected genes after that filter, calculated independently for each comparison. A glmer random-intercept hurdle model tests condition, with the H-component LRT P value used for DEG BH adjustment; missing tests remain missing and the number of filtered genes is the BH multiplicity. The MAST `logFC` component supplies coefficient and z-score. GO BP uses BH <=0.05 DEG lists split by coefficient sign, minimum five input genes, a background of genes with finite hurdle tests, and GO nominal P <=0.05 and qvalue <=0.10. GSEA uses all finite logFC z-scores, independent of DEG thresholds, and 10–500 gene set sizes before and after intersection. `fgsea` and `msigdbr` versions are recorded; MSigDB data may require a first-time network download.&#x20;

