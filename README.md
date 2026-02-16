# Clonotype identification tools
Evaluating TRUST4 and a custom MiXCR pipeline against the standard MiXCR workflow: IGV‑validated comparison 

Here’s a clear, professional `README.md` section for GitHub summarizing your workflow:


## Analysis Workflow Summary

This project compares immune repertoire analysis tools by validating results against visual inspection in IGV. Below is the step‑by‑step workflow.

### 1. Input Data & Reference Selection
- We started with **major clones** (dominant clonotypes) identified by the **standard MiXCR pipeline**.
- These clones represent high‑frequency B‑cell receptor (BCR) sequences in the sample.

### 2. IGV‑Based Validation of Ig Genes
- For each major clone, we inspected the corresponding **BAM files** in **IGV (Integrative Genomics Viewer)**.
- Goal: confirm the presence and structure of **immunoglobulin (Ig) genes** (e.g., V/D/J segments) supporting each clonotype.
- Criteria checked:
  - Read coverage across V(D)J regions
  - Splice junctions and alignment quality
  - Absence of artifacts or misassemblies


### 3. Re‑Analysis with Alternative Tools
Using the same input data (raw reads or BAM files), we ran:
- **TRUST4**: to reconstruct BCR clonotypes independently
- **MiXCR** (custom run): to compare against the standard pipeline’s output


### 4. Comparison & Validation
We then:
- Compared clonotypes from **TRUST4** and **custom MiXCR** with the original *major clones* from the standard MiXCR
- Cross‑referenced all calls with **IGV visualizations** to assess:
  - Concordance of V/D/J assignments
  - Accuracy of CDR3 boundaries
  - Support by aligned reads


### Key Objectives
- Evaluate how TRUST4 and a custom MiXCR run perform relative to the standard MiXCR pipeline
- Use IGV as a manual validation layer to identify potential false positives/negatives
- Understand tool‑specific biases in clonotype inference


### Outputs
Results include:
- Lists of clonotypes from each tool
- IGV screenshots highlighting key regions
- Concordance tables (overlap between tools)
- Notes on discrepancies and their visual validation status


### Results
<img width="4000" height="6000" alt="heatmap_v2" src="https://github.com/user-attachments/assets/a3ec288b-00d3-4556-a30d-3f5590905823" />

Precision, %
Precision_BCR                    62.5
Precision_TCR                    83.1
Precision_TRUST4_BAM             32.1
Precision_TRUST4_FASTQ           32.1
Precision_MiXCR                  55.0
