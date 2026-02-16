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

  mixcr analyze exome-seq -f --species hs -t 10 $R1 $R2 \
    "${work_dir}/${mixcr_res_dir}/${SAMPLE}/${SAMPLE}" > "${work_dir}/${mixcr_res_dir}/${SAMPLE}/${SAMPLE}.log" 2>&1
  BAM File Preparation
To focus the analysis on immunoglobulin (Ig) gene regions, we performed targeted extraction using GENCODE annotations:

Region selection:

Extracted genomic intervals corresponding to immunoglobulin genes (e.g., IGH, IGK, IGL) from the GENCODE annotation file (v43, GRCh38).

Generated a BED file listing coordinates of all Ig‑related exons and introns.

Targeted BAM subsetting:

Used samtools view to filter the original BAM files, retaining only reads mapped to Ig regions:

bash
samtools view -b -L ig_regions.bed input.bam -o ig_subset.bam
Indexed the subsetted BAM:

bash
samtools index ig_subset.bam
FASTQ conversion:

Converted the filtered BAM file to FASTQ format using Picard’s SamToFastq tool (v2.27.4):

bash
java -jar picard.jar SamToFastq \
  I=ig_subset.bam \
  O=output_R1.fastq \
  O2=output_R2.fastq
This produced paired‑end FASTQ files (output_R1.fastq, output_R2.fastq) for downstream analysis.


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
