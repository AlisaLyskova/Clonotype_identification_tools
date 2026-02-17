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

Using the same input data (BAM files), we ran:
- **TRUST4**: to reconstruct BCR clonotypes independently
- **MiXCR** (custom run): to compare against the standard pipeline’s output


#### BAM File Preparation for TRUST4 and MiXCR

To focus the analysis on immunoglobulin gene regions, we performed targeted extraction using GENCODE annotations:

Region selection:
- Extracted genomic intervals corresponding to immunoglobulin genes (e.g., IGH, IGK, IGL) from the GENCODE annotation file (v43, GRCh38).
- Generated a BED file listing coordinates of all Ig‑related exons and introns.

Targeted BAM subsetting:
- Used samtools view to filter the original BAM files, retaining only reads mapped to Ig regions:
<pre>
samtools view -b -L ig_regions.bed input.bam --write-index -o ig_subset.bam
</pre>

FASTQ conversion:

Converted the filtered BAM file to FASTQ format using Picard’s SamToFastq tool (v2.27.4):

<pre>
java -jar picard.jar SamToFastq -I input.bam -F input_R1.fastq.gz -F2 input_R2.fastq.gz --VALIDATION_STRINGENCY SILENT
</pre>

#### MiXCR run

We used the standard  command in exome-seq mode

<pre>
mixcr analyze exome-seq -f --species hs -t 10 input_R1.fastq.gz input_R2.fastq.gz dir/sample_name
</pre>

#### TRUST4 run

We ran in two modes: with BAM file and FASTQ files

<pre>
./run-trust4 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa -1 example/example_1.fq -2 example/example_2.fq -o TRUST_example
</pre>

<pre>
./run-trust4 -b example/example.bam -f hg38_bcrtcr.fa --ref human_IMGT+C.fa
</pre>

where hg38_bcrtcr.fa is created using the command (or you can use given by TRUST4)
perl BuildDatabaseFa.pl reference.fa annotation.gtf bcr_tcr_gene_name.txt > bcrtcr.fa
bcr_tcr_gene_name.txt - file with genes list (you can use given by TRUST4)
annotation.gtf - file with annotation in gtf format


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
<a href="images/heatmap.png" target="_blank">
  <img src="images/heatmap.png" alt="heatmap.png">
</a>


Precision, %
Precision_BCR                    62.5
Precision_TCR                    83.1
Precision_TRUST4_BAM             32.1
Precision_TRUST4_FASTQ           32.1
Precision_MiXCR                  55.0

Зависимость от покрытия
<a href="images/scatterplot.png" target="_blank">
  <img src="images/scatterplot.png" alt="scatterplot.png">
</a>
