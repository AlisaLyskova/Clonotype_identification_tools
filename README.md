# Clonotype identification tools
Evaluating TRUST4 and a custom MiXCR pipeline against the standard MiXCR workflow: IGV‑validated comparison 

## Analysis Workflow Summary

This project compares immune repertoire analysis tools by validating results against visual inspection in IGV. Below is the step‑by‑step workflow.

### 1. Input Data & Reference Selection
- We started with major clones identified by the standard MiXCR pipeline.

### 2. IGV‑Based Validation of Ig Genes
- For each major clone, we inspected the corresponding BAM files in IGV (Integrative Genomics Viewer).
- Goal: confirm the presence and structure of immunoglobulin (Ig) genes (e.g., V/D/J segments) supporting each clonotype.
- Criteria checked:
  - Read coverage across V(D)J regions
  - Splice junctions and alignment quality
  - Absence of artifacts or misassemblies


### 3. Re‑Analysis with Alternative Tools

Using the same input data (BAM files), we ran to compare against the standard pipeline’s output:
- **TRUST4**
- **MiXCR** (custom run)
- **Vidjil**


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

see https://github.com/liulab-dfci/TRUST4

We ran in two modes: with BAM file and FASTQ files

<pre>
./run-trust4 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa -1 example/example_1.fq -2 example/example_2.fq -o TRUST_example
</pre>

<pre>
./run-trust4 -b example/example.bam -f hg38_bcrtcr.fa --ref human_IMGT+C.fa
</pre>

where hg38_bcrtcr.fa is created using the command (or a file from TRUST4 can be used)
<pre>
perl BuildDatabaseFa.pl reference.fa annotation.gtf bcr_tcr_gene_name.txt > hg38_bcrtcr.fa
</pre>

- bcr_tcr_gene_name.txt - file with genes list (a file from TRUST4 was used)
- annotation.gtf - GENCODE annotation file (v43, GRCh38) in gtf format

Then we selected only funcitonal clonotypes with the command:
<pre>
cat "${sample}_cdr3.out" | awk '$NF == 1' > "${sample}_productive_cdr3.out"
perl trust-simplerep.pl "${sample}_productive_cdr3.out" > "${sample}_report_clean.tsv"
</pre>

#### Vidjil run

To run vidjil-algo we first combined reads in one file with the command:

<pre>
seqtk mergepe $R1 $R2 | pigz -c -p 8 > $merged_reads
</pre>

Then we ran vidjil-algo with the command:

<pre>
vidjil-algo -g germline/homo-sapiens.g -o vidjil_res $merged_reads
</pre>

### 4. Comparison & Validation
We then:
- Cross‑referenced *major clones* from the standard MiXCR with **IGV visualizations**
- Compared clonotypes from **TRUST4**, **Vidjil** and **custom MiXCR** with the original *major clones* from the standard MiXCR

### 5. Results

The comparison information is provided in the table [results/samples_clones.xlsx](results/samples_clones.xlsx)

The precision score is shown in the table [results/mixcr_igv_trust4_statistics.csv](results/mixcr_igv_trust4_statistics.csv)

The figure below shows the precision of the predicted major clonotypes using the standard MiXCR pipeline relative to those found in IGV, as well as the precision of the detected clonotypes using the TRUST4  running in two modes (with a bam file and with fastq files) and custom MiXCR command relative the predicted major clonotypes from the standard MiXCR pipeline.

<a href="images/heatmap.png" target="_blank">
  <img width="1500" height="2100" src="images/heatmap.png" alt="heatmap.png">
</a>

**Overall precision assessment**
Comparison groups | Precision, %
--- | ---
Precision BCR genes from the standard MiXCR pipeline relative to those found in IGV | 62.5
Precision TCR genes from the standard MiXCR pipeline relative to those found in IGV | 83.1

**Metrics Comparison Across Tools (relative proven in IGV clonotypes from the standard MiXCR pipeline)**

*The term Conjunction denotes the logical conjunction (AND operation) of the results from three underlying methods. The term Disjunction - the logical disjunction*

| metric , %  | Accuracy | F1   | Precision | Recall | Specificity |
|-------------|----------|------|-----------|--------|-------------|
| Conjunction | 41.38    | 29.01 | 48.85     | 23.74  | 56.03       |
| Disjunction | 73.10    | 79.37 | 87.87     | 80.06  | 37.07       |
| MiXCR       | 60.69    | 66.61 | 83.62     | 62.76  | 37.07       |
| TRUST4_BAM  | 51.03    | 45.68 | 66.09     | 38.45  | 51.44       |
| Vidjil      | 61.38    | 59.72 | 80.46     | 52.41  | 52.01       |

We also assessed the dependence of precision on sample coverage

<a href="images/scatterplot.png" target="_blank">
  <img width="1000" height="1500" src="images/scatterplot.png" alt="scatterplot.png">
</a>


#### Key Findings

1. **IGV Validation**  
   - Not all clonotypes identified by the standard MiXCR pipeline were confirmed via IGV visual inspection.  

2. **MiXCR Performance**  
   - No significant difference observed between BCR and TCR detection using the standard MiXCR pipeline.  

3. **Tools comparison**
   - Disjunction leads with the highest F1‑score (79.37 %), showing the best balance of precision and recall.
   - The highest precision is achieved by MiXCR (83.62 %). This tools produce fewer false positives, which is critical in tasks where the cost of an error is high.
   - TRUST4 is the least effective
     - F1‑score: 45.68 % (lowest)
     - Recall: 38.45 % (low)
     - Accuracy: 51.03 % (below average)

5. **TRUST4 Modes Comparison**  
   - BAM and FASTQ modes produced nearly identical results.  
   - BAM mode ran **twice as fast** as FASTQ mode.
  
6. **Coverage vs. Precision Trend**  
   - There is a clear tendency: **higher coverage correlates with improved precision in clonotype detection**.  
  
### Summary
- **MiXCR** (standard/custom) remains the most reliable tool for clonotype detection.  
- **TRUST4** in BAM mode offers a faster alternative with slightly lower accuracy.  
- **IGV validation** is strongly recommended for critical clonotypes to ensure confidence in results.
- **Coverage matters**: Higher sequencing coverage generally leads to more precision clonotype detection.

### Scripts
Script path | Description
--- | ---
[main.sh](main.sh) | Runnig the programs MiXCR and TRUST4 for samples and adding found clones to the xlsx file
[stat_table_plot.py](stat_table_plot.py) | Creates table with precision from [results/samples_clones.xlsx](results/samples_clones.xlsx), files with samples ids and samples coverage are required and creates heatmap and scatter plot
