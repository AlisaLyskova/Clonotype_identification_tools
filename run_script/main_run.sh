#!/usr/bin/env bash

source ./vdj_clonotypes_config.cfg

function process_vdj_tsv()
{
  python3 -c '
import csv
import sys
import itertools
import pandas as pd
import re
import numpy as np

# mixcr_mode: 1 - process mixcr tsv file, 0 - process tsv files from other tools
# separator: separator for splitting columns to new rows ex. "," for mixcr and vidjil tsv file, "|" for other tools
def process_vdj_tsv(inp, outp, mixcr_mode, separator):
    with open(inp, newline="", encoding="utf-8") as f_in, open(outp, "w", newline="", encoding="utf-8") as f_out:
        r = csv.reader(f_in, delimiter="\t")
        w = csv.writer(f_out, delimiter="\t", lineterminator="\n")

        header = next(r, None)
        if header is not None:
            w.writerow(header)

        # split columns with commas to new rows
        for row in r:
            if not row:
                continue
            parts_per_col = [[x.strip() for x in col.split(separator)] if separator in col else [col] for col in row]
            for expanded in itertools.product(*parts_per_col):
                w.writerow(expanded)

    df = pd.read_csv(outp, delimiter="\t")
    if mixcr_mode == 1:
        new_columns_names = ["MiXCR_V_score", "MiXCR_D_score", "MiXCR_J_score"]
        for i in range(3):
            col_name = df.columns[i]
            new_col_name = new_columns_names[i]

            # add score value to new column
            df[new_col_name] = df[col_name].apply(
                lambda cell: "|".join(re.findall(r"\(([^)]*)\)", str(cell)))
                if pd.notna(cell) else float("nan")
            )

    # remove allele and score value
    def process_gene_string(cell_value):
        if pd.isna(cell_value) or cell_value == float("nan"):
            return float("nan")
        cleaned = re.sub(r"\*.*", "", str(cell_value))
        return cleaned


    df.iloc[:, :3] = df.iloc[:, :3].applymap(process_gene_string)

    # Convert columns 4, 5, 6, 7 to float type if they exist
    for col_idx in range(3, 7):
        if col_idx < len(df.columns):
            df[df.columns[col_idx]] = pd.to_numeric(df[df.columns[col_idx]], errors="coerce")

    # sort by all columns
    df = df.sort_values(by=list(df.columns), ascending=False)

    # drop duplicates using first 3 columns as key (0,1,2)
    df = df.drop_duplicates(subset=df.columns[:3], keep="first")

    df.to_csv(outp, sep="\t", index=False)


process_vdj_tsv(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
' $1 $2 $3 $4 
}
export -f process_vdj_tsv


function mixcr_run()
{
  sample=$1
  work_dir=$2
  out_dir="${work_dir}/mixcr"
  mkdir $out_dir
  log=$3
  
  echo "RUN MiXCR $(date)" >> $log

  R1=${work_dir}/${sample}_R1.fastq.gz
  R2=${work_dir}/${sample}_R2.fastq.gz

  source activate $mixcr_env
  mixcr analyze exome-seq -f --species hs -t 10 $R1 $R2 \
    "${out_dir}/${sample}" > $log 2>&1
  
  conda activate base
    
  # merge tsv files
  res_file_tmp="${out_dir}/${sample}_mixcr_clonotypes_tmp.tsv"
  res_file="${out_dir}/${sample}_mixcr_clonotypes.tsv"
  echo -e "V\tD\tG\tMiXCR_readcount" > $res_file_tmp

  check=$(find $out_dir -name "*clones*.tsv")
  if [[ -z $check ]]; then
    echo "No tsv files from MiXCR for $sample" >> $log
  fi

  for file in $(find $out_dir -name "*clones*.tsv")
  do
    cat $file | sed '1d' | awk -F "\t" -v OFS="\t" '{print $6,$7,$8,$2}' >> $res_file_tmp
  done 
  
  ## process tsv file to format V, D, J, readcount, score for v, d, j, genes, delete duplicates and split to rows with commas
  process_vdj_tsv $res_file_tmp $res_file 1 ","
  rm $res_file_tmp 
}
export -f mixcr_run


function vidjil_run()
{
  sample=$1
  work_dir=$2
  out_dir="${work_dir}/vidjil"
  mkdir $out_dir
  log=$3

  echo "RUN Vidjil $(date)" >> $log

  R1=${work_dir}/${sample}_R1.fastq.gz
  R2=${work_dir}/${sample}_R2.fastq.gz
  
  # merge reads
  merged_reads=${work_dir}/${sample}_reads_merged.fastq.gz
  $seqtk mergepe $R1 $R2 | pigz -c -p 8 > $merged_reads

  # run vidjil
  $vidjil -g "${vidjil_dir}/germline/homo-sapiens.g" -o "${out_dir}" $merged_reads > $log 2>&1
  
  vidjil_res="${out_dir}/${sample}_reads_merged.fastq.tsv"
  res_file_tmp="${out_dir}/${sample}_vidjil_clonotypes_tmp.tsv"
  res_file="${out_dir}/${sample}_vidjil_clonotypes.tsv"

  echo -e "V\tD\tG\tVidjil_readcount" > $res_file_tmp
  cat $vidjil_res | sed '1d' | awk -v OFS="\t" -F "\t" '{print $3,$4,$5,$2}' | awk -F "\t" 'NR > 1 {count = 0; if ($2 != "") count++; if ($3 != "") count++; if ($4 != "") count++; if (count >= 2) print $0}' >> $res_file_tmp
  process_vdj_tsv $res_file_tmp $res_file 0 ","
  rm $res_file_tmp 
}
export -f vidjil_run


function trust4_run()
{
  sample=$1
  work_dir=$2
  out_dir="${work_dir}/trust4"
  mkdir $out_dir
  log=$3
  
  echo "RUN TRUST4 $(date)" >> $log

  bam="${work_dir}/${sample}_roi.processed.bed.bam"
  
  $trust4 -b $bam -f "${trust4_dir}/bcrtcr_gencode.v49.fa" --ref "${trust4_dir}/human_IMGT+C.new.fa" --od $out_dir -o $sample -t $per_sample_threads
  
  # filtering unfunctional clones
  cat "${out_dir}/${sample}_cdr3.out" | awk '$NF == 1' > "${out_dir}/${sample}_productive_cdr3.out"
  perl "${trust4_dir}/trust-simplerep.pl" "${out_dir}/${sample}_productive_cdr3.out" > "${out_dir}/${sample}_report_clean.tsv"
  
  res_file="${out_dir}/${sample}_trust4_clonotypes.tsv"
  res_file_tmp="${out_dir}/${sample}_trust4_clonotypes_tmp.tsv"
  echo -e "V\tD\tG\tTRUST4_readcount" > $res_file_tmp
  cat "${out_dir}/${sample}_report_clean.tsv" | sed '1d' | awk -v OFS="\t" -F "\t" '{print $5,$6,$7,$1}' >> $res_file_tmp
  process_vdj_tsv $res_file_tmp $res_file 0 "|"
  rm $res_file_tmp 
}
export -f trust4_run


function detect_clonotypes()
{
  sample_bam=$1
  new_bam="${SAMPLE_DIR}/${SAMPLE_NAME}_processed.bam"
  LOG="${SAMPLE_DIR}/${SAMPLE_NAME}_VDJ_clonotypes.log"
  RUN_LOG="${SAMPLE_DIR}/${SAMPLE_NAME}_VDJ_clonotypes_run.log"
  touch "$LOG"
  touch "$RUN_LOG"

  # files preprocess
  echo "INFO: RUN files preprocess $(date)" >> $RUN_LOG
  R1=${SAMPLE_DIR}/${SAMPLE_NAME}_R1.fastq.gz
  R2=${SAMPLE_DIR}/${SAMPLE_NAME}_R2.fastq.gz
  U=${SAMPLE_DIR}/${SAMPLE_NAME}_U.fastq.gz

  vdj_bam="${SAMPLE_DIR}/${SAMPLE_NAME}_VDJ.bam"
  # 1. Get reads overlapping VDJ regions
  $samtools view -P --threads $per_sample_threads -M -L $vdj_regions $sample_bam -o $vdj_bam 2>> "$LOG"

  # 2. Extract XA/SA coordinates from those reads, keep only ones in VDJ regions
  $samtools view $vdj_bam 2>> "$LOG" \
    | awk 'BEGIN{OFS="\t"} {
        for(i=12;i<=NF;i++) {
          if($i~/^XA:Z:/) {
            split(substr($i,6), alns, ";")
            for(j in alns) { split(alns[j], f, ","); if(f[1]!="") print f[1], f[2]<0?-f[2]:f[2] }
          }
          if($i~/^SA:Z:/) {
            split(substr($i,6), alns, ";")
            for(j in alns) { split(alns[j], f, ","); if(f[1]!="") print f[1], f[3] }
          }
        }
      }' \
    | awk 'NR==FNR{chr[$1]=$1; start[$1]=$2; end[$1]=$3; next}
          ($1 in chr) && ($2 >= start[$1]) && ($2 <= end[$1])' \
      $vdj_regions - \
    | awk '{print $1":"$2"-"($2+1)}' \
    | sort -u > "${SAMPLE_DIR}/xa_sa_regions.txt" 2>> "$LOG"

  # 3. Extract those reads from the original BAM
  xargs -a "${SAMPLE_DIR}/xa_sa_regions.txt" $samtools view -P --threads $per_sample_threads -b $sample_bam > "${SAMPLE_DIR}/xa_sa.bam" 2>> "$LOG"

  # 4. Merge with VDJ reads, sort and index
  $samtools merge -f "${SAMPLE_DIR}/merged.bam" "$vdj_bam" "${SAMPLE_DIR}/xa_sa.bam" 2>> "$LOG"
  $samtools sort --threads $per_sample_threads -o $new_bam "${SAMPLE_DIR}/merged.bam" 2>> "$LOG"
  $samtools index $new_bam 2>> "$LOG"

  echo "INFO: RUN picard $(date)" >> $RUN_LOG
  java -jar $picard SamToFastq -I $new_bam -FU $U  -F $R1 -F2 $R2 --VALIDATION_STRINGENCY SILENT >> "$LOG" 2>&1

  # run vdj tools
  echo "INFO: RUN VDJ tools parallel $(date)" >> $RUN_LOG
  
  parallel -j 3 ::: \
  "trust4_run $SAMPLE_NAME $SAMPLE_DIR $LOG" \
  "vidjil_run $SAMPLE_NAME $SAMPLE_DIR $LOG" \
  "mixcr_run $SAMPLE_NAME $SAMPLE_DIR $LOG" 2>> "$LOG"
  
  python3 -c '
import sys
import pandas as pd

trust4_file = sys.argv[1]
vidjil_file = sys.argv[2]
mixcr_file = sys.argv[3]
outfile = sys.argv[4]

trust4_df = pd.read_csv(trust4_file, sep="\t").fillna(".")
vidjil_df = pd.read_csv(vidjil_file, sep="\t").fillna(".")
mixcr_df = pd.read_csv(mixcr_file, sep="\t").fillna(".")

key_columns = trust4_df.columns[:3].tolist()

result = trust4_df.merge(vidjil_df, on=key_columns, how="outer").merge(mixcr_df, on=key_columns, how="outer")

readcount_cols = [c for c in result.columns if c.endswith("_readcount")]

def detect_tools(row):
    tools = []
    for col in readcount_cols:
        val = row[col]
        if pd.notna(val) and val != ".":
            tools.append(col.replace("_readcount", ""))
    return ",".join(sorted(tools)) if tools else "."

result["tools"] = result.apply(detect_tools, axis=1)

result.to_csv(outfile, sep="\t", index=False, na_rep=".")
' "${SAMPLE_DIR}/trust4/${SAMPLE_NAME}_trust4_clonotypes.tsv" "${SAMPLE_DIR}/vidjil/${SAMPLE_NAME}_vidjil_clonotypes.tsv" "${SAMPLE_DIR}/mixcr/${SAMPLE_NAME}_mixcr_clonotypes.tsv" "${SAMPLE_DIR}/${SAMPLE_NAME}_clonotypes.tsv" 2>> "$LOG"

  if [ -f "${SAMPLE_DIR}/${SAMPLE_NAME}_clonotypes.tsv" ]; then
    echo "INFO: file with clonotypes was created" >> $RUN_LOG
  else
    echo "ERROR: file with clonotypes was not created" >> $RUN_LOG
  fi
}
export -f detect_clonotypes

detect_clonotypes $BAM
