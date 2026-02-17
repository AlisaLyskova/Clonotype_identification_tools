export work_dir=""
export picard=
export ref=
export threads=

export trust4_dir=""
export trust4=="${trust4_dir}/run-trust4"

export ref_trust="${trust4_dir}/human_IMGT+C.fa"
export bcrtcr_gencode_fa="${trust4_dir}/bcrtcr_gencode.v49.fa"

export res_xlsx="samples_clones.xlsx"
export mixcr_pipeline_res="mixcr_pipeline_results"



function create_suppl_files_trust4()
{
	perl "${trust4_dir}/BuildDatabaseFa.pl" $ref $gencode "${trust4_dir}/human_vdjc.list" > $bcrtcr_gencode_fa
}
export -f create_suppl_files_trust4


function run_trust4()
{
	sample=$1
	sample_bam="${sample}.bam"
	trust4_output_dir_bam=${work_dir}/${sample}_trust4_bam
        trust4_output_dir_fastq=${work_dir}/${sample}_trust4_fastq
	R1=${work_dir}/${sample}_R1.fastq.gz
	R2=${work_dir}/${sample}_R2.fastq.gz

	java -jar $picard SamToFastq -I $sample_bam -F $R1 -F2 $R2 --VALIDATION_STRINGENCY SILENT

	mkdir -p $trust4_output_dir_fastq
	perl $trust4 -1 $R1 -2 $R2 -f ${trust4_dir}/hg38_bcrtcr.fa --ref $ref_trust --od $trust4_output_dir_fastq -o $sample -t $threads
	# filtering unfunctional clones
	cat "${trust4_output_dir_fastq}/${sample}_cdr3.out" | awk '$NF == 1' > "${trust4_output_dir_fastq}/${sample}_productive_cdr3.out"
	perl "${trust4_dir}/trust-simplerep.pl" "${trust4_output_dir_fastq}/${sample}_productive_cdr3.out" > "${trust4_output_dir_fastq}/${sample}_report_clean.tsv"

	# run with bam file
        mkdir -p $trust4_output_dir_bam
	$trust4 -b $roi_processed_bam -f $bcrtcr_gencode_fa --ref $ref_trust_2 --od $trust4_output_dir_bam -o $sample -t $threads

        # filtering unfunctional clones
        cat "${trust4_output_dir_bam}/${sample}_cdr3.out" | awk '$NF == 1' > "${trust4_output_dir_bam}/${sample}_productive_cdr3.out"
	perl "${trust4_dir}/trust-simplerep.pl" "${trust4_output_dir_bam}/${sample}_productive_cdr3.out" > "${trust4_output_dir_bam}/${sample}_report_clean.tsv"
}
export -f run_trust4


function run_mixcr()
{
        sample=$1
        R1=${work_dir}/${sample}_R1.fastq.gz
        R2=${work_dir}/${sample}_R2.fastq.gz

	source activate mixcr_env
	mkdir "${work_dir}/${sample}/mixcr_res"
	mixcr analyze exome-seq -f --species hs -t 10 $R1 $R2 \
		"${work_dir}/${sample}/mixcr_res/${sample}" > "${work_dir}/${sample}/mixcr_res/${sample}.log" 2>&1
}
export -f run_mixcr


function merge_results_mixcr()
{
  SAMPLE=$1
  sample_res_dir="${work_dir}/${SAMPLE}/mixcr_res"
  res_file="${sample_res_dir}/${SAMPLE}_results_merged.tsv"
  res_file_2="${sample_res_dir}/${SAMPLE}_results_merged_2.tsv"

  check=$(find $sample_res_dir -name "*.tsv")
  if [ -z $check ]; then
    echo $SAMPLE
    return 1
  fi

  for file in $(find $sample_res_dir -name "*.tsv")
  do
    cat $file | awk -F "\t" -v OFS="\t" '{print $2,$6,$7,$8,$9}' | head -n 1 > $res_file
    cat $file | awk -F "\t" -v OFS="\t" '{print $2,$6,$7,$8,$9}' | sed '1d' >> $res_file_2
  done

  cat $res_file_2 | sort -n -r -k 1 >> $res_file
  rm $res_file_2

}
export -f merge_results_mixcr


function add_res_to_xlsx()
{
	sample=$1
        res_tsv_trust_bam="${work_dir}/${sample}/${sample}_trust4_bam/${sample}_report_clean.tsv"
        res_tsv_trust_fastq="${work_dir}/${sample}/${sample}_trust4_fastq/${sample}_report_clean.tsv"
	res_tsv_mixcr="${work_dir}/${sample}/mixcr_res/${sample}_results_merged.tsv"

	# creating file with clones only
	cat $res_tsv_trust_bam | sed '1d' | awk -v OFS="\t" -F "\t" '{print $1, $5,$6,$7}' > "${work_dir}/${sample}/${sample}_trust4_bam/${sample}_clones.csv"
        cat $res_tsv_trust_fastq | sed '1d' | awk -v OFS="\t" -F "\t" '{print $1, $5,$6,$7}' > "${work_dir}/${sample}/${sample}_trust4_fastq/${sample}_clones.csv"

	# count common values in mixcr and trust
	id=$(cat "${mixcr_pipeline_res}/all_mixcr_vs_genome.csv" | grep "$sample" | awk -F, '{print $2}')
	echo $id
	pipeline_res="${mixcr_pipeline_res}/${id}_mixcr.txt"
	if [ ! -f $pipeline_res ]; then
		echo "No mixcr pipeline results for $sample"
		return 1
	fi

	new_column_trust_bam=(TRUST_BAM)
        new_column_trust_fastq=(TRUST_FASTQ)
        new_column_mixcr=(MiXCR)
	for row in {1..5}
	do
		expr_vdj=$(cat $pipeline_res | awk -F "\t" '{print $5,$6,$7}' | awk '{for(i=1;i<=3;i++) sub(/\*.*/, "", $i); print}' | sed '1d' | awk -F " " -v OFS=".*" '{print $1,$2,$3}' | head -n $row | tail -n 1)
		check_1=$(cat $res_tsv_trust_bam | grep -E "$expr_vdj")
                check_2=$(cat $res_tsv_trust_fastq | grep -E "$expr_vdj")
                check_3=$(cat $res_tsv_mixcr | grep -E "$expr_vdj")
		if [[ ! -z $check_1 ]]; then 
			new_column_trust_bam+=("+")
		else
                        new_column_trust_bam+=("-")
		fi
                if [[ ! -z $check_2 ]]; then
                        new_column_trust_fastq+=("+")
                else
                        new_column_trust_fastq+=("-")
                fi
                if [[ ! -z $check_3 ]]; then
                        new_column_mixcr+=("+")
                else
                        new_column_mixcr+=("-")
                fi
	done

#<<"xlsx"
	python3 -c '
import sys
import pandas as pd
from openpyxl import load_workbook

file_path = sys.argv[1]

wb = load_workbook(file_path)
sheet_name = sys.argv[2]
clones_file_trust_bam = sys.argv[3]
clones_file_trust_fastq = sys.argv[4]
clones_file_mixcr = sys.argv[5]

new_column_trust_bam = sys.argv[6:12]
new_column_trust_fastq = sys.argv[12:18]
new_column_mixcr = sys.argv[18:24]

clones_df = pd.read_csv(clones_file_trust_bam, sep="\t", header=None)
clones_trust_bam = clones_df.values.tolist()

clones_df = pd.read_csv(clones_file_trust_fastq, sep="\t", header=None)
clones_trust_fastq = clones_df.values.tolist()

clones_df = pd.read_csv(clones_file_mixcr, sep="\t").fillna(".")
clones_df = clones_df.iloc[:, :-1]
clones_mixcr = clones_df.values.tolist()

ws = wb[sheet_name]

# add clones from tools
## TRUST BAM
last_row = ws.max_row

### add new header - source of clones
ws.cell(row=last_row + 2, column=1, value="Clones obtained by the TRUST4 with a bam file")

### write clones by row from file
for row_data in clones_trust_bam:
    ws.append(row_data)

## TRUST FASTQ
last_row = ws.max_row

### add new header - source of clones
ws.cell(row=last_row + 2, column=1, value="Clones obtained by the TRUST4 with a fastq files")

### write clones by row from file
for row_data in clones_trust_fastq:
    ws.append(row_data)

## MiXCR
last_row = ws.max_row

### add new header - source of clones
ws.cell(row=last_row + 2, column=1, value="Clones obtained by the MiXCR")

### write clones by row from file
for row_data in clones_mixcr:
    ws.append(row_data)


# add new columns with comparing results (+ and -)
## TRUST BAM
col_num = 6
ws.insert_cols(col_num)
for i, value in enumerate(new_column_trust_bam, start=1):
    ws.cell(row=i, column=col_num, value=value)

## TRUST FASTQ
col_num = 7
ws.insert_cols(col_num)
for i, value in enumerate(new_column_trust_fastq, start=1):
    ws.cell(row=i, column=col_num, value=value)

## MiXCR
col_num = 8
ws.insert_cols(col_num)
for i, value in enumerate(new_column_mixcr, start=1):
    ws.cell(row=i, column=col_num, value=value)

wb.save(file_path)
' $res_xlsx $sample "${work_dir}/${sample}/${sample}_trust4_bam/${sample}_clones.csv" "${work_dir}/${sample}/${sample}_trust4_fastq/${sample}_clones.csv" $res_tsv_mixcr "${new_column_trust_bam[@]}" "${new_column_trust_fastq[@]}" "${new_column_mixcr[@]}"
#xlsx
}
export -f add_res_to_xlsx


# creating fasta file with bcr/tcr regions from gencode gtf
#create_suppl_files_trust4

#parallel -j $jobs run_trust4  :::: ./samples.txt

#parallel -j $jobs run_mixcr :::: ./samples.txt
#parallel -j $jobs merge_results_mixcr :::: ./samples.txt

#parallel -j 1 add_res_to_xlsx :::: ./samples.txt
