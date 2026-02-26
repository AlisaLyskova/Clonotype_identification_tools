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

function vidjil_run()
{
	sample=$1
    R1=${work_dir}/${sample}_R1.fastq.gz
    R2=${work_dir}/${sample}_R2.fastq.gz
	merged_reads=${work_dir}/${sample}_reads_merged.fastq.gz
	seqtk mergepe $R1 $R2 | pigz -c -p 8 > $merged_reads
	mkdir "${work_dir}/vidjil_res"
	$vidjil -g "${vidjil_dir}/germline/homo-sapiens.g" -o "${work_dir}/vidjil_res" $merged_reads > "${work_dir}/vidjil_res/vidjil.log" 2>&1
}
export -f vidjil_run

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
	res_tsv_vidjil="${work_dir}/${sample}/vidjil_res/${sample}_reads_merged.fastq.tsv"

	# creating file with clones only
	cat $res_tsv_trust_bam | sed '1d' | awk -v OFS="\t" -F "\t" '{print $1, $5,$6,$7}' > "${work_dir}/${sample}/${sample}_trust4_bam/${sample}_clones.csv"
	cat $res_tsv_trust_fastq | sed '1d' | awk -v OFS="\t" -F "\t" '{print $1, $5,$6,$7}' > "${work_dir}/${sample}/${sample}_trust4_fastq/${sample}_clones.csv"
	cat $res_tsv_vidjil | awk -v OFS="\t" -F "\t" '{print $2,$3,$4,$5}' | awk -F "\t" 'NR > 1 {count = 0; if ($2 != "") count++; if ($3 != "") count++; if ($4 != "") count++; if (count >= 2) print $0}' | sort -r -n -k 1 > "${work_dir}/${sample}/vidjil_res/${sample}_clones.tsv"

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
	new_column_vidjil=(Vidjil)
	for row in {1..5}
	do
		expr_vdj=$(cat $pipeline_res | awk -F "\t" '{print $5,$6,$7}' | awk '{for(i=1;i<=3;i++) sub(/\*.*/, "", $i); print}' | sed '1d' | awk -F " " -v OFS=".*" '{print $1,$2,$3}' | head -n $row | tail -n 1)
		check_1=$(cat $res_tsv_trust_bam | grep -E "$expr_vdj")
		check_2=$(cat $res_tsv_trust_fastq | grep -E "$expr_vdj")
		check_3=$(cat $res_tsv_mixcr | grep -E "$expr_vdj")
		check_4=$(cat $res_tsv_vidjil | grep -E "$expr_vdj")
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
		if [[ ! -z $check_4 ]]; then
			new_column_vidjil+=("+")
		else
			new_column_vidjil+=("-")
		fi
	done

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
clones_file_vidjil = sys.argv[6]

new_column_trust_bam = sys.argv[7:13]
new_column_trust_fastq = sys.argv[13:19]
new_column_mixcr = sys.argv[19:25]
new_column_vidjil = sys.argv[25:31]

clones_df = pd.read_csv(clones_file_trust_bam, sep="\t", header=None)
clones_trust_bam = clones_df.values.tolist()

clones_df = pd.read_csv(clones_file_trust_fastq, sep="\t", header=None)
clones_trust_fastq = clones_df.values.tolist()

clones_df = pd.read_csv(clones_file_mixcr, sep="\t").fillna(".")
clones_df = clones_df.iloc[:, :-1]
clones_mixcr = clones_df.values.tolist()

clones_df = pd.read_csv(clones_file_vidjil, sep="\t").fillna(".")
clones_vidjil = clones_df.values.tolist()

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
    
## Vidjil
last_row = ws.max_row

### add new header - source of clones
ws.cell(row=last_row + 2, column=1, value="Clones obtained by the Vidjil")

### write clones by row from file
for row_data in clones_vidjil:
    ws.append(row_data)



# add new column with + or - of clonotype summary from IGV
new_col = ["True clonotypes"]
for row in ws.iter_rows(min_row=2, max_row=6):
    igv_check = row[4].value
    count_plus = igv_check.strip().count("+")
    if count_plus >= 2:
        new_col += "+"
    else:
        new_col += "-"
col_num = 6
ws.insert_cols(col_num)
for i, value in enumerate(new_col, start=1):
    ws.cell(row=i, column=col_num, value=value)
        

# add new columns with comparing results (+ and -)
## TRUST BAM
col_num = 7
ws.insert_cols(col_num)
for i, value in enumerate(new_column_trust_bam, start=1):
    ws.cell(row=i, column=col_num, value=value)

## TRUST FASTQ
col_num = 8
ws.insert_cols(col_num)
for i, value in enumerate(new_column_trust_fastq, start=1):
    ws.cell(row=i, column=col_num, value=value)

## MiXCR
col_num = 9
ws.insert_cols(col_num)
for i, value in enumerate(new_column_mixcr, start=1):
    ws.cell(row=i, column=col_num, value=value)
    
## Vidjil
col_num = 10
ws.insert_cols(col_num)
for i, value in enumerate(new_column_vidjil, start=1):
    ws.cell(row=i, column=col_num, value=value)
    

## MiXCR and TRUST4 and Vidjil
bool_list1 = [x == "+" for x in new_column_mixcr[1:]]
bool_list2 = [x == "+" for x in new_column_trust_bam[1:]]
bool_list3 = [x == "+" for x in new_column_vidjil[1:]]
col_num = 11
ws.insert_cols(col_num)
conjunction = ["+" if a and b and c else "-" for a, b, c in zip(bool_list1, bool_list2, bool_list3)]
new_col = ["MiXCR*TRUST4*Vidjil"]+conjunction
for i, value in enumerate(new_col, start=1):
    ws.cell(row=i, column=col_num, value=value)
    
## MiXCR or TRUST4 or Vidjil
col_num = 12
ws.insert_cols(col_num)
disjunction = ["+" if a or b or c else "-" for a, b, c in zip(bool_list1, bool_list2, bool_list3)]
new_col = ["MiXCR+TRUST4+Vidjil"]+disjunction
for i, value in enumerate(new_col, start=1):
    ws.cell(row=i, column=col_num, value=value)
    
wb.save(file_path)
' $res_xlsx $sample "${work_dir}/${sample}/${sample}_trust4_bam/${sample}_clones.csv" "${work_dir}/${sample}/${sample}_trust4_fastq/${sample}_clones.csv" "${work_dir}/${sample}/vidjil_res/${sample}_clones.tsv" $res_tsv_mixcr "${new_column_trust_bam[@]}" "${new_column_trust_fastq[@]}" "${new_column_mixcr[@]}" "${new_column_vidjil[@]}"

}
export -f add_res_to_xlsx


# creating fasta file with bcr/tcr regions from gencode gtf
#create_suppl_files_trust4

#parallel -j $jobs run_trust4  :::: ./samples.txt

#parallel -j $jobs run_mixcr :::: ./samples.txt
#parallel -j $jobs merge_results_mixcr :::: ./samples.txt

#parallel -j $jobs vidjil_run :::: ./samples.txt

#parallel -j 1 add_res_to_xlsx :::: ./samples.txt
