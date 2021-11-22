# Analysis Start Prompt
echo -e 'RNA-Seq Analysis Start\n'
# Create folder for log files
mkdir -p log
# 1.Assign all files needed for analysis
# 1.1.Assign raw fastq file path
fq_path='/localdisk/data/BPSM/AY21/fastq'

# 1.2.Genome file And Bed file
fa_file='/localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz'
bed_file='/localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed'
dir=$(pwd)
ls -l $fq_path | awk '{FS=" ";  {print $9;}}' | sort | grep "fq.gz" | awk '{print(substr($1,1,13))}' | sort | uniq >./sampleid.txt
echo -e "Loading fasta file from $fa_file......\n"
echo -e "Loading bed file from $bed_file......\n"
echo -e "Fastq files saved in $fq_path......\n"


# 2.Perperation for mapping
echo -e "Preparation for mapping: building index for reference genome.......\n"
mkdir -p genome
# copy fasta file to genome folder
cp $fa_file ./genome/Genome.fasta.gz
gzip -d ./genome/Genome.fasta.gz
# hisat2 index building, redirect output to index.log file in log folder
hisat2-build ./genome/Genome.fasta ./genome/Genome &> "./log/index.log"
index_path='./genome/Genome'

# 3. fastq to bam(.bai)
## Using a loop to analysis each pair of fastq data
## Using hisat2 to align reads, Using samtools to generate bam files and its index
mkdir -p map
mkdir -p fastqc
mkdir -p ./log/fastqc
echo -e "Step 1: Quality Contro AND Reads Alignment\n"
cat ${PWD}/sampleid.txt | while read id
do
    # 3.1 Perform quality control
    ## -t thread: 8
    ## output to fastqc folder
    echo -e "Analysing quality of sample ${id:5:8}\n"
    fastqc --extract ${fq_path}/${id}_1.fq.gz ${fq_path}/${id}_2.fq.gz -t 8 -o ${PWD}/fastqc/ &> "./log/fastqc/${id:5:8}.log"
    cat ${PWD}/fastqc/${id}_1_fastqc/summary.txt >> fastqc_summary.txt
    echo -e "\n" >> fastqc_summary.txt 
    cat ${PWD}/fastqc/${id}_2_fastqc/summary.txt >> fastqc_summary.txt
    echo -e "\n" >> fastqc_summary.txt 
    # 3.2 Perform reads alignment
    ## -p thread: 8
    ## sam file output to map folder
    ## log message redirects to map folder and save with filename mapeval(mapping-evaluation)
    echo -e "Aligning reads of sample ${id:5:8}\n"
    hisat2 -p 8 -x $index_path -1 ${fq_path}/${id}_1.fq.gz -2 ${fq_path}/${id}_2.fq.gz -S ./map/${id:5:8}.sam &> "./map/${id:5:8}.mapeval" 
    # 3.3 sam to bam
    echo -e "Converting the output to indexed "bam" format of sample ${id:5:8}\n"
    samtools sort -@ 8 ./map/${id:5:8}.sam -o ./map/${id:5:8}.bam
    # 3.4 bam index
    samtools index ./map/${id:5:8}.bam ./map/${id:5:8}.bam.bai
done

# 4,5. bam to count AND mean expression level
# Using a loop to analysis each pair of fastq data
# Using hisat2 to align reads, Using samtools to generate bam files and its index

# create a folder to save files that containing counts of samples
# counts of samples within a group are save in one file, the file name is the condition of each group  
mkdir -p group_count 
## create a folder to save mean count of each group
mkdir -p group_mean
cd ./map
## copy sample datail file
cp /localdisk/data/BPSM/AY21/fastq/100k.fqfiles .

## creating grouping file: the first column are expriment condition
## the second column lists the position of bam file regagrding to each group
echo -e "Creating file of grouping information......\n"
sed "1d" 100k.fqfiles | cut -f2,4,5 | sed 's/\t/_/g' > condition.txt
sed "1d" 100k.fqfiles | cut -f1 | paste condition.txt - | \
awk  '{a[$1]=a[$1]?a[$1]","$2:$2}END{for (i in a) print i,a[i]}' > grouping_file.txt

## Generation count file for each group and Calculation of mean count
for i in $(seq 1 $(cat grouping_file.txt | wc -l))
do
    ## extract sample name
    s=$(sed -n "${i}p" grouping_file.txt | awk '{print $2}')
    s=${s//,/.bam } # name of corresponding bam file
    ## generate count
    echo -e "Generating count file for condition $(sed -n "${i}p" grouping_file.txt | awk '{print $1}')\n"
    echo -e "$s sample(s) are conducted in the condition $(sed -n "${i}p" grouping_file.txt | awk '{print $1}')\n"
    echo -e "bedtools multicov generating......\n"
    bedtools multicov -bams ${s//,/ }.bam -bed $bed_file > $dir/group_count/$(sed -n "${i}p" grouping_file.txt | awk '{print $1}').txt
    ## adding header
    s=$(sed -n "${i}p" grouping_file.txt | awk '{print $2}')
    sed -i "1 s/^/$(echo -e "id\tstart\tend\tuni\tdescription\t${s//,/\t}")\n/" $dir/group_count/$(sed -n "${i}p" grouping_file.txt | awk '{print $1}').txt
    echo -e "Count file for $(sed -n "${i}p" grouping_file.txt | awk '{print $1}') have been saved.\n"
    ## calculate mean line by line
    echo -e "Calculating mean count for group $(sed -n "${i}p" grouping_file.txt | awk '{print $1}')\n"
    cat $dir/group_count/$(sed -n "${i}p" grouping_file.txt | awk '{print $1}').txt | cut -f 6- |\
    awk '{FS="\t"; sum=0;for(i=1;i<=NF;i++) sum+=$i; print sum/=NF}' > $dir/group_mean/$(sed -n "${i}p" grouping_file.txt | awk '{print $1}')_mean.txt 
done
cd ..


# 6. Group-wise comparisons
mkdir -p foldchange
## Creating group comparisons file
for ((i=1;i<12;i=i+3))
do  
    echo -e $(cat ./map/grouping_file.txt | awk '{IFS='\t';print $1}' | sed 's/_/\t/g' | sort -k2,2 -k3,3 | sed 's/\t/_/g' | sed -n ${i}p)'\t'\
    $(cat ./map/grouping_file.txt | awk '{IFS='\t';print $1}' | sed 's/_/\t/g' | sort -k2,2 -k3,3 | sed 's/\t/_/g' | sed -n $(( $i + 1 ))p) >> foldchange_group.txt;
    echo -e $(cat ./map/grouping_file.txt | awk '{IFS='\t';print $1}' | sed 's/_/\t/g' | sort -k2,2 -k3,3 | sed 's/\t/_/g' | sed -n ${i}p) '\t'\
    $(cat ./map/grouping_file.txt | awk '{IFS='\t';print $1}' | sed 's/_/\t/g' | sort -k2,2 -k3,3 | sed 's/\t/_/g' | sed -n $(( $i + 2 ))p)'\t'>> foldchange_group.txt;
    echo -e $(cat ./map/grouping_file.txt | awk '{IFS='\t';print $1}' | sed 's/_/\t/g' | sort -k2,2 -k3,3 | sed 's/\t/_/g' | sed -n $(( $i + 1 ))p)'\t'\
    $(cat ./map/grouping_file.txt | awk '{IFS='\t';print $1}' | sed 's/_/\t/g' | sort -k2,2 -k3,3 | sed 's/\t/_/g' | sed -n $(( $i + 2 ))p) >> foldchange_group.txt;
done
col=$(ls ./group_count/ | sed -n '1p')
cat ./group_count/$col | cut -f 1-5 > descrpition.txt
sed -i "1d" descrpition.txt 

## read group comparisons file line by line, extract two groups' count and calculate fold-changes
while read groupA groupB; do
    # combine groupA and groupB counts into a same file
    paste ./group_mean/${groupA}_mean.txt ./group_mean/${groupB}_mean.txt > ./foldchange/${groupA}_vs_${groupB}.txt
    # calculate fold-changes(groupA_count/groupB_count)
    cat ./foldchange/${groupA}_vs_${groupB}.txt | \
    awk '{if ($2==0) print "NA"; else print $1/$2}' > ./foldchange/${groupA}_vs_${groupB}_foldchange.txt
    paste descrpition.txt ./foldchange/${groupA}_vs_${groupB}_foldchange.txt  > ./foldchange/${groupA}_vs_${groupB}.txt
    # sort file by fold-changes
    cat ./foldchange/${groupA}_vs_${groupB}.txt | sort -t$'\t' -k6,6nr > ./foldchange/${groupA}_vs_${groupB}_foldchange.txt  
done < foldchange_group.txt

## log2 fold change
mkdir -p log2foldchange
while read groupA groupB; do
    # combine groupA and groupB counts into a same file
    paste ./group_mean/${groupA}_mean.txt ./group_mean/${groupB}_mean.txt > ./log2foldchange/${groupA}_vs_${groupB}.txt
    # calculate fold-changes(groupA_count/groupB_count)
    cat ./log2foldchange/${groupA}_vs_${groupB}.txt | \
    awk '{for(i=1;i<=NF;i++){$i+=1}}{print $0}' > ./log2foldchange/${groupA}_vs_${groupB}_log2foldchange.txt
    cat ./log2foldchange/${groupA}_vs_${groupB}_log2foldchange.txt | \
    awk '{for(i=1;i<=NF;i++){$i=log($i)/log(2)}}{print $0}' > ./log2foldchange/${groupA}_vs_${groupB}.txt
    cat ./log2foldchange/${groupA}_vs_${groupB}.txt | \
    awk '{print $1-$2}' > ./log2foldchange/${groupA}_vs_${groupB}_log2foldchange.txt
    paste descrpition.txt ./log2foldchange/${groupA}_vs_${groupB}_log2foldchange.txt  > ./log2foldchange/${groupA}_vs_${groupB}_log2fc.txt
    # sort file by log2fold-changes
    cat ./log2foldchange/${groupA}_vs_${groupB}_log2fc.txt | sort -t$'\t' -k6,6nr > ./log2foldchange/${groupA}_vs_${groupB}_log2foldchange.txt 
done < foldchange_group.txt

## Pasting mean count file
paste descrpition.txt ./group_mean/*_mean.txt > mean_count.txt
ls ./group_mean/*_mean.txt > mean_name
cut -d "/" -f 3 mean_name > mn
cut -d "." -f 1 mn > mean_name
n=$(cat mean_name | xargs -n 15)
sed -i "1 s/^/$(echo -e "chr\tstart\tend\tname\tdescription\t${n// /\t}")\n/"  mean_count.txt


## removing extra files
rm -rf descrpition.txt; rm -rf sampleid.txt; rm -rf mn; rm -rf mean_name
rm -rf ./log2foldchange/*log2fc.txt
rm -rf ./log2foldchange/*Induced.txt
rm -rf ./foldchange/*Induced.txt
rm -rf ./group_mean
echo -e "RNA-Seq Analysis finish."
