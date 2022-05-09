# Pre-processing
## 0. Software installation
https://github.com/Kinggerm/GetOrganelle/wiki/Installation https://github.com/RemiAllio/MitoFinder#get-and-install-mitofinder-linux
	
	sudo apt install fastp ncbi-entrez-direct spades bowtie2 ncbi-blast+ pigz automake autoconf gcc default-jre
	conda install -c bioconda getorganelle
	
	conda install -c bioconda seqkit
	

## 1. Merging the read files to a single file
To about 40% of the sequencing samples have have more than one readfiles, we first need to merge these files together. We use a loop that counts the number of files per sample. If that equals two (biderectional read files) it copies and renames the files to the folder merged. If we have more files, unzips the files per direction and zip it into a new merged container. 

		
	mkdir /fileserver/1merged/
	
	SAMPLES=/fileserver/sample_no.txt

	for i in $(cat $SAMPLES)
	   do
		if [ "$(ls /fileserver/0raw/*_"$i"_*.gz | wc -l)" -eq 2 ]
	 	then
	    		cp /fileserver/0raw/*_"$i"_*_R1*.gz /fileserver/1merged/"$i"_R1.fastq.gz
			cp /fileserver/0raw/*_"$i"_*_R2*.gz /fileserver/1merged/"$i"_R2.fastq.gz
		else
	    		zcat /fileserver/0raw/*_"$i"_*_R1*.gz | gzip > /fileserver/1merged/"$i"_R1.fastq.gz
	    		zcat /fileserver/0raw/*_"$i"_*_R2*.gz | gzip > /fileserver/1merged/"$i"_R2.fastq.gz
	 	fi
	   done

## 2. Quality assessment and trimming

Trimming the low quality ends of the sequences with _fastp_ (Chen et al. 2018).

	mkdir /fileserver/2fastp_trimmed/
	mkdir /fileserver/2fastp_trimmed/fastp_reports
	
	for i in $(cat $SAMPLES)
	   do
		fastp \
			-i /fileserver/1merged/"$i"_R1.fastq.gz \
			-I /fileserver/1merged/"$i"_R2.fastq.gz \
			-o /fileserver/2fastp_trimmed/"$i"_trimmed_R1.fastq.gz \
			-O /fileserver/2fastp_trimmed/"$i"_trimmed_R2.fastq.gz \
			-j /fileserver/2fastp_trimmed/fastp_reports/fastp_"$i".json \
			-h /fileserver/2fastp_trimmed/fastp_reports/fastp_"$i".html --verbose
 	   done

# De Novo Assembly

## 3. Normalizing
Normalizing the coverage by down-sampling the reads in high-depth areas of a genome (to an average coverage of 100x), and removing reads with a low coverage (less than 5x) with the BBNorm tool of BBMap. This accelarates the de novo assembly, reduces the memory needed for the assembly smaller, and the dataset more tractable for the assembler. Also it seems to improve the assembly quality (@@@ check the original reference @@@). BBMap automatically recognizes that the two imputfiles contain readpairs

	mkdir /fileserver/3bbmap_normalized/

	for i in $(cat $SAMPLES)
	   do
		/directory/bbnorm.sh \
			in=/fileserver/2fastp_trimmed/"$i"_trimmed_R1.fastq.gz \
			in2=/fileserver/2fastp_trimmed/"$i"_trimmed_R2.fastq.gz \
			out=/fileserver/3bbmap_normalized/"$i"_normalized.fq target=100 min=5
	   done
## 3. Downloading reference library
First we construct a local reference database. This have to be only done once. We download the  

	mkdir /mnt/e/2020_mtmozseq/blastdb

	taxa=Culicidae
	esearch -db nuccore -query "\"mitochondrion\"[All Fields] \
	AND (\"${taxa}\"[Organism]) AND (refseq[filter] \
	AND mitochondrion[filter] AND (\"12000\"[SLEN] : \"20000\"[SLEN]))" | \
	tee >(efetch -format gbwithparts > /mnt/e/2020_mtmozseq/blastdb/culicidae_mt_refseq.gb) | \
	efetch -format fasta > /mnt/e/2020_mtmozseq/blastdb/culicidae_mt_refseq.fasta
	
	makeblastdb -in /mnt/e/2020_mtmozseq/blastdb/culicidae_mt_refseq.fasta -parse_seqids -dbtype nucl

## 4. De Novo assembly

	mkdir /mnt/e/2020_mtmozseq/3getorganelle_assembly/
	
	for i in $(cat $SAMPLES)
	   do
		get_organelle_from_reads.py \
			-s /mnt/e/2020_mtmozseq/blastdb/culicidae_mt_refseq.fasta \
			-1 /mnt/e/2020_mtmozseq/2fastp_trimmed/"$i"_trimmed_R1.fastq.gz \
			-2 /mnt/e/2020_mtmozseq/2fastp_trimmed/"$i"_trimmed_R2.fastq.gz \
			-o /mnt/e/2020_mtmozseq/4getorganelle_assembly/"$i"_getorganelle \
			-R 5 -k 21,45,65,85,105 -F animal_mt
	   done

## 4. nBLAST results

	mkdir /mnt/e/2020_mtmozseq/5nblast_results
	mkdir /mnt/e/2020_mtmozseq/6mtgenomes
	
copy contigs from the assembly to new folder

	cp /fileserver/4getorganelle_assembly/*accession_number*_getorganelle/extended_spades/contigs.fasta \ 
	/fileserver/5nblast_results/*accession_number*_all_contigs.fasta
	
selecting only contigs larger than 12000 bps

	awk -v n=12000 '/^>/{ if(l>n) print b; b=$0;l=0;next } {l+=length;b=b ORS $0}END{if(l>n) print b }' \ 
	/fileserver/5nblast_results/*accession_number*_all_contigs.fasta > \
	/fileserver/5nblast_results/*accession_number*_contigs_l12000.fasta
	
blastn contigs against local database of mosquito mt sequences, both into a detailed table and a list with only contig names
	
	blastn \
		-query /fileserver/5nblast_results/*accession_number*_contigs_l12000.fasta \
		-db /fileserver/blastdb/culicidae_mt_refseq.fasta \
		-outfmt "10 qseqid sseqid ssciname pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
		-evalue 1e-30 \
		-out /fileserver/5nblast_results/*accession_number*_blast_results.csv
	blastn \
		-query /fileserver/5nblast_results/*accession_number*_contigs_l12000.fasta 
		-db /fileserver/blastdb/culicidae_mt_refseq.fasta 
		-outfmt "6 qseqid"  
		-evalue 1e-30 
		-max_target_seqs 1 | sort | uniq > /fileserver/5nblast_results/*accession_number*_blast_results.lst
		
extract contigs with a match in the database
	
	seqtk subseq \
		/fileserver/5nblast_results/*accession_number*_contigs_l12000.fasta \
		/fileserver/5nblast_results/*accession_number*_blast_results.lst > \
		/fileserver/6mtgenomes/*accession_number*_blasted_contigids.fasta

correct overlapping ends of with python script
	
	python /directory/simple_circularise.py \
		/fileserver/6mtgenomes/*accession_number*_blasted_contigids.fasta \
		/fileserver/6mtgenomes/*accession_number*_blasted_contigids_circularized.fasta

# Mapping against reference

## 1. Paired-end mapping assembly

	minimap2 \
		-ax sr \
		/fileserver/reference_sequence.fasta \
		/fileserver/2fastp_trimmed/*accession_number*_trimmed_R1.fastq.gz \
		/fileserver/2fastp_trimmed/*accession_number*_trimmed_R2.fastq.gz | \
	samtools view \
		-b -u -F 0x04 --threads 4 \
		-o /fileserver/4minimap/*accession_number*.bam
		
## 2. cleaning mapping results

Sort by read name

	samtools sort \
		-n -l 0 -m 3G --threads 4 \
		-o /fileserver/4minimap/*accession_number*.nsort.bam \
		/fileserver/4minimap/*accession_number*.bam

Correct mate pairs

	samtools fixmate \
		-r -m --threads 4 \
		/fileserver/4minimap/*accession_number*.nsort.bam \
		/fileserver/4minimap/*accession_number*.fixmate.bam

Sort numerically

	samtools sort \
		-l 0 -m 3G --threads 4 \
		-o /fileserver/4minimap/*accession_number*.sort.bam \
		/fileserver/4minimap/*accession_number*.fixmate.bam

Remove duplicates

	samtools markdup \
		-r --threads 4 \
		/fileserver/4minimap/*accession_number*.sort.bam \
		/fileserver/4minimap/*accession_number*.markdup.bam

## 3. variant calling

	bcftools mpileup \
		-Ou -f \
		/fileserver/reference_sequence.fasta \
		/fileserver/4minimap/*accession_number*.markdup.bam | \
	bcftools call \
		-mv -Oz \
		-o /fileserver/4minimap/*accession_number*.calls.vcf.gz

Indexing the vcf.gz file

	bcftools index \
		/fileserver/4minimap/*accession_number*.calls.vcf.gz

## 4. create consensus sequence
cat /media/jordy/mosqdisk/2020_mtmozseq/final_genomes/done/RMNH.INS.1271376_Ae_albopictus.fasta | bcftools consensus /media/jordy/mosqdisk/2020_mtmozseq/4minimap/103943-047-079.calls.vcf.gz > /media/jordy/mosqdisk/2020_mtmozseq/4minimap/103943-047-079.consensus.fasta

reformat.sh in=/media/jordy/mosqdisk/2020_mtmozseq/4minimap/103943-047-079.consensus.fasta out=/media/jordy/mosqdisk/2020_mtmozseq/final_genomes/done/consensus_sequences/103943-047-079.consensus.fasta fastawrap=0 tuc

## Annotation

SAMPLES=/media/jordy/mosqdisk/2020_mtmozseq/sample_name.txt

for i in $(cat $SAMPLES)
do 
mitofinder -j ${i} -a /media/jordy/mosqdisk/2020_mtmozseq/7final_fasta_genomes/${i}.fasta -t "mitfi" -r /media/jordy/mosqdisk/2020_mtmozseq/blastdb/culicidae_mt_refseq.gb -o 5
done

SAMPLES=/media/jordy/mosqdisk/2020_mtmozseq/sample_name.txt

for i in $(cat $SAMPLES)
do 
awk -v seq="$i@COX1" -v RS='>' '$1 == seq {print RS $0}' /media/jordy/mosqdisk/2020_mtmozseq/8annotation/"$i"/"$i"_MitoFinder_mitfi_Final_Results/"$i"_final_genes_NT.fasta > /media/jordy/mosqdisk/2020_mtmozseq/cox1_sequences/"$i"_COX1.fasta
done
