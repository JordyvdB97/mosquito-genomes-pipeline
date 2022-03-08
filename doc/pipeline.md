# Pre-processing
## 0. Software installation
https://github.com/Kinggerm/GetOrganelle/wiki/Installation
	
	sudo apt install fastp ncbi-entrez-direct spades bowtie2 ncbi-blast+ pigz
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

# Assembly

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


## 3. Indexing
Indexing the forward and reverse sequence. As well as the reference genome. Indexing of the reference genome needs to be to be done only once to efficiently map sequences to it. 

	# fw reads
        IN1=./merged/"$i"_trimmed_R1.fastq.gz
	# fw reads
        IN2=/fileserver/*accession_2*.fastq

	# cpDNA fasta refseq, indexed with minimap2
	# https://github.com/naturalis/tomatogenome-en-tibi/blob/master/doc/Pipeline.md#3-indexing-the-reference
        REF=/fileserver/Solanum_lycopersicum_NC_007898.fasta

	# base name of the output file
        BASE=/fileserver/*accession*	

## 4. Paired-end mapping assembly
We use _minimap2_ for mapping the reads against our reference genome. Output should be SAM-format (-a) on four cores (-t 4). Than use _samtools view_ to convert to a uncompressed (-u) BAM-file (-b), excluding unmapped reads (-F 0x04).

	minimap2 -ax sr -a -t 4 $REF \
	        ${IN1} ${IN2} | samtools view \
	       -b -u -F 0x04 --threads 4 -o ${BASE}.bam -

## 5. Sort by read name
We sort the reads by name using _samtools sort_ as most downstraem correction of the matepairs: https://www.biostars.org/p/365882/#365887.

    samtools sort -n -l 0 -m 3G --threads 4 \
	    -o ${BASE}.nsort.bam ${BASE}.bam

## 6. Correct mate pairs
A proofread step. We need to correct for any flaws in read-pairing that may have been introduced by the minimap2. It ensures that imporatnt fields are correct and internally consistent.

    samtools fixmate -r -m  --threads 4 \
    	    ${BASE}.nsort.bam ${BASE}.fixmate.bam

## 7. Sort numerically
We sort now nummerically, i.e. location of the read on the genome. 

    samtools sort -l 0 -m 3G --threads 4 \
	    -o ${BASE}.sort.bam ${BASE}.fixmate.bam

## 8. Remove duplicates

	samtools markdup -r --threads 4 \
	        ${BASE}.sort.bam ${BASE}.markdup.bam

## 9. Create consensus sequence for each individual

       bcftools mpileup -Ou -f Solanum_lycopersicum_NC_007898.fa -o *accession_number*.pileup.bcf *accession_number*.markdup.bam

## 10. Convert BCF to FASTA

       cd /home/vasia.kakakiou/samtools/bin/PGD/

       for file in /home/vasia.kakakiou/samtools/bin/PGD/*.bcf
       do
       java -Xmx5g -jar /home/vasia.kakakiou/samtools/bin/PGD/PGDSpider2-cli.jar -inputfile "$file" -outputfile "${file%.bcf}.fasta" -inputformat BCF -outputformat FASTA -spid BCF_to_FASTA.spid
       done

# Gene splicing and SNP calling

## 12. Allignment and SNP calling
#Combination of make mpileup and calls
#mpileup: -Ou flags for specifying output type as uncompressed BCF, -f to specify reference file
#call: -m multiallelic-caller, -v to output variant sites only , -Oz output type as compressed VCF, -o specifying output file name

	bcftools mpileup -Ou -f /location_of_reference/reference.fasta /location_of_input/sorted_mapped_[SAMPLE].bam | bcftools call -mv -Oz -o called_[SAMPLE].vcf.gz

#Indexing VCF files

	bcftools index called_[SAMPLE].vcf.gz

## 11. Gene splicing
Genen uit de consensus sequences splicen. 

Dacht deze code:

	samtools faidx reference.fasta [ref_fa_header]:[START coord – STOP coord] | bcftools consensus called_[SAMPLE].vcf.gz > [GENE]_[SAMPLE].fasta


