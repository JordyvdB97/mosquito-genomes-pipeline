# Mitogenome UNIX pipeline
This R script accompanies the unpublished manuscript:  

Van der Beek, J. G., Ibáñez-Justicia, A., Biesmeijer, J. C., Lizarazo-Forero, E., Stroo, A., van de Vossenberg, B. T. L. H., Warbroek, T., & Schrama, M. J. J. (n.d.). The differentiating power of mitochondrial genes: complete mitogenome sequences of 27 mosquito species present in Europe. [Unpublished manuscript].

- 📁 **All raw sequencing data can be downloaded [here](https://www.ncbi.nlm.nih.gov/sra/PRJNA1219649)** (pre-merged read files)
- 📁 **Assembled genomes are available here**

For questions or further inquiries, please contact: Jordy van der Beek (jordy.vanderbeek@naturalis.nl).

# Pre-processing
## 0. Software installation
To begin the analysis, you need to install the necessary software packages. Below are the commands and links for installation instructions:

- **GetOrganelle**: Installation guide available [here](https://github.com/Kinggerm/GetOrganelle/wiki/Installation)
- **MitoFinder**: Installation instructions available [here](https://github.com/RemiAllio/MitoFinder#get-and-install-mitofinder-linux)
- **Seqkit**: Installation instructions available [here](https://bioinf.shenwei.me/seqkit/)

Run the following commands to install the other essential software: 
	
	sudo apt install fastp ncbi-entrez-direct spades bowtie2 ncbi-blast+ pigz automake autoconf gcc default-jre
	

## 1. Merging the read files to a single file
Some sequencing samples may have multiple read files per sample. This step merges these files to ensure that subsequent processing steps operate on consolidated data.

Creating directory + defining file:
- <code>mkdir /fileserver/1merged/</code>: Create a directory for merged read files.
- <code>SAMPLES=/fileserver/sample_no.txt</code>: Define a file containing sample identifiers.

Merging process:		

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

This loop checks whether a sample has two read files. If true, it copies and renames them. If there are more, it decompresses, merges, and compresses them into a new file.

## 2. Quality assessment and trimming

To enhance data quality, remove low-quality bases from sequencing reads using [fastp](https://github.com/OpenGene/fastp).

Creating directories:
- <code>mkdir /fileserver/2fastp_trimmed/</code>: Directory for trimmed reads.
- <code>mkdir /fileserver/2fastp_trimmed/fastp_reports</code>: Directory for reports.

Trimming process:
	
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

fastp generates JSON and HTML reports and creates trimmed FASTQ files for both read directions.

# De Novo Assembly

## 3. Normalizing
Normalize read coverage with _BBNorm_ to reduce data complexity and improve assembly efficiency. This accelarates the de novo assembly, reduces the memory needed for the assembly smaller, and the dataset more tractable for the assembler. 

Creating directory:
- <code>mkdir /fileserver/3bbmap_normalized/</code>: Directory for normalized reads.

Normalization process:

	mkdir /fileserver/3bbmap_normalized/

	for i in $(cat $SAMPLES)
	   do
		/directory/bbnorm.sh \
			in=/fileserver/2fastp_trimmed/"$i"_trimmed_R1.fastq.gz \
			in2=/fileserver/2fastp_trimmed/"$i"_trimmed_R2.fastq.gz \
			out=/fileserver/3bbmap_normalized/"$i"_normalized.fq target=100 min=5
	   done

The command targets an average coverage depth of 100x and filters out regions with less than 5x coverage.

## 3. Downloading reference library
Construct a local reference database from mitochondrial genome sequences. 

Making directory:
- <code>mkdir /mnt/e/2020_mtmozseq/blastdb</code>: Directory for the database.

Downloading process:

	taxa=Culicidae
	esearch -db nuccore -query "\"mitochondrion\"[All Fields] \
	AND (\"${taxa}\"[Organism]) AND (refseq[filter] \
	AND mitochondrion[filter] AND (\"12000\"[SLEN] : \"20000\"[SLEN]))" | \
	tee >(efetch -format gbwithparts > /mnt/e/2020_mtmozseq/blastdb/culicidae_mt_refseq.gb) | \
	efetch -format fasta > /fileserver/blastdb/culicidae_mt_refseq.fasta
	makeblastdb -in /fileserver/blastdb/culicidae_mt_refseq.fasta -parse_seqids -dbtype nucl

As the Refseq dataset is continously updated we placed the uploaded a copy of the database as it was in January 2022 (139 mitogenomes - 139 species; 17 genera) in this Github repoistory: [link](https://github.com/JordyvdB97/mosquito-genomes-pipeline/blob/5eb2f5d19e5e043b5f727edfb5274728bbfae4c5/Culicidae_mt_GenBank_refseq_January_2022.gb).

## 4. De Novo assembly
Assemble mitochondrial genomes using [GetOrganelle](https://github.com/Kinggerm/GetOrganelle/).

Creating directory:
- <code>mkdir /fileserver/4getorganelle_assembly/</code>: Directory for assembly results.

Assembly process:
	
	for i in $(cat $SAMPLES)
	   do
		get_organelle_from_reads.py \
			-s /fileserver/blastdb/culicidae_mt_refseq.fasta \
			-1 /fileserver/2fastp_trimmed/"$i"_trimmed_R1.fastq.gz \
			-2 /fileserver/2fastp_trimmed/"$i"_trimmed_R2.fastq.gz \
			-o /fileserver/4getorganelle_assembly/"$i"_getorganelle \
			-R 5 -k 21,45,65,85,105 -F animal_mt
	   done

This command assembles mitochondrial sequences by specifying k-mer sizes and output directories.

## 4. nBLAST results

Creating directories:
- <code>mkdir /fileserver/5nblast_results</code>
- <code>mkdir /fileserver/6mtgenomes</code>
	
### Copy contigs from the assembly to new folder

	cp /fileserver/4getorganelle_assembly/*accession_number*_getorganelle/extended_spades/contigs.fasta \ 
	/fileserver/5nblast_results/*accession_number*_all_contigs.fasta
	
### Selecting Contigs

Extract contigs larger than 12,000 base pairs and filter matches from the reference database.

	awk -v n=12000 '/^>/{ if(l>n) print b; b=$0;l=0;next } {l+=length;b=b ORS $0}END{if(l>n) print b }' \ 
	/fileserver/5nblast_results/*accession_number*_all_contigs.fasta > \
	/fileserver/5nblast_results/*accession_number*_contigs_l12000.fasta

### BLAST Analysis
Perform local BLAST against the reference database and extract matching sequences.

Creating a detailed table:

	blastn \
		-query /fileserver/5nblast_results/*accession_number*_contigs_l12000.fasta \
		-db /fileserver/blastdb/culicidae_mt_refseq.fasta \
		-outfmt "10 qseqid sseqid ssciname pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
		-evalue 1e-30 \
		-out /fileserver/5nblast_results/*accession_number*_blast_results.csv

Creating a list with only contig names:
  
	blastn \
		-query /fileserver/5nblast_results/*accession_number*_contigs_l12000.fasta 
		-db /fileserver/blastdb/culicidae_mt_refseq.fasta 
		-outfmt "6 qseqid"  
		-evalue 1e-30 
		-max_target_seqs 1 | sort | uniq > /fileserver/5nblast_results/*accession_number*_blast_results.lst
		
Extract relevant contigs:

	seqtk subseq \
		/fileserver/5nblast_results/*accession_number*_contigs_l12000.fasta \
		/fileserver/5nblast_results/*accession_number*_blast_results.lst > \
		/fileserver/6mtgenomes/*accession_number*_blasted_contigids.fasta

### Correct overlapping ends
Correct overlapping ends with the python script [SympleCircularise](https://github.com/Kzra/Simple-Circularise/tree/master):
	
	python /directory/simple_circularise.py \
		/fileserver/6mtgenomes/*accession_number*_blasted_contigids.fasta \
		/fileserver/6mtgenomes/*accession_number*_blasted_contigids_circularized.fasta

# Mapping against reference

## 1. Paired-end mapping assembly
Creating directories:

- <code>mkdir /fileserver/4minimap</code>
- <code>mkdir /fileserver/7final_fasta_genomes</code>

Align paired-end reads to a reference sequence using [minimap2](https://github.com/lh3/minimap2):

	minimap2 \
		-ax sr \
		/fileserver/reference_sequence.fasta \
		/fileserver/2fastp_trimmed/*accession_number*_trimmed_R1.fastq.gz \
		/fileserver/2fastp_trimmed/*accession_number*_trimmed_R2.fastq.gz | \
	samtools view \
		-b -u -F 0x04 --threads 4 \
		-o /fileserver/4minimap/*accession_number*.bam
		
## 2. cleaning mapping results
Sort, correct mate pairs, and remove duplicates using [samtools](https://www.htslib.org/):

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
Call variants using [bcftools](https://samtools.github.io/bcftools/):

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
Generate a consensus sequence:

	cat \
		/fileserver/reference_sequence.fasta | \
	bcftools consensus \
		/fileserver/4minimap/*accession_number*.calls.vcf.gz > \
		/fileserver/4minimap/*accession_number*.consensus.fasta

	reformat.sh \
		in=/fileserver/4minimap/*accession_number*.consensus.fasta \
		out=/fileserver/7final_fasta_genomes/*accession_number*.consensus.fasta \
		fastawrap=0 tuc

# Annotation
Annotate mitochondrial genomes using [MitoFinder](https://github.com/RemiAllio/MitoFinder):

	mitofinder \
		-j *sample_name* \
		-a /fileserver/7final_fasta_genomes/*accession_number*.consensus.fasta \
		-t "mitfi" \
		-r /media/jordy/mosqdisk/2020_mtmozseq/blastdb/culicidae_mt_refseq.gb \
		-o 5
