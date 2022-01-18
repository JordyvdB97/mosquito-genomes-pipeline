# Pre-processing

## 1. Merging the read files to a single file

	zcat *read_file_name1*_R1.fastq.gz *read_file_name2*_R1.fastq.gz | gzip > *accession_1*_R1.fastq.gz
	zcat *read_file_name1*_R2.fastq.gz *read_file_name2*_R2.fastq.gz | gzip > *accession_1*_R2.fastq.gz

## 2. Quality assessment and trimming

Trimming the low quality ends of the sequences with _fastp_ (Chen et al. 2018).

	fastp \
                -i $IN1 \
                -I $IN2 \
                -o ${IN1}.qt \
                -O ${IN2}.qt \
                -j fastp.json -h fastp.html --verbose

# Assembly

## 3. Indexing
Indexing the forward and reverse sequence. As well as the reference genome. Indexing of the reference genome needs to be to be done only once to efficiently map sequences to it. 

	# fw reads
        IN1=/fileserver/*accession_1*.fastq
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
	        ${IN1}.qt ${IN2}.qt | samtools view \
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


