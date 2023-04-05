# Introduction

## Dependencies

# Stepwise guide
## Removing the parentheses in filename are not accepted. 
When uploading the tracefiles to the Barcode of Life Database (https://boldsystems.org/) parentheses are not allowed. Therefore we start by removing these from the filenames. 
	
	for i in *'('*')'*
      do
        mv "$i" `echo $i | sed 's/\(.*\) (.*).\(.*\)/\1.\2/'`
      done
    
## de novo assembly of electropherograp files from Sanger sequencing
bi-directional sequencing. First reverse, for 3'-5' direction output. For the assembly we use the open-source software [Tracy](https://github.com/gear-genomics/tracy) (Rausch et al. 2020).
	
	/fileserver/tracy \
      assemble \
        /fileserver/raw/SAMPLE-NUMBER_COIIR*.ab1 \
        /fileserver/raw/SAMPLE-NUMBER_COIIF*.ab1 \
        -o /fileserver/intermediate/SAMPLE-NUMBER

## clipping primer sequences
We use the software [Cutadapt](https://github.com/marcelm/cutadapt) (Martin, 2011)
	
	cutadapt \
      -a "ATGGCAGATTAGTGCAATGA...CAAGTACTGGTCTCTTAAAC;max_error_rate=0.15;" \
      --discard-untrimmed \
      -o /fileserver/intermediate/SAMPLE-NUMBER.cons.trimmed.fa \
      /fileserver/intermediate/SAMPLE-NUMBER.cons.fa
    
## rename sequence headers
	
	sed \
      "s/>.*/>SAMPLE-NUMBER/" \
      /fileserver/intermediate/SAMPLE-NUMBER.cons.trimmed.fa > \
      /fileserver/final/SAMPLE-NUMBER.clean.fa
    
 ## combine fasta files into single file
	
	cat /fileserver/final/*.fa > /fileserver/final/combined.fa
 
# Pipeline

	
	cd /fileserver/
	
	mkdir ./raw/
	mkdir ./intermediate/
	mkdir ./final/
	
move the succesfull electropherograp files (.ab1) files to the "raw" folder
  
	cd ./raw/
  
	for i in *'('*')'*
      do
        mv "$i" `echo $i | sed 's/\(.*\) (.*).\(.*\)/\1.\2/'`
      done
   
	cd ..
   
	SAMPLES=sample_no.txt
   
	for i in $(cat $SAMPLES)
	    do
	        /tracy/bin/tracy assemble ./raw/"$i"_COIIR*.ab1 ./raw/"$i"_COIIF*.ab1 -o ./intermediate/"$i"
        
          cutadapt -a "ATGGCAGATTAGTGCAATGA...CAAGTACTGGTCTCTTAAAC;max_error_rate=0.15;" --discard-untrimmed -o ./intermediate/"$i".cons.trimmed.fa ./intermediate/"$i".cons.fa
	      
          sed "s/>.*/>${i}/" ./intermediate/"$i".cons.trimmed.fa > ./final/"$i".clean.fa
	    done
      
	cat ./final/*.fa > ./final/combined.fa

# References
- Rausch, T., Fritz, M.HY., Untergasser, A. et al. Tracy: basecalling, alignment, assembly and deconvolution of sanger chromatogram trace files. BMC Genomics 21, 230 (2020). https://doi.org/10.1186/s12864-020-6635-8
- Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], 17(1), 10-12, may 2011. https://doi.org/10.14806/ej.17.1.200.
