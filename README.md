# Mites-assembly
Assembly Ppr, Nps, Hga with Pacbio HIFI, Tell-seq, Hi-c
## 1.hifiasm assembly the Pacbio HIFI data
  Ppr for example
  
    /home/shangao/Software/Assembly/hifiasm/hifiasm -o test.asm -t 32 /home/shangao/Data/mites/reads/PacBio/Ppr/m64093_200831_134054.Q20.fastq.gz

## 2.Tell-seq tools for the linked reads
  Using the docker version of Tell-reads, Tell-link, Tell-sort 
  
### 2.1 GenerateGenomeIndex

    /home/shangao/Software/TELLSeq/prerelease/tellread-release/generateGenomeIndexBed.sh ./Ppr/Ppr.fa
    
### 2.2 Tell-read filter the raw data

    /home/shangao/Software/TELLSeq/prerelease/tellread-release/run_tellread.sh \
	    -i /home/shangao/Data/bcl_Tell-seq/raw_data_new/201211_A00685_0102_BHWLNMDRXX_Dezember4/ \
	    -o /home/shangao/Data/hifiasm_tell-seq/test_Ppr_new \
	    -f /home/shangao/Scratch/hifiasm/genome \
	    -s T501 \
	    -g Ppr_new
  
### 2.3 Convert the reads th linked Reads

    /home/shangao/Software/TELLSeq/conversion_tool/ust10x \
    	-sz 16000000 \
    	-i1 /home/shangao/Data/gaoshan/hifiasm_tell-seq/test_Ppr_new/Full/test_Ppr_I1_T501.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
    	-r1 /home/shangao/Data/gaoshan/hifiasm_tell-seq/test_Ppr_new/Full/test_Ppr_R1_T501.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
    	-r2 /home/shangao/Data/gaoshan/hifiasm_tell-seq/test_Ppr_new/Full/test_Ppr_R2_T501.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
	    -wl /home/shangao/Software/TELLSeq/conversion_tool/4M-with-alts-february-2016.txt
  
    mv R1_sl.fastq.gz.4tenx.fastq T501_R1_sl.fastq.gz.4tenx.fastq
    mv R2_sl.fastq.gz.4tenx.fastq T501_R2_sl.fastq.gz.4tenx.fastq
    pigz *_sl.fastq.gz.4tenx.fastq

## 3. Scaff10x improve the primary assembly

### 3.1 Scaff_reads convert the linked reads

    /home/shangao/Software/Assembly/Scaff10X/src/scaff_reads -nodes 30 ./input.dat genome-BC_1.fastq.gz genome-BC_2.fastq.gz
    
### 3.2 Scaff10X improved the assembly

    /home/shangao/Software/Assembly/Scaff10X/src/scaff10x \
	    -nodes 30 -longread 1 -gap 100 -matrix 2000 -reads 10 -score 10 -edge 50000 -link 8 -block 50000 -plot Ppr_result/barcode_lengtg.png \
	    /home/shangao/Scratch/gaoshan/hifiasm/genome/Ppr_new/Ppr_new.fa \
	    genome-BC_1.fastq.gz \
	    genome-BC_2.fastq.gz \
	    output_scaffolds.fasta
      
## 4. bbstats for N50 .etc

    ~/Software/tools/bbmap/bbstats.sh in=output_scaffolds.fasta > bbstats.out
    
## 5. Busco 

    conda activate Busco

    busco -i .fa -m genome -o busco --auto-lineage-euk

## 6. Blobtools2 

    conda activate btk_env

### 6.1 Create 

    /home/shangao/Software/blobtoolkit/blobtools2/blobtools create \
	    --fasta .fa \
	    /home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/blobtools2/dataset/Ppr
      
### 6.2 Add cov, hits, busco results

    /home/shangao/Software/blobtoolkit/blobtools2/blobtools add \
    	--taxrule bestsumorder \
    	--taxdump /RAID/Data/databases/taxdump/ \
    	--cov /home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/Ppr_Pacbio1.bam \
    	--hits /home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/Ppr.ncbi.blastn.out \
      --busco /home/shangao/Scratch/hifiasm/test_Ppr/busco_p/run_arachnida_odb10/full_table.tsv \
	    --busco /home/shangao/Scratch/hifiasm/test_Ppr/busco_p/run_eukaryota_odb10/full_table.tsv \
    	/home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/blobtools2/dataset/Ppr
      
### 6.3 Filter

    /home/shangao/Software/blobtoolkit/blobtools2/blobtools filter \
 	    --param length--Min=1000 \
    	--summary-rank species \
     	--fasta ../test.p_ctg.fa \
    	/home/shangao/Scratch/gaoshan/hifiasm/test_Ppr/blobtools2/dataset/Ppr
      
### 6.4 View

    /home/shangao/Software/blobtoolkit/blobtools2/blobtools view --remote /home/shangao/Scratch/hifiasm/test_Ppr/blobtools2/dataset/Ppr

## 7. Hi-c

    conda activate hi-c

### 7.1 Dovetail_tools to get valid.pairs reads

    /home/shangao/Software/dovetail_tools/omni-c_qc.bash ~/Scratch/gaoshan/scaff10x/Ppr_result/output_scaffolds.fasta /home/shangao/Scratch/gaoshan/hi-c/Omni-C/X204SC21021794-Z01-F001/raw_data/Ppr1/Ppr1_FKDL210026677-1a_HWLYJDSXY_L4_1.fq.gz /home/shangao/Scratch/gaoshan/hi-c/Omni-C/X204SC21021794-Z01-F001/raw_data/Ppr1/Ppr1_FKDL210026677-1a_HWLYJDSXY_L4_2.fq.gz Ppr Ppr1 40
    
    /home/shangao/Software/dovetail_tools/contact_map.sh ./Ppr.PT.pairs.gz ~/Scratch/gaoshan/scaff10x/Ppr_result/output_scaffolds.fasta Ppr 50

### 7.2 SALSA for improve assembly

    python run_pipeline.py -a contigs.fasta -l contigs.fasta.fai -b alignment.bed -e DNASE -o scaffolds 

#### 8.HGT
	
	conda activate BRAKER
	
### 8.1 get pro.seq from genome by gff

	gffread ../01braker/braker/braker.gtf -g /home/shangao/Scratch/hi-c/sala2/Ppr2/scaffolds/scaffolds_FINAL.fasta -y Ppr.pep

### 8.2 prepare the db and tax files

	bd: /home/shangao/Data/HGT_db/
	perl -lane 'if(/^>(\w+)\s.+TaxID\=(\d+)/){print "$1 $2"}' <(zcat uniref90.fasta.gz) | gzip > uniref90.fasta.taxlist.gz
	diamond blastp --sensitive --index-chunks 1 -k 500 -e 1e-5 -p 32 -q Ppr.pep -d /home/shangao/Data/HGT_db/uniref90.dmnd -a Ppr_diamond_results
	cat <(zcat /home/shangao/Data/HGT_db/uniref90.fasta.taxlist.gz) <(diamond view -a Ppr_diamond_results.daa)|perl -lane 'if(@F==2){$tax{$F[0]}=$F[1];}else{if(exists($tax{$F[1]})){print join("\t",@F,$tax{$F[1]});} else {print join("\t",@F,"NA");}}'| gzip > diamond_results.daa.taxid.gz
	
### 8.3 get results

	/home/shangao/script/hgt/diamond_to_HGT_candidates.pl -i diamond_results.daa.taxid.gz -f Ppr.pep -p /home/shangao/Data/HGT_db/taxdump
	
The program outputs 3 files, suffixed with the tags:

  HGT_results: hU, AI and CHS scores for all query proteins.
  HGT_candidates: All queries which show evidence from both hU (or AI) and CHS above the specifed thresholds (defaults are >=30 for hU, >=90% for CHS).
  HGT_warnings: Any non-fatal warnings detected during the run. Worth checking, but most warnings can probably safely be ignored.
  
 when you meet some perl lib problem get out of conda
 
 ### 8.4 HGT_candidates_to_fasta
 
 	/home/shangao/script/hgt/HGT_candidates_to_fasta.pl -i diamond_results.daa.taxid.gz -c diamond_results.daa.taxid.gz.HGT_candidates.Metazoa.hU30.CHS90.txt -u /home/shangao/Data/HGT_db/uniref90.fasta -f Ppr.pep -p /home/shangao/Data/HGT_db/taxdump
	
	mkdir ../mafft_alns
	for f in *.fasta; do echo $f; mafft --auto --quiet --thread 8 $f > ../mafft_alns/${f}.mafft; done
	
	mkdir processed_files
	COUNT=1;
	
	## iqtree commands
	for file in *mafft;
	   do echo $COUNT: $file;
	   iqtree-omp -s $file -st AA -nt 16 -quiet -bb 1000 -m TESTNEW -msub nuclear
	   mv $file processed_files/;
	   COUNT=$[COUNT+1];
	done
	
	mkdir iqtree treefiles
	mv *treefile treefiles/
	mv *bionj *gz *contree *iqtree *log *mldist *model *nex iqtree/
	
## 9. JCVI_Synteny
		https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)
		https://github.com/tanghaibao/jcvi
	
	
	
	
