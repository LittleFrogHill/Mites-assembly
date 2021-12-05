# Mites-assembly
Assembly Ppr, Nps, Hga with Pacbio HIFI, Tell-seq, Hi-c
## 1.hifiasm assembly the Pacbio HIFI data
	Ppr for example
	  HIFI+Hi-c data
	  
	  	/home/shangao/software/hifiasm-0.16.1/hifiasm -o test -t 50 --h1 /home/shangao/Scratch/hi-c/dovetail/Ppr_new/merge/Ppr1_1.fq.gz --h2 /home/shangao/Scratch/hi-c/dovetail/Ppr_new/merge/Ppr1_2.fq.gz /home/shangao/Data/../mites/reads/PacBio/Ppr/m64093_200831_134054.Q20.fastq.gz
	 	
	  Only HIFI
		
		/home/shangao/software/hifiasm-0.16.1/hifiasm -o test -t 50 /home/shangao/Data/../mites/reads/PacBio/Ppr/m64093_200831_134054.Q20.fastq.gz

## 2.Tell-seq tools for the linked reads
  Using the docker version of Tell-reads, Tell-link, Tell-sort 
  
### 2.1 GenerateGenomeIndex

    /home/shangao/Software/TELLSeq/prerelease/tellread-release/generateGenomeIndexBed.sh ./Ppr/Ppr.fa
    
### 2.2 Tell-read filter the raw data
	    
	/home/shangao/Software/TELLSeq/prerelease/tellread-release/run_tellread.sh \
	-i /RAID/Data/mites/reads/TELLSeq/bcl_Tell-seq/raw_data01/201211_A00685_0102_BHWLNMDRXX_Dezember4/ \
	-o /RAID/Data/shangao/hifiasm_tell-seq/Ppr/test/01 \
	-f /home/shangao/Data/hifiasm/test_Ppr/new_hifiasm/0.16_new/only_hifi/genome \
	-s T502 \
	-g Ppr
	
	/home/shangao/Software/TELLSeq/prerelease/tellread-release/run_tellread.sh \
	-i /RAID/Data/mites/reads/TELLSeq/bcl_Tell-seq/raw_data02/JB01/210813_A00620_0167_AHHFG5DRXY_August6/ \
	-o /RAID/Data/shangao/hifiasm_tell-seq/Ppr/test/02 \
	-f /home/shangao/Data/hifiasm/test_Ppr/new_hifiasm/0.16_new/only_hifi/genome \
	-s T502 \
	-g Ppr

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
    
    /home/jbast/anaconda3/envs/EDTA/bin/blastn -db /RAID/Data/databases/nt/nt \
       -query ../Ppr_instagrall.polished.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 32 \
       -out Ppr.ncbi.blastn.out
      
    /home/shangao/Software/minimap2/minimap2 -ax map-pb \
	-t 40 ../Ppr_instagrall.polished.fa \
	/home/shangao/Data/../mites/reads/PacBio/Ppr/m64093_200831_134054.Q20.fastq.gz \
	| samtools sort -@16 -O BAM -o Ppr_Pacbio1.bam -

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

### 7.2 improve assembly
	SALSA
    python run_pipeline.py -a contigs.fasta -l contigs.fasta.fai -b alignment.bed -e DNASE -o scaffolds 

	3D-DNA
	
	python /NVME/Software/3d-dna/juicer/misc/generate_site_positions.py HindIII Ppr_instagrall Ppr_instagrall.fa
	
        awk 'BEGIN{OFS="\t"}{print $1, $NF}' Ppr_instagrall_HindIII.txt > Ppr_instagrall.chrom.sizes
	
	/NVME/Software/3d-dna/juicer/scripts/juicer.sh -g Ppr -s none \
	-z reference/Ppr.fa -t 40 \
	-p reference/Ppr.genome.chrom.size \
	-D /NVME/Software/3d-dna/juicer 
	
	/NVME/Software/3d-dna/run-asm-pipeline.sh -r 0 reference/Ppr.fa aligned/merged_nodups.txt

	seqkit grep -v -f filter.list Ppr.FINAL.sort.fasta > Ppr.FINAL.fa
	seqkit fx2tab Ppr.FINAL.fa -l -g -n -i -H
	![image](https://user-images.githubusercontent.com/34407101/144658065-e1b75f1e-97ff-4be4-808c-c085ce8a5efc.png)


#### 8. BRAKER
### 8.1 mapping
	
		/home/shangao/Software/STAR-2.7.8a/source/STAR \
		    --runThreadN 40 \
		    --runMode genomeGenerate \
		    --genomeDir STAR \
		    --genomeSAindexNbases 12 \
		    --genomeFastaFiles /home/shangao/Data/EDTA/Ppr/Ppr_instagrall/Ppr_instagrall.polished.fa.mod.MAKER.masked 
		
		/home/shangao/Software/STAR-2.7.8a/source/STAR \
			--genomeDir STAR \
			--runThreadN 40 \
			--readFilesIn /RAID/Data/mites/transcriptomes/RNA_annotation/150360_R1.fq.gz, /RAID/Data/mites/transcriptomes/RNA_annotation/150360_R2.fq.gz \
			--readFilesCommand zcat \
			--outFileNamePrefix Ppr \
			--outSAMtype BAM SortedByCoordinate \
			--outBAMsortingThreadN 10 \
			--outSAMstrandField intronMotif \
			--outFilterIntronMotifs RemoveNoncanonical 
		
		#ulimit -n 10000
		
		/home/shangao/Software/STAR-2.7.8a/source/STAR \
		        --genomeDir STAR \
		        --runThreadN 40 \
		        --readFilesIn /RAID/Data/mites/transcriptomes/RNA_annotation/clean_data/SRR4039022_1_val_1.fq.gz, /RAID/Data/mites/transcriptomes/RNA_annotation/clean_data/SRR4039022_2_val_2.fq.gz \
		        --readFilesCommand zcat \
		        --outFileNamePrefix Ppr_down \
		        --outSAMtype BAM SortedByCoordinate \
		        --outBAMsortingThreadN 10 \
		        --outSAMstrandField intronMotif \
		        --outFilterIntronMotifs RemoveNoncanonical

### 8.2 BRAKER

	conda activate BRAKER	
	
	braker1:
	braker.pl --cores 40 --species=Ppr --genome=/home/shangao/Data/EDTA/Ppr/Ppr_instagrall/Ppr_instagrall.polished.fa.mod.MAKER.masked \
	 --softmasking --bam=/home/shangao/Scratch/breaker/000rna-seq/Ppr/Ppr_instagrall/PprAligned.sortedByCoord.out.bam \
	 --gff3 --useexisting --GENEMARK_PATH=/home/shangao/Software/BrAKER/gmes_linux_64/
	
	braker2:
	braker.pl --cores 40 --species=Ppr1 --genome=/home/shangao/Data/EDTA/Ppr/Ppr_instagrall/Ppr_instagrall.polished.fa.mod.MAKER.masked \
	 --softmasking \
	 --gff3 --GENEMARK_PATH=/home/shangao/Software/BrAKER/gmes_linux_64/ \
	 --prot_seq /RAID/Data/databases/braker_db/proteins.fasta --epmode --workingdir=braker2_out --PROTHINT_PATH=/home/shangao/software/ProtHint-2.6.0/bin/

	TSEBRA:
	#/home/shangao/Software/TSEBRA/bin/fix_gtf_ids.py --gtf ../../test1/braker/braker.gtf --out braker1.fix.gtf
	
	#/home/shangao/Software/TSEBRA/bin/fix_gtf_ids.py --gtf ../braker2_out/braker.gtf --out braker2.fix.gtf
	
	the hyperparameters themselves.  Evidence sources in a standard BRAKER1 and BRAKER2 output138are:  protein database (P), EST database (E), combined EST/protein database (C), and manual139anchored (M). The default weights for these arewP= 0.1,wE= 10,wC= 5 andwM= 1.  A140transcript has low evidence support in this default setting if the fraction of supported introns is less141than 0.75 and the supported start/stop-codon fraction is less than 1.0.  The score specific thresholds142are•1= 0,•2= 0.5,•3= 25,•4= 10.  We have shown that TSEBRA using default parameters143performs with high accuracy across several species, see Results and discussion  
	
	/home/shangao/Software/TSEBRA/bin/tsebra.py -g braker1.fix.gtf,braker2.fix.gtf \
		-c /home/shangao/Software/TSEBRA/config/default.cfg \
		-e ../../test1/braker/hintsfile.gff,../braker2_out/hintsfile.gff \
        	-o braker1+2_combined.gtf
	
	cat default.cfg
		# Weight for each hint source
		# Values have to be >= 0
		P 0.1
		E 20
		C 5
		M 1
		# Required fraction of supported introns or supported start/stop-codons for a transcript
		# Values have to be in [0,1]
		intron_support 0.75
		stasto_support 1
		# Allowed difference for each feature 
		# Values have to be in [0,1]
		e_1 0
		e_2 0.5
		# Values have to be >0
		e_3 25
		e_4 10
	
#### 9.HGT  https://github.com/reubwn/hgt
	
	conda activate BRAKER
	
### 9.1 get pro.seq from genome by gff

	gffread ../01braker/braker/braker.gtf -g /home/shangao/Scratch/hi-c/sala2/Ppr2/scaffolds/scaffolds_FINAL.fasta -y Ppr.pep

### 9.2 prepare the db and tax files

	bd: /home/shangao/Data/HGT_db/
	perl -lane 'if(/^>(\w+)\s.+TaxID\=(\d+)/){print "$1 $2"}' <(zcat uniref90.fasta.gz) | gzip > uniref90.fasta.taxlist.gz
	diamond blastp --sensitive --index-chunks 1 -k 500 -e 1e-5 -p 32 -q Ppr.pep -d /home/shangao/Data/HGT_db/uniref90.dmnd -a Ppr_diamond_results
	cat <(zcat /home/shangao/Data/HGT_db/uniref90.fasta.taxlist.gz) <(diamond view -a Ppr_diamond_results.daa)|perl -lane 'if(@F==2){$tax{$F[0]}=$F[1];}else{if(exists($tax{$F[1]})){print join("\t",@F,$tax{$F[1]});} else {print join("\t",@F,"NA");}}'| gzip > diamond_results.daa.taxid.gz
	
### 9.3 get results

	/home/shangao/script/hgt/diamond_to_HGT_candidates.pl -i diamond_results.daa.taxid.gz -f Ppr.pep -p /home/shangao/Data/HGT_db/taxdump
	
The program outputs 3 files, suffixed with the tags:

  HGT_results: hU, AI and CHS scores for all query proteins.
  HGT_candidates: All queries which show evidence from both hU (or AI) and CHS above the specifed thresholds (defaults are >=30 for hU, >=90% for CHS).
  HGT_warnings: Any non-fatal warnings detected during the run. Worth checking, but most warnings can probably safely be ignored.
  
 when you meet some perl lib problem get out of conda
 
 ### 9.4 HGT_candidates_to_fasta
 
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
	
	GFF format
	chr1    AUGUSTUS        gene    2843    3493    .       +       .       ID=Ppr_g1
	chr1    AUGUSTUS        transcript      2843    3493    .       +       .       Parent=Ppr_g1;ID=Ppr_g1;transcriptID=P
	chr1    AUGUSTUS        start_codon     2843    2845    .       +       0       Parent=Ppr_g1;ID=Ppr_g1;transcriptID=P
	chr1    AUGUSTUS        exon    2843    3493    .       +       .       Parent=Ppr_g1;ID=Ppr_g1;transcriptID=Ppr_g1.t1
	chr1    AUGUSTUS        CDS     2843    3493    0.56    +       0       Parent=Ppr_g1;ID=Ppr_g1;transcriptID=Ppr_g1.t1
	chr1    AUGUSTUS        stop_codon      3491    3493    .       +       0       Parent=Ppr_g1;ID=Ppr_g1;transcriptID=P
	chr1    AUGUSTUS        gene    3581    4192    .       +       .       ID=Ppr_g2

	/home/shangao/script/hgt/get_location_of_HGT_candidates.pl -i diamond_results.daa.taxid.gz.HGT_results.Metazoa.txt -g  /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/test6_TSEBRA/TSEBRA/braker1+2_combined_rmHiC_scaffold_10_changeHiC_scaffold_11.sort.rename.gff -f Ppr.pep

	Rscript /home/shangao/script/hgt/analyse_trees.R -q anno -p ~/Scratch/breaker/04hgt/Ppr/instagraal/TSEBRA/mafft_alns/treefiles

	/home/bonkis/anaconda2/pkgs/python-2.7.16-h9bab390_7/bin/python /home/shangao/software/gff3sort/Gtf2GFF.py braker1+2_combined_rmHiC_scaffold_10_changeHiC_scaffold_11.sort.rename.gtf
	
## 10. JCVI_Synteny
		https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)
		https://github.com/tanghaibao/jcvi

	/home/shangao/lastz-distrib/bin/lastz ../test.p_ctg.fa[multiple] ../test.asm.a_ctg.fa --gfextend --chain --gapped --format=blastn > lastzfomat6.out.1
	
## 11. circos 

![image](https://user-images.githubusercontent.com/34407101/144664967-9b099105-f9b1-47cb-94c8-2d6dfd39dd85.png)
> require(RIdeogram)
Loading required package: RIdeogram
> human_karyotype <- read.table("mites_karyotype", sep = "\t", header = T, stringsAsFactors = F)
> gene_density <- GFFex(input = "/home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/test6_TSEBRA/TSEBRA/braker1+2_combined_rmHiC_scaffold_10_changeHiC_scaffold_11.sort.rename.gtf", karyotype = "mites_karyotype", feature = "gene", window = 100000)
> gene_density <- GFFex(input = "/home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/test6_TSEBRA/TSEBRA/braker1+2_combined_rmHiC_scaffold_10_changeHiC_scaffold_11.sort.rename.gtf", karyotype = "mites_karyotype", feature = "gene", window = 100000)
> ideogram(karyotype = human_karyotype, overlaid = gene_density)
> convertSVG("chromosome_100k.svg", device = "png")
Error in read_data(svg) : 
  Argument 'svg' or 'css' must be a file path, url, or raw vector.
> ideogram(karyotype = human_karyotype, overlaid = gene_density)
> convertSVG("chromosome_100k.svg", device = "png")
Error in read_data(svg) : 
  Argument 'svg' or 'css' must be a file path, url, or raw vector.
> ideogram(karyotype = human_karyotype, overlaid = gene_density)
> convertSVG("chromosome.svg", device = "png")


## 12. SV
