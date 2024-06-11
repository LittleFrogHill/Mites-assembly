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
#### build reference(chrom size, assembly and bwa index) fastq(_R1.fastq.gz) two folder
#### prepare genome files
	python /NVME/Software/3d-dna/juicer/misc/generate_site_positions.py HindIII Ppr_instagrall Ppr_instagrall.fa
        awk 'BEGIN{OFS="\t"}{print $1, $NF}' Ppr_instagrall_HindIII.txt > Ppr_instagrall.chrom.sizes
	
	~/Software/QC/TrimGalore-0.6.5/trim_galore -j 30 -q 30 --fastqc --paired --output_dir fastq ./Hga_1.fq.gz ./Hga_2.fq.gz
	
	/NVME/Software/3d-dna/juicer/scripts/juicer.sh -g Ppr -s none \
	-z reference/Ppr.fa -t 40 \
	-p reference/Ppr.genome.chrom.size \
	-D /NVME/Software/3d-dna/juicer 
	
	/NVME/Software/3d-dna/run-asm-pipeline.sh -r 0 reference/Ppr.fa aligned/merged_nodups.txt

	seqkit grep -v -f filter.list Ppr.FINAL.sort.fasta > Ppr.FINAL.fa
	seqkit fx2tab Ppr.FINAL.fa -l -g -n -i -H
	![image](https://user-images.githubusercontent.com/34407101/144658065-e1b75f1e-97ff-4be4-808c-c085ce8a5efc.png)

### 7.3 Kat and genomescope
	kat comp -t 8 -o 022812 /RAID/Data/Mites/Reads/PacBio/Ppr/238/m64093_220206_022812.hifi_reads.fastq.gz purged.fa
	
	/NVME/Software/KMC/kmc -k21 -t10 -m64 -ci1 -cs10000 -fm 110050 /RAID/Data/Mites/Reads/PacBio/Ppr/238/m64093_220131_110050.hifi_reads.fastq.gz ./tmp/
	for i in 022812 110050
	do
	/NVME/Software/KMC/kmc_tools transform ${i}_csi histogram ${i} -cx10000

	L=$(/NVME/Software/smudgeplot/exec/smudgeplot.py cutoff ${i}_csi_k21.hist L)
	U=$(/NVME/Software/smudgeplot/exec/smudgeplot.py cutoff ${i}_csi_k21.hist U)
	echo $L $U

	/NVME/Software/KMC/kmc_tools transform ${i}_csi -ci"$L" -cx"$U" dump -s ${i}_csi_L"$L"_U"$U".dump
	/NVME/Software/smudgeplot/exec/smudgeplot.py hetkmers -o ${i}_csi_L"$L"_U"$U" < ${i}_csi_L"$L"_U"$U".dump

	/NVME/Software/smudgeplot/exec/smudgeplot.py plot ${i}_csi_L"$L"_U"$U"_coverages.tsv
	done
	
### 7.4 GenomeMask
	RepeatMasker -lib ath-families.fa -e ncbi -dir ../RepeatModeler /home/shangao/Scratch/hi-c/sala2/Ppr2/scaffolds/scaffolds_FINAL.fasta
	
	perl ~/Software/EDTA/EDTA.pl --genome /home/shangao/Scratch/hi-c/juicer/Ppr/Ppr_withoutchange/review/Ppr.FINAL.sort.fasta --sensitive 1 --anno 1  --threads 50 --overwrite 1

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
	
#### using hardmasking assembly
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
#### using softmasking assembly
	cat braker.sh 
	mkdir braker1_out
	braker.pl --cores 40 --species=Ppr --genome=/RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.softmask.fa \
		 --softmasking --bam=/home/shangao/Scratch/breaker/000rna-seq/Ppr/Ppr_instagrall/changed_chr/PprAligned.sortedByCoord.out.bam \
		 --gff3 --useexisting --GENEMARK_PATH=/home/shangao/Software/BrAKER/gmes_linux_64/ --workingdir=braker1_out
	
	mkdir braker2_out
	
		braker.pl --cores 40 --species=Ppr_gene --genome=/RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.softmask.fa \
		 --softmasking \
		 --gff3 --GENEMARK_PATH=/home/shangao/Software/BrAKER/gmes_linux_64/ \
		 --prot_seq /RAID/Data/databases/braker_db/proteins.fasta --epmode --workingdir=braker2_out --PROTHINT_PATH=/home/shangao/software/ProtHint-2.6.0/bin/
	
	/home/shangao/Software/TSEBRA/bin/fix_gtf_ids.py --gtf braker1_out/braker.gtf --out braker1.fix.gtf
	
	/home/shangao/Software/TSEBRA/bin/fix_gtf_ids.py --gtf braker2_out/braker.gtf --out braker2.fix.gtf
	
	/home/shangao/Software/TSEBRA/bin/tsebra.py -g braker1.fix.gtf,braker2.fix.gtf \
		-c /home/shangao/Software/TSEBRA/config/default.cfg \
		-e braker1_out/hintsfile.gff,braker2_out/hintsfile.gff \
	        -o braker1+2_combined.gtf

#### Add UTR region
	mkdir braker_UTR_out
	braker.pl --cores 40 --species=Ppr_UTR --genome=/RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.softmask.fa \
	 --softmasking --bam=/home/shangao/Scratch/breaker/000rna-seq/Ppr/Ppr_instagrall/changed_chr/PprAligned.sortedByCoord.out.bam \
	  --UTR=on --useexisting --GENEMARK_PATH=/home/shangao/Software/BrAKER/gmes_linux_64/ --workingdir=braker_UTR_out
	  
	cat ../../braker2_out/add_UTR.sh
	braker.pl --genome=/RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.softmask.fa --addUTR=on --softmasking \
    --bam=/home/shangao/Scratch/breaker/000rna-seq/Ppr/Ppr_instagrall/changed_chr/PprAligned.sortedByCoord.out.bam --workingdir=./ \
    --AUGUSTUS_hints_preds=augustus.hints.gtf --cores=8 \
    --skipAllTraining --species=Ppr_gene

	/home/shangao/Software/TSEBRA/bin/fix_gtf_ids.py --gtf ../braker_utr.gtf --out braker1.fix.gtf
	
	/home/shangao/Software/TSEBRA/bin/tsebra.py -g braker1.fix.gtf,/home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/softmask/braker2_out/augustus.hints_utr.gtf \
	-c /home/shangao/Software/TSEBRA/config/default.cfg \
	-e ../hintsfile.gff,../../braker2_out/hintsfile.gff \
        -o braker1+2_combined.gtf

#### change prefix
	/NVME/Software/TSEBRA/bin/rename_gtf.py --gtf braker1+2_combined.gtf --prefix hap1 --out hap1.gtf

#### gtf2gff
	/home/bonkis/anaconda2/pkgs/python-2.7.16-h9bab390_7/bin/python /home/shangao/software/gff3sort/Gtf2GFF.py braker1+2_combined_rmHiC_scaffold_10_changeHiC_scaffold_11.sort.rename.gtf
	sed -i 's/ID/name/g' hap1.gff
	sed -i 's/ID/name/g' hap2.gff
	
	java -cp /home/shangao/software/TBtools-1.098745/TBtools_JRE1.6.jar biocjava.bioIO.FastX.FastaIndex.FastaLongestRepresentater --inFasta Ppr_hap2.pep --outFasta Ppr_hap2_longest_protein.fa
	

	

#### eggmaper
	cat eggmaper.sh 
	export EGGNOG_DATA_DIR=/RAID/Data/databases/eggnog-mapper-data
	#emapper.py -i Ppr.pep -o test
	#emapper.py -m mmseqs --itype CDS --translate -i Ppr.cds -o test --decorate_gff Ppr.gff --decorate_gff_ID_field ID --resume
	emapper.py -m mmseqs --itype CDS --translate -i Ppr.cds -o cds --decorate_gff Ppr.gff --decorate_gff_ID_field ID
	emapper.py -i Ppr.pep --itype proteins -o pro --decorate_gff Ppr.gff --decorate_gff_ID_field ID > emapper.pro.log
	
##### add UTR eggmaper

	emapper.py -m mmseqs --itype CDS -i Ppr.cds -o cds --decorate_gff Ppr.gff --decorate_gff_ID_field ID --resume
	
#### clusterprofilter
	https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz
	http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
	https://guangchuangyu.github.io/cn/2017/09/ko-id-conversion/
	#利用注释信息构建数据库进行富集分析
	#motoko上R包的安装，注意是哪一个R，在conda环境外安装
	
	require(DOSE)
	require(clusterProfiler)
	library(enrichplot)
	library(ggplot2)
	library(stringr)
	
	![image](https://user-images.githubusercontent.com/34407101/161494391-93458d13-1e89-49fc-94b6-8573de74f331.png)

	
##### Go database build https://blog.csdn.net/liuninghua521/article/details/111030304
	#protein database
	egg<-read.delim("pro.emapper.annotations.test")
	#cds database
	egg<-read.delim("cds.emapper.annotations.test")
		
	gene_ids <- egg$query
	eggnog_lines_with_go <- egg$GOs!='-'
	eggnog_annoations_go <- str_split(egg[eggnog_lines_with_go,]$GOs, ",")
	gene_to_go <- data.frame(gene = rep(gene_ids[eggnog_lines_with_go], times = sapply(eggnog_annoations_go, length)), term = unlist(eggnog_annoations_go))
	term2gene1 <- gene_to_go[, c(2, 1)]
	term2gene <- buildGOmap(term2gene1)
	go2ont <- go2ont(term2gene$GO)
	go2term <- go2term(term2gene$GO)
	
##### kegg database build https://www.jieandze1314.com/post/cnposts/208/
	gene2ko <- egg %>%
			dplyr::select(GID = query, KO = KEGG_ko) %>%
			na.omit()
	download 'ko00001.json' and creat 'kegg_info.RData' , https://www.genome.jp/kegg-bin/get_htext?ko00001
	
	if(!file.exists('kegg_info.RData')){
	   library(jsonlite)
	   library(purrr)
	   library(RCurl)
	   
	   update_kegg <- function(json = "ko00001.json",file=NULL) {
	     pathway2name <- tibble(Pathway = character(), Name = character())
	     ko2pathway <- tibble(Ko = character(), Pathway = character())
	     
	     kegg <- fromJSON(json)
	     
	     for (a in seq_along(kegg[["children"]][["children"]])) {
	 	      A <- kegg[["children"]][["name"]][[a]]
	       
	       for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
	         B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
	         
	         for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
	           pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
	           
	           pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
	           pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
	           pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
	           
	           kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
	           
	           kos <- str_match(kos_info, "K[0-9]*")[,1]
	           
	           ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
	         }
	       }
	     }
	     
	     save(pathway2name, ko2pathway, file = file)
	   }
	   
	   update_kegg(json = "ko00001.json",file="kegg_info.RData")
	   
	+ }
	#要load一遍，数据存在当前目录
	load('./kegg_info.RData')
	
	colnames(ko2pathway)=c("KO",'Pathway')
	library(stringr)
	gene2ko$KO=str_replace(gene2ko$KO,"ko:","")
	
	library(dplyr)
	
	gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "KO") %>% 
	   dplyr::select(GID, Pathway) %>%
	   na.omit()
	 head(gene2pathway)
	
#### kegg database with python script /home/shangao/script/python/divide_eggnoger_results2clusterprofliter.py
	kegg <- read.delim("pro.emapper.annotations_Ko.term")
	ko2name(kegg$keggId) -> y

#### Reactome 
	awk -vFS='\t' -vOFS='\t'  '$15!="-"&&$15!=""{print$1,$15}' hap0_UTR_longest.f.pep_1.tsv.test|sort |uniq > hap0_UTR_longest.f.pep_1.tsv.testforReactomeBuild
	python /home/shangao/script/python/01dea/10Build_ReactomeDB_from_interproscan.py -t hap0_UTR_longest.f.pep_1.tsv.testforReactomeBuild -o hap0_UTR_longest.f.pep_1.tsv.testforReactomeBuild.result
	sort hap0_UTR_longest.f.pep_1.tsv.testforReactomeBuild.result|uniq|grep 'Rea' > hap0_UTR_longest.f.pep_1.tsv.testforReactomeBuild.result.Reactom

	egg<-read.delim("ReactomePathways.txt",header=F)
	rea2gene<-read.delim("hap0_UTR_longest.f.pep_1.tsv.testforReactomeBuild.result.Reactome.head")
	rea2descrption<-read.delim("ReactomePathways.txt",header=F)
	colnames(rea2descrption) <- c("reaID", "pathway", "animals")
	head(rea2descrption)
	head(rea2gene)
	hgt_rea <- enricher(hgt_gene, TERM2GENE = rea2gene[,c("items","geneID")], TERM2NAME = rea2descrption[,c("reaID","pathway")], pAdjustMethod = "BH",pvalueCutoff  = 0.05, qvalueCutoff  = 0.2)
	write.table(as.data.frame(ego),"go_enrich.csv",sep="\t",row.names =F,quote=F)


#### enrich
	hgt<-read.delim("/home/shangao/Scratch/breaker/04hgt/Ppr/softmask/diamond_results.daa.taxid.gz.HGT_candidates.Metazoa.hU30.CHS90.txt.genelist")
	hgt<-as.matrix(hgt)
	x = enricher(hgt, TERM2GENE=kegg, TERM2NAME=y)
	
	dea<-read.delim("/home/shangao/Scratch/breaker/06homology_gene2haps/deseq_prep/divgent_fq/stringtie/DEseq.adj.gene.list.inhap0")
	dea<-as.matrix(dea)

	#R直出结果
	hgt_go_enrich<-enricher(hgt, TERM2GENE = term2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff = 0.05)
	dea_go_enrich<-enricher(dea, TERM2GENE = term2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff = 0.05)
	
	new_cols<-c("Pathway","GID")
	gene2pathway<- gene2pathway[,new_cols]
	hgt_kegg_enrich<-enricher(hgt, TERM2GENE = gene2pathway, TERM2NAME = pathway2name, pvalueCutoff = 1, qvalueCutoff = 0.05)
	dea_kegg_enrich<-enricher(dea, TERM2GENE = gene2pathway, TERM2NAME = pathway2name, pvalueCutoff = 1, qvalueCutoff = 0.05)
	
	pdf('results_dea_protein.pdf',width=10)
	dotplot(dea_kegg_enrich, showCategory=30) + ggtitle("Dotplot for dea KEGG in protein database")
	dotplot(dea_go_enrich, showCategory=30) + ggtitle("dotplot for dea GO in protein database")
	dev.off()
	
	pdf('results_hgt_protein.pdf',width=10)
        dotplot(hgt_kegg_enrich, showCategory=30) + ggtitle("Dotplot for hgt KEGG in protein database")
        dotplot(hgt_go_enrich, showCategory=30) + ggtitle("Dotplot for hgt GO in protein database")
        dev.off()
	
	hgt_dea<-read.delim("/home/shangao/Scratch/breaker/06homology_gene2haps/deseq_prep/divgent_fq/stringtie/hgt2Dea/dea2hgt.list")
	hgt_dea<-as.matrix(hgt_dea)
	hgt_dea_go_enrich<-enricher(hgt_dea, TERM2GENE = term2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff = 0.05)
	hgt_dea_kegg_enrich<-enricher(hgt_dea, TERM2GENE = gene2pathway, TERM2NAME = pathway2name, pvalueCutoff = 1, qvalueCutoff = 0.05)
	pdf('results_both_dea_hgt_protein.pdf',width=10)
	dotplot(hgt_dea_kegg_enrich, showCategory=30) + ggtitle("Dotplot for both dea and hgt KEGG in protein database")
	dotplot(hgt_dea_go_enrich, showCategory=30) + ggtitle("dotplot for both dea and hgt GO in protein database")
	dev.off()

### 9.HGT  https://github.com/reubwn/hgt
	
	conda activate BRAKER
	
#### 9.1 get pro.seq from genome by gff

	gffread ../01braker/braker/braker.gtf -g /home/shangao/Scratch/hi-c/sala2/Ppr2/scaffolds/scaffolds_FINAL.fasta -y Ppr.pep

#### 9.2 prepare the db and tax files

	bd: /home/shangao/Data/HGT_db/
	perl -lane 'if(/^>(\w+)\s.+TaxID\=(\d+)/){print "$1 $2"}' <(zcat uniref90.fasta.gz) | gzip > uniref90.fasta.taxlist.gz
	#注意版本，要用旧版
	diamond blastp --sensitive --index-chunks 1 -k 500 -e 1e-5 -p 32 -q Ppr.pep -d /home/shangao/Data/HGT_db/uniref90.dmnd -a Ppr_diamond_results
	cat <(zcat /home/shangao/Data/HGT_db/uniref90.fasta.taxlist.gz) <(diamond view -a Ppr_diamond_results.daa)|perl -lane 'if(@F==2){$tax{$F[0]}=$F[1];}else{if(exists($tax{$F[1]})){print join("\t",@F,$tax{$F[1]});} else {print join("\t",@F,"NA");}}'| gzip > diamond_results.daa.taxid.gz
	
#### 9.3 get results

	/home/shangao/script/hgt/diamond_to_HGT_candidates.pl -i diamond_results.daa.taxid.gz -f Ppr.pep -p /home/shangao/Data/HGT_db/taxdump
	
The program outputs 3 files, suffixed with the tags:

  HGT_results: hU, AI and CHS scores for all query proteins.
  HGT_candidates: All queries which show evidence from both hU (or AI) and CHS above the specifed thresholds (defaults are >=30 for hU, >=90% for CHS).
  HGT_warnings: Any non-fatal warnings detected during the run. Worth checking, but most warnings can probably safely be ignored.
  
 when you meet some perl lib problem get out of conda
 
 #### 9.4 HGT_candidates_to_fasta
 
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
	
	/NVME/Software/TSEBRA/bin/rename_gtf.py --gtf braker1+2_combined.gtf --prefix hap1 --out hap1.gtf
	
## 10. RIdeogram syteny
### 1.prepare assembly
	seqkit sort -lr output_scaffolds100.fasta> output_scaffolds100.sort.fasta
	seqkit fx2tab output_scaffolds100.sort.fasta -l -g -n -i -H

	cat change.sh 
	for i in $(cat order.list)
	do
	echo "sed -i 's/$i//g' output_scaffolds100.sort.fasta" >> changeOrder2.sh
	done
	
	#change name to the order beginning with the biggest one
	
	minimap2 -t 40 -x asm5 /RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.fa output_scaffolds100.sort.fasta >minimap_hap1_haploid0_scaffolds100.paf
	
	minimap2 -t 40 -x asm5 ~/Data/hifiasm_tell-seq/Ppr/hap1/test1/output_scaffolds100.sort.fasta output_scaffolds.sort.fasta > hap2_v_hap1.paf

	awk -v OFS='\t' '{print $1,"1",$2,"969696","hap2","12","252525"}' hap2_v_hap1.paf |sort|uniq >/home/shangao/Scratch/breaker/05RIdeogram_plot/synteny/hap1tohap2/hap2.karyotype

	awk -v OFS='\t' '{print $6,"1",$7,"969696","hap1","12","252525"}' hap2_v_hap1.paf |sort|uniq > /home/shangao/Scratch/breaker/05RIdeogram_plot/synteny/hap1tohap2/hap1.karyotype

	awk -v OFS='\t' '{print $1,$3,$4,$6,$8,$9,"cccccc"}' hap2_v_hap1.paf > /home/shangao/Scratch/breaker/05RIdeogram_plot/synteny/hap1tohap2/synteny
	awk '$3-$2>=100000 {print$0}' synteny > synteny.f



	cat karyotype1 ../hap1/hap1.karyotype.sort  > karyotype
	
	awk NF ../hap2/new_syteny|awk -v OFS="\t" '{print$0,"1"}' > 1
	
	awk -v OFS="\t" '{print$4,$5,$6,$1,$2,$3,$7,"2"}' ../hap1/synteny_h1_h0.color >2
	
	awk -v OFS="\t" '{print$0,"3"}' ../hap1tohap2/synteny.f.color > 3

	cat 1 2 3 > 4

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

## 11. circos mcscanx
	gffread /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap1/braker1+2_combined.gtf -g /RAID/Data/mites/genomes/Ppr/version03/hap1/hap1.fa -y Ppr_hap1.pep

	#gffread /RAID/Data/mites/genomes/Ppr/version03/Ppr.gtf -g /RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.fa -y Ppr.pep
	
	awk -v OFS='\t' '$3=="gene"{print$1,$9,$4,$5}' /RAID/Data/mites/genomes/Ppr/version03/Ppr.gtf > Ppr_hap0.gff
	
	diamond makedb --in Ppr.pep -d Ppr.pep.aa
	
	diamond blastp -e 1e-5 -p 8 -q Ppr_hap1.pep -d Ppr.pep -a Ppr_hap1tohap0.blastp
	
	diamond view -a Ppr_hap1tohap0.blastp -o Xyz_hapi1.blast
	
	/home/jbast/anaconda3/envs/BRAKER/bin/getAnnoFasta.pl --seqfile=/RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.fa /RAID/Data/mites/genomes/Ppr/version03/Ppr.gtf
	
	perl -lane 'print join("\t",$F[0],$F[8],$F[3],$F[4]) if ($F[2]eq"transcript")' ~/Scratch/breaker/04hgt/Ppr/instagraal/TSEBRA/Ppr.gff3 > Xyz.gff
	### get simple gff and blast results then use mcscanx get results
	https://genomevis.usask.ca/
	/home/shangao/Scratch/breaker/02hap2assembly
	
## 12. SV
	conda activate leviathan
	
	cat Data/gaoshan/hifiasm_tell-seq/tell_sort/leviathan/test.sh
	for i in Ppr_T501 \
	Ppr_T502 \
	Ppr_T505 \
	Ppr_T506 
	do
	LRez index bam -p -b ../tell_sort_out/$i/${i}_temp/$i.sorted.bam -o ../tell_sort_out/$i/${i}_temp/$i.barcodeIndex.bci
	LEVIATHAN -b ../tell_sort_out/$i/${i}_temp/$i.sorted.bam -i ../tell_sort_out/$i/${i}_temp/$i.barcodeIndex.bci -g /RAID/Data/mites/genomes/Ppr/version03/Ppr_instagrall.polished.FINAL.fa -o $i.SV.vcf
	done

## 13. KAKS
### get representive transcription

	####cat rename1.sh
	agat_convert_sp_gff2gtf.pl --gff ~/Scratch/breaker/09maker/$i/EVM.all.raw.gff -o ~/Scratch/breaker/09maker/$i/EVM.all.gtf
	/NVME/Software/TSEBRA/bin/rename_gtf.py --gtf ~/Scratch/breaker/09maker/$i/EVM.all.gtf --prefix $i --out ~/Scratch/breaker/09maker/$i/$i.gtf
	agat_convert_sp_gxf2gxf.pl -g ~/Scratch/breaker/09maker/$i/$i.gtf -o ~/Scratch/breaker/09maker/$i/$i.gff3
	agat_convert_sp_gff2gtf.pl --gff ~/Scratch/breaker/09maker/$i/$i.gff3 -o ~/Scratch/breaker/09maker/$i/$i.gtf
	mv ~/Scratch/breaker/09maker/$i/EVM.all.gff ~/Scratch/breaker/09maker/$i/EVM.all.raw.gff
	cp ~/Scratch/breaker/09maker/$i/$i.gff3 ~/Scratch/breaker/09maker/$i/EVM.all.gff3
	
	####cat rename1.sh
	agat_convert_sp_gff2gtf.pl --gff /home/shangao/Scratch/breaker/09maker/$i/add_UTR/$i.gff3 -o /home/shangao/Scratch/breaker/09maker/$i/add_UTR/$i.gtf
	#agat_convert_sp_gxf2gxf.pl -g /home/shangao/Scratch/breaker/09maker/Ppr_alt1/add_UTR/$i.gff3 -o /home/shangao/Scratch/breaker/09maker/Ppr_alt1/add_UTR/$i.gtf
	/NVME/Software/TSEBRA/bin/rename_gtf.py --gtf /home/shangao/Scratch/breaker/09maker/$i/add_UTR/$i.gtf --prefix $i --out /home/shangao/Scratch/breaker/09maker/$i/add_UTR/$i.named.gtf
	agat_convert_sp_gxf2gxf.pl -g /home/shangao/Scratch/breaker/09maker/$i/add_UTR/$i.named.gtf -o /home/shangao/Scratch/breaker/09maker/$i/add_UTR/$i.named.gff3

	
	java -cp /home/shangao/software/TBtools-1.098745/TBtools_JRE1.6.jar biocjava.bioIO.FastX.FastaIndex.FastaLongestRepresentater --inFasta Ppr_alt1.cds --outFasta Ppr_alt1_longest.cds

	/home/shangao/software/Ruminant_pip/Annotation/Gene_annotation/script/cds2aa.pl Ppr_alt1_longest.cds > Ppr_alt1_longest.pep


	java -cp /home/shangao/software/TBtools-1.098745/TBtools_JRE1.6.jar biocjava.bioIO.FastX.FastaIndex.FastaLongestRepresentater --inFasta /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap1/Ppr_hap1.pep --outFasta Ppr_hap1_longest_protein.fa
	java -cp /home/shangao/software/TBtools-1.098745/TBtools_JRE1.6.jar biocjava.bioIO.FastX.FastaIndex.FastaLongestRepresentater --inFasta /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap2/Ppr_hap2.pep --outFasta Ppr_hap2_longest_protein.fa
### run OrthoFinder
	/home/shangao/software/OrthoFinder_source/orthofinder.py -f ./
	
	grep -v ',' hap1_hap2/OrthoFinder/Results_Apr11/Orthogroups/Orthogroups.tsv > hap1_hap2.homo

	awk -vOFS='\t' '{print$6,$7}' ../hap1_vs_hap2_filtered_by_hap0 > hap1_hap2.homo
	
	python /home/shangao/script/python/01dea/05filter_hap1_hap2_hap0_orthfinder.py -s hap0_hap1.homo -t hap0_hap2.homo -l hap1_hap2.homo -o hap1_vs_hap2_filtered_by_hap0

### run paraAT
	ParaAT.pl -h hap1_vs_hap2.homologs -n hap1_vs_hap2.cds -a hap1_vs_hap2.pep -p proc -m muscle -f axt -g -k -o hap1_vs_hap2_results
	awk -v OFS='-' '{print$6,$7}' hap1_vs_hap2_withoutUTR_longest.filtered_by_hap0 > kaks_filtered_name
### get Kaks
	for i in $(cat kaks_filtered_name)
	do
	cp ../hap1_hap2_kaks/hap1_vs_hap2_longest.results/${i}.cds_aln.axt.kaks ./kaks
	done

	cat kaks/* > hap1_hap2_filtered.kaks
	grep -v 'Seq' hap1_hap2_filtered.kaks|awk -vOFS='\t' '$4>0.001&&$5!="NA"{print$1,$5,$6}'> hap1_hap2_rmNA_ks0.001.kaks
	
#### dot plot
	#awk -v OFS='\t' '{print$1,$5,$6}' hap1_hap2_filtered_rmNAkaksless10.kaks > ../plot/dot/hap1_hap2_filtered_rmNAkaksless10.kaks.simple
	
	pdf("dot.pdf")
	plot(data$Ka.Ks, -log(data$P.Value.Fisher.),xlab = "Ka/Ks",ylab = "-log10(P)", type="p",xlim=c(0,10),pch = 20,cex = 0.8,main="Ka/Ks distribution for allelic genes")
	dev.off()
#### volin plot
	python /home/shangao/script/python/01dea/07list2list4dea2kaks.py -s /home/shangao/Scratch/breaker/07kaks/without_UTR/hap1_vs_hap2_withoutUTR_longest.filtered_by_hap0 -t /home/shangao/Scratch/breaker/06homology_gene2haps/deseq_prep/divgent_fq/featureCount/without_UTR/DEseq.adj.gene -l hap1_hap2_filtered_rmNAkaksless10.kaks.simple -o hap1_hap2_filtered_rmNAkaksless10.kaks.simple.withdea
	library(ggpubr)
	library(ggplot2)
	kaks_dealess3 <-read.table('hap1_hap2_filtered_rmNAkaksless10.kaks.simple.withdea', row.names='Sequence',header=T)
	my_comparisons <-list( c("1", "2"))
	kaks_dealess3$group<-as.factor(kaks_dealess3$group)
	kaks_dealess10<-kaks_dealess3[kaks_dealess3$Ka.Ks<3,]
	pdf("violin_plotless10.pdf")
	ggplot(kaks_dealess10, aes(x = group, y = Ka.Ks))+ geom_violin(trim = FALSE, fill = "steelblue")+geom_boxplot(width = 0.2)+stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y =3)
	dev.off()
	
## 14. build the shared folder
### Hap0
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/softmask/braker_UTR_out/test_TSEBRA/Ppr.unigene.pep ./Ppr.pep
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/softmask/braker_UTR_out/test_TSEBRA/cds.emapper.gff ./Ppr_cds_emapper.gff
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/softmask/braker_UTR_out/test_TSEBRA/Ppr.gtf .
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/softmask/braker_UTR_out/test_TSEBRA/Ppr.gff .
### Hap2
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap2/hap2.gff .	
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap2/hap2.gtf .
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap2/cds.emapper.gff ./Ppr_hap2_emapper.gff
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap2/Ppr_hap2_unigene.pep .
### Hap1
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap1/hap1.gff .
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap1/hap1.gtf .
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap1/cds.emapper.gff ./Ppr_hap1_emapper.gff
	cp /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap1/Ppr_hap1_unigene.pep .

## 15. dea

	/home/shangao/software/subread-2.0.3-Linux-x86_64/bin/featureCounts -T 40 -p -B -M -a /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap1/hap1.gtf -o featureCounts_hap1.txt /home/shangao/Scratch/breaker/000rna-seq/Ppr/divergent_fq4hap12/map_indentpent/Ppr_part1Aligned.sortedByCoord.out.bam /home/shangao/Scratch/breaker/000rna-seq/Ppr/divergent_fq4hap12/map_indentpent/Ppr_part2_1Aligned.sortedByCoord.out.bam
	
	python /home/shangao/script/python/dea/01merge_featureCounts_result_tohap1.py -s featureCounts_hap1.txt -t /home/shangao/Scratch/breaker/07kaks/hap1_vs_hap2_unigene/unigene_homologs/hap1_vs_hap2_unigene.homologs.1 -l featureCounts_hap2.txt -o merged_feature_count_hap1_hap2.txt
	
	R
	library(DESeq2)
	data <- read.table("merged_feature_count_hap1_hap2_title.txt", header=TRUE, quote="\t")
	sampleNames <- c("Ppr_part1", "Ppr_part2_1", "Ppr_part2", "Ppr_part1_2")
	database <- data.frame(name=sampleNames, condition=c("hap1", "hap1","hap2", "hap2"))
	names(data)[7:10] <- sampleNames
	countData <- as.matrix(data[7:10])
	rownames(countData) <- data$Geneid
	rownames(database) <- sampleNames
	dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
	dds <- dds[ rowSums(counts(dds)) > 1, ]
	dds <- DESeq(dds)
	res <- results(dds)
	write.csv(res, "res_des_output.csv")
	resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
	write.csv(resdata, "DEseq.result.gene", row.names=FALSE)


	head -n 1 DEseq.result.gene |sed 's/,/\t/g' > gene.title
	sed 's/,/\t/g' DEseq.result.gene |awk '$7<= 0.05 {print$0}' > gene
	cat gene.title gene > DEseq.adj.gene
### new version
	/home/shangao/software/subread-2.0.3-Linux-x86_64/bin/featureCounts -T 40 -p -B -M -a /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap1/without_UTR/hap1.gtf -o featureCounts_hap1.txt /home/shangao/Scratch/breaker/000rna-seq/Ppr/divergent_fq4hap12/map_indentpent/Ppr_part1Aligned.sortedByCoord.out.bam /home/shangao/Scratch/breaker/000rna-seq/Ppr/divergent_fq4hap12/map_indentpent/Ppr_part2_1Aligned.sortedByCoord.out.bam
	/home/shangao/software/subread-2.0.3-Linux-x86_64/bin/featureCounts -T 40 -p -B -M -a /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_hap2/without_UTR/hap2.gtf -o featureCounts_hap2.txt /home/shangao/Scratch/breaker/000rna-seq/Ppr/divergent_fq4hap12/map_indentpent/Ppr_part1_2Aligned.sortedByCoord.out.bam /home/shangao/Scratch/breaker/000rna-seq/Ppr/divergent_fq4hap12/map_indentpent/Ppr_part2Aligned.sortedByCoord.out.bam
	
	cat feature_count_results/featureCounts_hap1.txt feature_count_results/featureCounts_hap2.txt > merged_featureCounts.txt

	python /home/shangao/script/python/01dea/06list2list4hap12tohap0_usingOrthFinder_results.py -s /home/shangao/Scratch/breaker/07kaks/without_UTR/hap1_vs_hap2_withoutUTR_longest.filtered_by_hap0 -t merged_featureCounts.txt -o merged_featureCounts_hap1hap2 -l merged_featureCounts_hap1hap2.missed
	
	library(DESeq2)
	data <- read.table("merged_featureCounts_hap1hap2", header=TRUE, quote="\t")
	sampleNames <- c("Ppr_part1", "Ppr_part2_1", "Ppr_part2", "Ppr_part1_2")
	database <- data.frame(name=sampleNames, condition=c("hap1", "hap1","hap2", "hap2"))
	rownames(countData) <- data$Geneid
	dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
	dds <- dds[ rowSums(counts(dds)) > 1, ]
	dds <- DESeq(dds)
	res <- results(dds)
	write.csv(res, "res_des_output.csv")
	resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
	write.csv(resdata, "DEseq.result.gene", row.names=FALSE)

## 16. MAKER
### 16.1 BRAKER (ab initio)
	braker.pl --cores 40 --species=Ppr_2bam --genome=/RAID/Data/Mites/Genomes/Ppr/German_eiffel/hap0/Ppr.hap0.softmasked.fasta \
	 --softmasking --bam=/home/shangao/Scratch/breaker/000rna-seq/Ppr/Ppr_instagrall/changed_chr/PprAligned.sortedByCoord.out.bam,/home/shangao/Scratch/breaker/000rna-seq/Ppr/Ppr_instagrall/changed_chr/Ppr_downAligned.sortedByCoord.out.bam \
	  --gff3 --useexisting --GENEMARK_PATH=/home/shangao/Software/BrAKER/gmes_linux_64/ --workingdir=braker_UTR_out3
	  
### 16.2 Trinity+PASA
#### Trinity
	LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
	export LD_LIBRARY_PATH
	Trinity --seqType fq --max_memory 60G --CPU 20  \
	  	--left /RAID/Data/mites/transcriptomes/RNA_annotation/150360_R1.fq.gz,/RAID/Data/mites/transcriptomes/RNA_annotation/clean_data/SRR4039022_1_val_1.fq.gz \
	  	--right /RAID/Data/mites/transcriptomes/RNA_annotation/150360_R2.fq.gz,/RAID/Data/mites/transcriptomes/RNA_annotation/clean_data/SRR4039022_2_val_2.fq.gz \
		--output ./trinity.out1
		
	
#### PASA (conda activate pasa)
	/home/jbast/anaconda3/envs/funannotate_env/opt/pasa-2.5.2/Launch_PASA_pipeline.pl \
	-c /home/shangao/Scratch/breaker/03RNA_assembly/denovo_assembly/trinity.out/pasa/alignAssembly.config \
	-C -R -g /RAID/Data/Mites/Genomes/Ppr/German_eiffel/hap0/Ppr.hap0.softmasked.fasta \
	-t /home/shangao/Scratch/breaker/03RNA_assembly/denovo_assembly/trinity.out/pasa/../Trinity.fasta \
	--ALIGNERS blat --CPU 10
	
	cat alignAssembly.config
	## templated variables to be replaced exist as <__var_name__>

	# database settings
	DATABASE=/home/shangao/Scratch/breaker/03RNA_assembly/denovo_assembly/trinity.out/pasa/pasa.sqlite
	
	#######################################################
	# Parameters to specify to specific scripts in pipeline
	# create a key = "script_name" + ":" + "parameter" 
	# assign a value as done above.
	
	#script validate_alignments_in_db.dbi
	validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
	validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80
	
	#script subcluster_builder.dbi
	subcluster_builder.dbi:-m=50
	
### 16.3 StringTie
	stringtie -p 50 -o Ppr_hap0_PprAligned_stringtie.gtf /home/shangao/Scratch/breaker/000rna-seq/Ppr/Ppr_instagrall/changed_chr/PprAligned.sortedByCoord.out.bam 
	stringtie -p 50 -o Ppr_hap0_Ppr_downAligned_stringtie.gtf /home/shangao/Scratch/breaker/000rna-seq/Ppr/Ppr_instagrall/changed_chr/Ppr_downAligned.sortedByCoord.out.bam
	stringtie --merge -p 50  -o merged.gtf Ppr_hap0_PprAligned_stringtie.gtf Ppr_hap0_Ppr_downAligned_stringtie.gtf
	/home/jbast/anaconda3/envs/funannotate_env/opt/pasa-2.5.2/misc_utilities/cufflinks_gtf_genome_to_cdna_fasta.pl merged.gtf /RAID/Data/Mites/Genomes/Ppr/German_eiffel/hap0/Ppr.hap0.softmasked.fasta > Ppr_merged.fasta
	/home/jbast/anaconda3/envs/funannotate_env/opt/pasa-2.5.2/misc_utilities/cufflinks_gtf_to_alignment_gff3.pl merged.gtf > transcripts.gff3

	/home/shangao/.conda/envs/Trinity/bin/TransDecoder.LongOrfs -t Ppr_merged.fasta > TransDecoder.LongOrfs.log
	/home/shangao/.conda/envs/Trinity/bin/TransDecoder.Predict -t Ppr_merged.fasta > TransDecoder.Predict.log 

	/home/jbast/anaconda3/envs/funannotate_env/bin/cdna_alignment_orf_to_genome_orf.pl Ppr_merged.fasta.transdecoder.gff3 transcripts.gff3 Ppr_merged.fasta > Ppr_merged.fasta.transdecoder.genome.gff3
	
### 16.4 EVidenceModeler 
	ln /home/shangao/Scratch/breaker/03RNA_assembly/denovo_assembly/trinity.out/pasa/pasa.sqlite.pasa_assemblies.gff3 transcript_alignments.gff3
	cat /home/shangao/Scratch/breaker/03RNA_assembly/genome_guide/Ppr/stringtie/hap0/Ppr_merged.fasta.transdecoder.genome.gff3 /home/shangao/Scratch/breaker/01braker/Ppr/Ppr_instagrall/softmask/braker_UTR_out3/augustus.hints.gff3 > gene_predictions.gff3
	
	export genome=/RAID/Data/Mites/Genomes/Ppr/German_eiffel/hap0/Ppr.hap0.softmasked.fasta
	export prefix=Ppr_hap0
	
	/home/shangao/software/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome $genome --gene_predictions `pwd`/gene_predictions.gff3 --transcript_alignments `pwd`/transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
	
	/home/shangao/software/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl \
		--genome $genome --weights /home/shangao/Scratch/breaker/09maker/hap0/single_bam/weights.txt \
		--gene_predictions `pwd`/gene_predictions.gff3 \
		--transcript_alignments `pwd`/transcript_alignments.gff3 \
		--output_file_name evm.out  --partitions partitions_list.out >  commands.list
	parallel --jobs 10 < commands.list
	
	/home/shangao/software/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
	
	/home/shangao/software/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $genome
	
	find . -regex ".*evm.out.gff3" -exec cat {} \; | bedtools sort -i - > EVM.all.gff
	
### 16.5 fix GFF and convert GTF, pep,cds
	export genome=/RAID/Data/Mites/Genomes/Ppr/German_eiffel/hap0/Ppr.hap0.softmasked.fasta
	export prefix=Ppr_hap0
	
	/home/jbast/anaconda3/bin/gffread EVM.all.gff -g $genome -y EVM.all.pep
	
	agat_convert_sp_gff2gtf.pl --gff EVM.all.gff -o EVM.all.gtf
	
	/NVME/Software/TSEBRA/bin/rename_gtf.py --gtf EVM.all.gtf --prefix $prefix --out $prefix.gtf 
	
	agat_convert_sp_gxf2gxf.pl -g $prefix.gtf -o $prefix.gff
	
	awk -v OFS="\t" '$5-$4>3 {print$0}' $prefix.gff > $prefix.gff.1
	
	/home/jbast/anaconda3/bin/gffread $prefix.gff.1 -g $genome -y $prefix.$prefix.gff.1.pep
	
	/home/jkirangw/miniconda3/envs/genome_asm/bin/bioawk -c fastx 'length($seq) < 50 {print $name}' $prefix.$prefix.gff.1.pep | cut -d '=' -f 2 > short_aa_gene_list.txt
	
	grep -v -w -f short_aa_gene_list.txt $prefix.gff.1 > ${prefix}_filter_50aa_exon.gff
	
	rm -f $prefix.$prefix.gff.1.pep $prefix.gff.1 
	
	/home/jbast/anaconda3/bin/gffread ${prefix}_filter_50aa_exon.gff -g $genome -y ${prefix}_filter_50aa_exon.pep
	/home/jbast/anaconda3/bin/gffread ${prefix}_filter_50aa_exon.gff -g $genome -x ${prefix}_filter_50aa_exon.cds
	agat_convert_sp_gff2gtf.pl --gff ${prefix}_filter_50aa_exon.gff -o ${prefix}_filter_50aa_exon.gtf	
	
	
	rename new_version
	agat_convert_sp_gff2gtf.pl --gff ~/Scratch/breaker/09maker/$i/EVM.all.gff -o ~/Scratch/breaker/09maker/$i/EVM.all.gtf
	/NVME/Software/TSEBRA/bin/rename_gtf.py --gtf ~/Scratch/breaker/09maker/$i/EVM.all.gtf --prefix $i --out ~/Scratch/breaker/09maker/$i/$i.gtf
	agat_convert_sp_gxf2gxf.pl -g ~/Scratch/breaker/09maker/$i/$i.gtf -o ~/Scratch/breaker/09maker/$i/$i.gff3
	agat_convert_sp_gff2gtf.pl --gff ~/Scratch/breaker/09maker/$i/$i.gff3 -o ~/Scratch/breaker/09maker/$i/$i.gtf
	mv ~/Scratch/breaker/09maker/$i/EVM.all.gff ~/Scratch/breaker/09maker/$i/EVM.all.raw.gff
	cp ~/Scratch/breaker/09maker/$i/$i.gff3 ~/Scratch/breaker/09maker/$i/EVM.all.gff3
	
### 16.6 Add UTR(run rename.sh before)
	for i in Russia_hap0 \
	Russia_hapA \
	Russia_hapB
	do
	mkdir $i/add_UTR
	echo "DATABASE=/home/shangao/Scratch/breaker/09maker/$i/pasa/test.sqlite" > /home/shangao/Scratch/breaker/09maker/$i/add_UTR/annotCompare.config
	echo '
	genome=
	/home/jbast/anaconda3/pkgs/pasa-2.5.2-h87f3376_0/opt/pasa-2.5.2/scripts/Load_Current_Gene_Annotations.dbi -c  annotCompare.config -g $genome -P ../EVM.all.gff
	
	/home/jbast/anaconda3/envs/funannotate_env/opt/pasa-2.5.2/Launch_PASA_pipeline.pl -c  annotCompare.config -A -g $genome -t ~/Scratch/breaker/03RNA_assembly/denovo_assembly/trinity.out/pasa/Trinity.fasta.clean
	' > /home/shangao/Scratch/breaker/09maker/$i/add_UTR/test.sh
	
	done
## 17 effect size
	https://rpubs.com/cly/compare2condition_independent
	conda activate R2.11
 	R
  	library(dplyr)      #data_frame, %>%, filter, summarise, group_by
	library(rstatix)    # cohen's d
	library(ggplot2)    # ggplot, stat_…, geom_…, etc
	library(coin) 
 	kaks_dealess3 <-read.table('exon.list')
  	head(kaks_dealess3)
	    V1 V2  V3
	1 4.25  3 dea
	2 2.00  3 dea
 	exon_SNP <-read.table('exon.list')
	wilcox_effsize(V1 ~ V3, data = exon_SNP)
  	t_effectSize <- cohens_d(V1 ~ V3, data = exon_SNP)


