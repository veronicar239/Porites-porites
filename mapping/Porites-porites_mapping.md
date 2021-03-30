---
layout: post
title: Porites porites transcriptome data processing
date: '2020-08-28'
categories: Protocols
tags: [RNASeq, Bioinformatics]
---

# Coral transcriptome data processing
## ***Porites porites***

### Analyzed by Veronica Radice, Barshis Lab, Old Dominion University

### 2-year transplant of ***Porites porites***
- Puerto Morelos, Mexico
- from the reef to low pH ojo site (submarine discharge spring) and control site
- Ana Martinez (PhD candidate) & Adina Paytan (PI) & Dan Barshis (PI)

### Overall transplant design (information from Ana)
- Corals were taken from 3 different origins: within the ojo (center), outside the ojo (control) and in the reef.
- Corals were then fixed in 3 transplant sites: within the ojo (Ojo LAJA Center and Ojo NORTE Center) and in a control site (Ojo Laja CONTROL).
- From each colony, 3 replicates (or individual cores) were taken

- This species ***Porites porites*** is only found on the reef, not in the ojo lagoons (Martinez et al. 2019 Proc R Soc B)
- 10 samples total
    - 7 samples R-C (Reef Origin transplanted to Control site)
        - 80% survival during transplantation (Martinez et al. 2019 Proc R Soc B)
    - 3 samples R-La (Reef Origin transplanted to Ojo Laja site)
        - 33% survival during transplantation
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

----------------------------------------------------------------------------------------------
## Data
- raw data (year 2 of transplant, fastq g-zipped) backed up on RC drive
- contains data for multiple species:  *Porites astreoides*, *Porites porites*, *Siderastrea siderea*
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2017_April/gslserver.qb3.berkeley.edu/170419_50SR_HS4K2A/Paytan

- other raw data from Paytan (year 1 data?) backed up on RC drive:
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2014_October/AdinaOA_2014_October_rawdata.tar.gz

### Files:
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites

----------------------------------------------------------------------------------------------
### Copy ***Porites porites*** files
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq/

```
nano cp-P.porites-nof.sh
```

```
#!/bin/bash -l

#SBATCH -o cp-P.porites-nof.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=cp-P.porites-nof

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq/QCFastqs/nofilter/*_Pp_yr2_R1_clippedtrimmed_nofilter.fastq /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/
```

```
sbatch cp-P.porites-nof.sh
```

--------------------------------------------------------------------------------------------
## QA-QC
- check adapter clipping stats for each file
- run trim clipped filtered stats Python script on stats output file 
    - Schafran_trimstatstable_advbioinf.py
- file:  Sid_trimclipstats_out.txt  
    - local folder on computer
- file includes Siderastrea and Porites porites stats
- I had run trimming-clipping on batch of files that included Siderastrea siderea and Porites porites


--------------------------------------------------------------------------------------------
## Test ***Porites porites*** mapping to ***Porites astreoides*** reference transcriptome
- Currently there is no existing ***Porites porites*** reference transcriptome
- We will try mapping to the ***Porites astreoides*** reference transcriptome
- Test alignment with Kenkel et al. 2014 full Past assembly 
- File:  past.fasta
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/

## Check quality of assemblies
avg_cov_len_fasta_advbioinf.py

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py past.fasta
```

```
The total number of sequences is 30740
The average sequence length is 550
The total number of bases is 16907062
The minimum sequence length is 100
The maximum sequence length is 8171
The N50 is 661
Median Length = 653
contigs < 150bp = 899
contigs >= 500bp = 11640
contigs >= 1000bp = 3274
contigs >= 2000bp = 325
```

--------------------------------------------------------------------------------------------
## ***Porites astreoides*** reference transcriptome

## Make transcriptome mappable
- need bowtie build module
- creates 6 files for mapping

```
enable_lmod
module load bowtie2/2.2.4
```

#### Bowtie for ***Porites porites*** mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/past.fasta P_por
```

## Test Kenkel 2014 ***Porites astreoides*** transcriptome reference

Mapping to reference
```
nano mapreads_Ppor.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Ppor.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Ppor

enable_lmod
module load bowtie2/2.2.4

bowtie2 --local -x P_por -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/R_033_C_Pp_yr2_R1_clippedtrimmed_nofilter.fastq -S R_033_C_Pp_nof_Ppor.sam -k 5\n
bowtie2 --local -x P_por -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/R_045_C_Pp_yr2_R1_clippedtrimmed_nofilter.fastq -S R_045_C_Pp_nof_Ppor.sam -k 5\n
bowtie2 --local -x P_por -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/R_049_La_Pp_yr2_R1_clippedtrimmed_nofilter.fastq -S R_049_La_Pp_nof_Ppor.sam -k 5\n
bowtie2 --local -x P_por -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/R_052_La_Pp_yr2_R1_clippedtrimmed_nofilter.fastq -S R_052_La_Pp_nof_Ppor.sam -k 5\n
bowtie2 --local -x P_por -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/R_057_C_Pp_yr2_R1_clippedtrimmed_nofilter.fastq -S R_057_C_Pp_nof_Ppor.sam -k 5\n
bowtie2 --local -x P_por -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/R_058_La_Pp_yr2_R1_clippedtrimmed_nofilter.fastq -S RR_058_La_Pp_nof_Ppor.sam -k 5\n
```

```
sbatch mapreads_Ppor.sh
```

```
cat bowtie2_Ppor.txt 
```

#### NOT FILTERED
```
32983430 reads; of these:
  32983430 (100.00%) were unpaired; of these:
    7682241 (23.29%) aligned 0 times
    11945490 (36.22%) aligned exactly 1 time
    13355699 (40.49%) aligned >1 times
76.71% overall alignment rate

37536994 reads; of these:
  37536994 (100.00%) were unpaired; of these:
    6807559 (18.14%) aligned 0 times
    14119380 (37.61%) aligned exactly 1 time
    16610055 (44.25%) aligned >1 times
81.86% overall alignment rate

10103615 reads; of these:
  10103615 (100.00%) were unpaired; of these:
    2984174 (29.54%) aligned 0 times
    3970308 (39.30%) aligned exactly 1 time
    3149133 (31.17%) aligned >1 times
70.46% overall alignment rate

12053447 reads; of these:
  12053447 (100.00%) were unpaired; of these:
    7536244 (62.52%) aligned 0 times
    2765493 (22.94%) aligned exactly 1 time
    1751710 (14.53%) aligned >1 times
37.48% overall alignment rate

34541453 reads; of these:
  34541453 (100.00%) were unpaired; of these:
    7296321 (21.12%) aligned 0 times
    13727296 (39.74%) aligned exactly 1 time
    13517836 (39.14%) aligned >1 times
78.88% overall alignment rate

20866419 reads; of these:
  20866419 (100.00%) were unpaired; of these:
    6264048 (30.02%) aligned 0 times
    8359922 (40.06%) aligned exactly 1 time
    6242449 (29.92%) aligned >1 times
69.98% overall alignment rate
```


--------------------------------------------------------------------------------------------
### Rename fasta file based on # of contigs
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/

```
nano copy-ref.sh
```

```
#!/bin/bash -l

#SBATCH -o copy-ref.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=copy-ref

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/past.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/
```

```
sbatch copy-ref.sh
```

### Rename file based on # of contigs
```
mv past.fasta 30740_past.fasta
```


--------------------------------------------------------------------------------------------
## Add suffix to fasta reference sequence names
addsuffixtofastaseqnames.py
> /cm/shared/courses/dbarshis/barshislab/referenceomes/
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/addsuffixtofastaseqnames.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/

#### Add suffix Ppor to host reference
```
nano addsuffixtofastaseqnames_Ppor2.sh
```

```
#!/bin/bash -l

#SBATCH -o addsuffixtofastaseqnames_Ppor2.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=addsuffixtofastaseqnames_Ppor2

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/addsuffixtofastaseqnames.py Ppor /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/30740_past.fasta
```

```
sbatch addsuffixtofastaseqnames_Ppor2.sh
```

Output:
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/30740_past_suffixed.fasta


--------------------------------------------------------------------------------------------
## Identify symbiont type of Siderastrea siderea
- based on internal transcribed spacer 2 (ITS2) region, a multicopy genetic marker commonly used to analyse Symbiodinium diversity
- test ITS2 fasta files 
- Arif et al. 2014
- [https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869)
- file:  ArifITS2_mec12869-sup-0001-FileS1.txt

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/ArifITS2_mec12869-sup-0001-FileS1.txt

#### Check assembly details
Arif_ITS2
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py ArifITS2_mec12869-sup-0001-FileS1.txt
```

```
The total number of sequences is 433
The average sequence length is 279
The total number of bases is 120855
The minimum sequence length is 152
The maximum sequence length is 376
The N50 is 283
Median Length = 280
contigs < 150bp = 0
contigs >= 500bp = 0
contigs >= 1000bp = 0
contigs >= 2000bp = 0
```

--------------------------------------------------------------------------------------------
## Map all ***Siderastrea siderea*** sequences to Arif ITS2 (all clades) reference
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

##### Make file mappable
Need bowtie build module - creates 6 files for mapping
```
enable_lmod
module load bowtie2/2.2.4
```

Bowtie files for Arif ITS2 clade mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/ArifITS2_mec12869-sup-0001-FileS1.txt Arif_ITS2
```

##### Mapping to Arif_ITS2 reference
```
nano mapreads_Arif_ITS2_Ppor.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Arif_ITS2_Ppor.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Arif_ITS2_Ppor

enable_lmod
module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Arif_ITS2 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Arif_ITS2.sam -k 5\n; done
```

```
sbatch mapreads_Arif_ITS2_Ppor.sh
```


--------------------------------------------------------------------------------------------
## Count expression - all reads mapped to Arif_ITS2 symbiont reference

```
nano countexpression_Arif_ITS2.sh 
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Arif_ITS2_Ppor.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_Arif_ITS2_Ppor

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/*_Arif_ITS2.sam
```

```
sbatch countexpression_Arif_ITS2.sh
```

#### copy output into local folder on computer
scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/match_counts.txt ./

File:  match_counts_Arif_Ppor.txt


--------------------------------------------------------------------------------------------
## Merge all *nof_Arif_ITS2_counts.txt files into one big table 
- first need to add column to each .txt file with unique sample id
- so we can later identify which sample has which ITS2 sequences

Using regular expressions:

for f in file1 file2 file3; do sed -i "s/$/\t$f/" $f; done

- For each file in the list, this will use sed to append to the end of each line a tab and the filename
- Using the -i flag with sed to perform a replacement in-place, overwriting the file
- Perform a substitution with s/PATTERN/REPLACEMENT/. 
- In this example PATTERN is $, the end of the line, and REPLACEMENT is \t (= a TAB), 
- and $f is the filename, from the loop variable. 
- The s/// command is within double-quotes so that the shell can expand variables
- sed is most practical for pattern substitution and saving in-place. 


```
nano append_filename.sh
```

```
#!/bin/bash -l

#SBATCH -o append_filename.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=append_filename

for f in *_Arif_ITS2_counts.txt; do sed -i "s/$/\t$f/" $f; done
```

```
sbatch append_filename.sh
```

#### Concatenate files
```
cat *_nof_Arif_ITS2_counts.txt > merged_Ppor_ITS2_counts.txt
```

scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/merged_Ppor_ITS2_counts.txt ./

## Outcome symbiont clade mapping
- 82% reads mapped to I sequences
  - BLAST confirmed 100% match to former Symbiodinium sp. partial rRNA genes and ITS2
- After removing mapping to I sequences
  - 75% mapped to C
  - 21% mapped to A
- Need to test the few different C and A references
- Then pick the two with the best mapping rate, merge them 
- Then re-map to see how much we lose to multiply mapping contigs

### Notes from the literature
- ITS2 type C10 reported from *Porites porites* from Bahamas (LaJeunesse 2002)


--------------------------------------------------------------------------------------------
## C1 *Cladocopium goreaui* - Davies transcriptome
- from Assistant Professor Sarah Davies, Boston University [https://sites.bu.edu/davieslab/data-code/](https://sites.bu.edu/davieslab/data-code/)
- Assembled and annotated transcriptome for the symbiotic dinoflagellate algae Cladocopium hosted by Siderastrea siderea with all host contamination removed. 
- Data were generated using Illumina HiSeq2000 2*100bp reads and assembled using Trinity.
- Citaion: Davies SW, Ries JB, Marchetti A, Castillo KD (2018). Symbiodinium functional diversity in the coral Siderastrea siderea is influenced by thermal stress and reef environment, but not ocean acidification. Frontiers in Marine Science. FMARS-05-00150.
- File:  davies_cladeC_feb.fasta
[https://doi.org/10.3389/fmars.2018.00150](https://doi.org/10.3389/fmars.2018.00150)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py davies_cladeC_feb.fasta
```

```
The total number of sequences is 65838
The average sequence length is 1482
The total number of bases is 97581498
The minimum sequence length is 500
The maximum sequence length is 18168
The N50 is 1746
Median Length = 1297
contigs < 150bp = 0
contigs >= 500bp = 65838
contigs >= 1000bp = 40840
contigs >= 2000bp = 13246
```

### Rename fasta file based on # of contigs
```
mv davies_cladeC_feb.fasta 65838_davies_cladeC_feb.fasta
```


--------------------------------------------------------------------------------------------
## Davies Cladocopium transcriptome
### Make transcriptome mappable
- need bowtie build module
- creates 6 files for mapping

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

```
enable_lmod
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb.fasta davies_cladeC
```

### Test Davies Cladocopium reference

Mapping to reference
```
nano mapreads_davies_cladeC.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_davies_cladeC.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_davies_cladeC

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x davies_cladeC -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_davies_cladeC.sam -k 5\n; done
```

```
sbatch mapreads_davies_cladeC.sh
```

```
cat bowtie2_davies_cladeC.txt 
```

Outcome
- 1.8-6.1% overall alignment (avg. 4%)
- 1.6-4.3%% singly aligned (avg. 3.1%)

## Count expression 
```
nano countexpression_davies_cladeC.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_davies_cladeC.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_davies_cladeC

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/*_nof_davies_cladeC.sam
```

```
sbatch countexpression_davies_cladeC.sh
```


--------------------------------------------------------------------------------------------
## ***Cladocopium sp. C92*** reference
- Chen et al. 2019 "Revised genome sequences and annotations of six Symbiodiniaceae taxa"
- file:  Cladocopium_sp_C92.CDS.fna
[https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947](https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947)
[https://espace.library.uq.edu.au/view/UQ:8279c9a](https://espace.library.uq.edu.au/view/UQ:8279c9a)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/Cladocopium_sp_C92/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Cladocopium_sp_C92.CDS.fna
```

```
The total number of sequences is 33421
The average sequence length is 2201
The total number of bases is 73570391
The minimum sequence length is 93
The maximum sequence length is 28947
The N50 is 3198
Median Length = 2229
contigs < 150bp = 1
contigs >= 500bp = 30177
contigs >= 1000bp = 23167
contigs >= 2000bp = 12664
```

#### Rename fasta file based on # of contigs
```
mv Cladocopium_sp_C92.CDS.fna 33421_Cladocopium_sp_C92.CDS.fna
```

### Make transcriptome mappable
- need bowtie build module
- creates 6 files for mapping

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/33421_Cladocopium_sp_C92.CDS.fna Cladocopium_sp_C92
```

### Test Cladocopium_sp_C92 reference

Mapping to reference
```
nano mapreads_Cladocopium_sp_C92.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Cladocopium_sp_C92.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Cladocopium_sp_C92

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Cladocopium_sp_C92 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Cladocopium_sp_C92.sam -k 5\n; done
```

```
sbatch mapreads_Cladocopium_sp_C92.sh
```

```
cat bowtie2_Cladocopium_sp_C92.txt 
```

Outcome
- 0.1-1.7% overall alignment (avg. 0.7%)
- 0.1-1.5% singly aligned (avg. 0.6%)

## Count expression 
```
nano countexpression_Cladocopium_sp_C92.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Cladocopium_sp_C92.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_Cladocopium_sp_C92

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/*_nof_Cladocopium_sp_C92.sam
```

```
sbatch countexpression_Cladocopium_sp_C92.sh
```

--------------------------------------------------------------------------------------------
## ***Cladocopium goreaui*** reference
- Chen et al. 2019 "Revised genome sequences and annotations of six Symbiodiniaceae taxa"
- file:  Cladocopium_goreaui.CDS.fna
[https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947](https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947)
[https://espace.library.uq.edu.au/view/UQ:8279c9a](https://espace.library.uq.edu.au/view/UQ:8279c9a)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/Cladocopium_goreaui/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Cladocopium_goreaui.CDS.fna
```

```
The total number of sequences is 39006
The average sequence length is 1625
The total number of bases is 63419290
The minimum sequence length is 111
The maximum sequence length is 40317
The N50 is 2388
Median Length = 852
contigs < 150bp = 5
contigs >= 500bp = 31703
contigs >= 1000bp = 21315
contigs >= 2000bp = 9624
```

#### Rename fasta file based on # of contigs
```
mv Cladocopium_goreaui.CDS.fna 39006_Cladocopium_goreaui.CDS.fna
```

### Make transcriptome mappable
- need bowtie build module
- creates 6 files for mapping

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/
```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/39006_Cladocopium_goreaui.CDS.fna Cladocopium_goreaui
```

## Test Cladocopium_goreaui reference

Mapping to reference
```
nano mapreads_Cladocopium_goreaui.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Cladocopium_goreaui.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Cladocopium_goreaui

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Cladocopium_goreaui -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Cladocopium_goreaui.sam -k 5\n; done
```

```
sbatch mapreads_Cladocopium_goreaui.sh
```

```
cat bowtie2_Cladocopium_goreaui.txt 
```

Outcome
- 0.1-2.1% overall alignment (avg. 0.8%)
- 0.1-1.9% singly aligned (avg. 0.7%)

## Count expression 
```
nano countexpression_Cladocopium_goreaui.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Cladocopium_goreaui.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_Cladocopium_goreaui

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/*_nof_Cladocopium_goreaui.sam
```

```
sbatch countexpression_Cladocopium_goreaui.sh
```

--------------------------------------------------------------------------------------------
## Symbiont reference selected
- The Davies Cladocopium symbiont reference was selected among the 3 transcriptomes evaluated
  - had ~8 times the amount of singly aligned reads comapred to the other two references
  - had 5% quality aligned (to the ***holobiont*** sequences) compared to 1% in the other two references


--------------------------------------------------------------------------------------------
## Add suffix to fasta reference sequence names
addsuffixtofastaseqnames.py
> /cm/shared/courses/dbarshis/barshislab/referenceomes/
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/

#### Add suffix Sym to symbiont reference
```
nano addsuffixtofastaseqnames_Sym.sh
```

```
#!/bin/bash -l

#SBATCH -o addsuffixtofastaseqnames_Sym.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=addsuffixtofastaseqnames

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/addsuffixtofastaseqnames.py Sym /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/65838_davies_cladeC_feb.fasta
```

```
sbatch addsuffixtofastaseqnames_Sym.sh
```

Output:
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/65838_davies_cladeC_feb_suffixed.fasta


--------------------------------------------------------------------------------------------
## Hybrid reference transcriptome
Concatenate symbiont transcriptome with host transcriptome to make mapping file
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/hybridref/

```
cat /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/30740_past_suffixed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/65838_davies_cladeC_feb_suffixed.fasta > hybridreference_Past_cladeC.fasta
```

```
ls -alh
```

```
-rwxrwxrwx 1 vradice users 113M Sep 15 11:10 hybridreference_Past_cladeC.fasta
```

--------------------------------------------------------------------------------------------
## Map all samples to concatenated (host plus symbiont) hybrid reference transcriptome
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

#### Make files mappable

Create mapping files
```
enable_lmod
module load container_env trinity
crun bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/hybridref/hybridreference_Past_cladeC.fasta hybridref
```

#### Mapping to hybridreference
```
nano mapreads_hybridref.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_hybridref.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_hybridref

enable_lmod
module load container_env trinity

for i in *_clippedtrimmed_nofilter.fastq; do
        crun bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
                     --rg SM:${i%_clippedtrimmed_nofilter.fastq} \
                     --local -x hybridref -U $i \
                        > ${i%_clippedtrimmed_nofilter.fastq}_nof_hybridref.sam -k 5\n;
done
```

```
sbatch mapreads_hybridref.sh
```

Output:
```
cat bowtie2_hybridref.txt
```

Number of singly aligned reads:
```
grep "aligned exactly 1 time" bowtie2_hybridref.txt
```


--------------------------------------------------------------------------------------------
## seq2iso tables for host and symbiont
- Add Carly Kenkel's past_seq2iso.tab as a -g argument to the countxpression_SB_advbioinf.py script
- Sandrine (SB) re-wrote the script to be able to merge counts into isogroups at the same time as counting if you supply a seq2iso table

#### Use associated seq2iso table (host)
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/pastreoides_may2014

Show count of how many lines in table
```
wc -l past_seq2iso.tab 
grep "" -c past_seq2iso.tab 
```
30740 lines in past_seq2iso.tab

#### Add "_Ppor" suffix to contig names (column 1) and gene names (column 2) and change seq2iso table to tab delimited

In Visual Studio, open file past_seq2iso.tab

Use regular expressions:
> Find:
>
> (\w+) (\w+)
> 
> Replace:
>
> $1_Ppor\t$2_Ppor
>
> (Visual Studio Code uses .net syntax)
> 
> Unix syntax equivalent:
>
> \1_Ppor\t\2_Ppor

New file:  past_seq2iso_suffixed_Ppor.tab

scp past_seq2iso_suffixed_Ppor.tab vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/


#### Add suffix to symbiont seq2iso table
- Symbiont seq2iso file (associated with symbiont reference):  davies_cladeC_seq2iso.tab

scp davies_cladeC_seq2iso.tab vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/

Show count of how many lines in table
```
wc -l davies_cladeC_seq2iso.tab
grep "" -c davies_cladeC_seq2iso.tab
```
65838 lines in davies_cladeC_seq2iso.tab

#### Add "_Sym" suffix to contig names (column 1) and gene names (column 2)

Use regular expressions:
> Find:
>
> (\w+)\t(\w+)
> 
> Replace:
>
> $1_Sym\t$2_Sym
>
> (Visual Studio Code uses .net syntax)
> 
> Unix syntax equivalent:
>
> \1_Sym\t\2_Sym

New file:  davies_cladeC_seq2iso_suffixed.tab

scp davies_cladeC_seq2iso_suffixed.tab vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/


--------------------------------------------------------------------------------------------
## Concatenate host and symbiont seq2iso tables
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/hybridref/

```
cat /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/past_seq2iso_suffixed_Ppor.tab /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/davies_cladeC_seq2iso_suffixed.tab > hybrid_Past-cladeC_seq2iso.tab
```

--------------------------------------------------------------------------------------------
## Count expression 
Generate counts file for all samples
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/mapped_hybridref/

```
nano countexpression_hybridref.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_hybridref.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_hybridref

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/*_nof_hybridref.sam -g /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/refassembly/hybridref/hybrid_Past-cladeC_seq2iso.tab
```

```
sbatch countexpression_hybridref.sh
```

```
head -5 match_counts.txt
```

scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/mapped_hybridref/match_counts.txt ./


--------------------------------------------------------------------------------------------
## Parse expression to table

Copy file hybrid_Past-cladeC_seq2iso.tab to local machine.

#### Make genelist

Format should be:

GeneName |
gene1name |
gene2name |
gene3name |
gene4name |

##### Use regular expressions:

Open file hybrid_Past-cladeC_seq2iso.tab in Visual Studio:

> Find:
>
> (\w+)\t(\w+)
> 
> Replace:
>
> $2
>
> (Visual Studio Code uses .net syntax)
> 
> Unix syntax equivalent:
>
> \2

**Add column header "GeneName"**

New file:   genelist_Ppor_Final.txt

Secure copy from local machine.

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/mapped_hybridref/


Secure copy from local machine.
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

scp genelist_Ppor_Final.txt vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

```
nano ParseExpression_c.sh
```

```
#!/bin/bash -l

#SBATCH -o ParseExpression_c.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=ParseExpression_c

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/ParseExpression2BigTable_advbioinf.py genelist_Ppor_Final.txt merged_hybridref_counts_Ppor_b.txt nomatch *_nof_hybridref_counts.txt
```

```
sbatch ParseExpression_c.sh
```


#### Output 

Mapping details:  match_counts.txt

scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/mapped_hybridref/match_counts.txt ./


Holobiont gene expression:  merged_hybridref_counts_Ppor_b.txt

scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/merged_hybridref_counts_Ppor_b.txt ./


Host gene expression:
```
grep _Ppor merged_hybridref_counts_Ppor_b.txt > merged_hybridref_counts_Ppor_Host.csv
```

scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/merged_hybridref_counts_Ppor_Host.csv ./


Symbiont gene expression:
```
grep _Sym merged_hybridref_counts_Ppor_b.txt > merged_hybridref_counts_Ppor_Sym.csv
```

scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/merged_hybridref_counts_Ppor_Sym.csv ./


--------------------------------------------------------------------------------------------



***checking symbiont typing results with revised Arif ITS2 database***

## Identify dominant symbiont type
- based on internal transcribed spacer 2 (ITS2) region, a multicopy genetic marker commonly used to analyse Symbiodinium diversity
- map samples against ITS2 BLAST database
- Arif et al. 2014 [https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869)
    - Corrigendum, with updated ITS2 database file, published 17 July 2019 [https://onlinelibrary.wiley.com/doi/10.1111/mec.14956](https://onlinelibrary.wiley.com/doi/10.1111/mec.14956)
- File:  mec14956-sup-0001-files1_corrigendum.fasta

scp mec14956-sup-0001-files1_corrigendum.fasta vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

#### Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py mec14956-sup-0001-files1_corrigendum.fasta
```

```
The total number of sequences is 400
The average sequence length is 376
The total number of bases is 150400
The minimum sequence length is 376
The maximum sequence length is 376
The N50 is 376
Median Length = 376
contigs < 150bp = 0
contigs >= 500bp = 0
contigs >= 1000bp = 0
contigs >= 2000bp = 0
```

--------------------------------------------------------------------------------------------
## Map all samples to Arif ITS2 (all clades) reference
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

### Make file mappable
Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

Bowtie files for Arif ITS2 clade mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/mec14956-sup-0001-files1_corrigendum.fasta Arif_ITS2_corrigendum
```

##### Mapping to Arif_ITS2_corrigendum reference

```
nano mapreads_Arif_ITS2_corrigendum.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Arif_ITS2_corrigendum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Arif_ITS2_corrigendum

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Arif_ITS2_corrigendum -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Arif_ITS2_corrigendum.sam -k 5\n; done
```

```
sbatch mapreads_Arif_ITS2_corrigendum.sh
```

--------------------------------------------------------------------------------------------
## Count expression - all reads mapped to Arif_ITS2 symbiont reference

```
nano countexpression_Arif_ITS2_corrigendum.sh 
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Arif_ITS2_corrigendum_Ppor.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpression_Arif_ITS2_corrigendum_Ppor

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/*_Arif_ITS2_corrigendum.sam
```

```
sbatch countexpression_Arif_ITS2_corrigendum.sh
```

#### copy output into local folder on computer
scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/match_counts.txt ./

File:  Pporites_symbiont-reference_match-counts.xlsx


--------------------------------------------------------------------------------------------
## Merge all *nof_Arif_ITS2_counts.txt files into one big table 
- first need to add column to each .txt file with unique sample id
- so we can later identify which sample has which ITS2 sequences

Using regular expressions:

for f in file1 file2 file3; do sed -i "s/$/\t$f/" $f; done

- For each file in the list, this will use sed to append to the end of each line a tab and the filename
- Using the -i flag with sed to perform a replacement in-place, overwriting the file
- Perform a substitution with s/PATTERN/REPLACEMENT/. 
- In this example PATTERN is $, the end of the line, and REPLACEMENT is \t (= a TAB), 
- and $f is the filename, from the loop variable. 
- The s/// command is within double-quotes so that the shell can expand variables
- sed is most practical for pattern substitution and saving in-place. 


```
nano append_filename.sh
```

```
#!/bin/bash -l

#SBATCH -o append_filename.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=append_filename

for f in *_Arif_ITS2_corrigendum_counts.txt; do sed -i "s/$/\t$f/" $f; done
```

```
sbatch append_filename.sh
```

#### Concatenate files
```
cat *_Arif_ITS2_corrigendum_counts.txt > merged_Ppor_ITS2_counts_Arif-corrigendum.txt
```

scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/merged_Ppor_ITS2_counts_Arif-corrigendum.txt ./

## Outcome symbiont clade mapping
- 83% reads mapped to I sequences
  - BLAST confirmed 100% match to former Symbiodinium sp. partial rRNA genes and ITS2
- After removing mapping to I sequences
  - 70% mapped to C
  - 23% mapped to A

### Notes from the literature
- ITS2 type C10 reported from *Porites porites* from Bahamas (LaJeunesse 2002)


--------------------------------------------------------------------------------------------

**resume here**
?? 

- Need to test the few different C and A references
- Then pick the two with the best mapping rate, merge them 
- Then re-map to see how much we lose to multiply mapping contigs



> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_porites/

Host genes from reference:

grep -c _Ppor merged_hybridref_counts_Ppor_b.txt
29422