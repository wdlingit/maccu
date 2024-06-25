## Steps to build a co-expression database from DEE2 datasets

This document contains steps for downloading specified SRS metadata and classification of samples. Steps in this document were done in a Ubuntu 20 server with 128GB memory. For human and mouse data, some steps may take up to more than 500GB memory.

### Processing the metadata file and aggregate the read counts

The [metadata tables made by DEE2](https://dee2.io/metadata/) was used for the initial sample qualification. The following steps were done using Excel.

1. The metadata tables provide QC results of SRR accessions, which rather correspond to technical replicates. SRR's are the basic records in DEE2. To qualify biological replicates, i.e., SRS accessions, which also in the metadata tables, we collected SRS accessions where their corresponding SRR's were all PASS in the QC column.
2. To collect count data associated with SRR's from the DEE2 database, we collected SRR accessions (i) under SRS accessions collected in step 1, (ii) with `experiment_library_strategy` of `RNA-Seq`, and (iii) `experiment_library_selection` with `cDNA`, `RANDOM`, `PolyA`, or `Oligo-dT`.
3. Save SRR-SRS mapping (two columns) into a tab-delimited text file.

Download [the count file](https://dee2.io/mx/) and use the perl oneliner command like the following example to extract and aggregate read counts into biological replicates, i.e., SRS accessions.
```
wdlin@comp04:SOMEWHERE/ath$ head ath_SRR_20240529.txt
DRR008476       DRS007600
DRR008477       DRS007601
DRR008478       DRS007602
DRR016112       DRS014211
DRR016113       DRS014211
DRR016114       DRS014212
DRR016115       DRS014212
DRR016116       DRS014212
DRR021335       DRS030798
DRR021336       DRS030797

wdlin@comp04:SOMEWHERE/ath$ bzip2 -dc athaliana_se.tsv.bz2 | perl -ne 'if($.==1){ open(FILE,"<ath_SRR_20240529.txt"); while($line=<FILE>){ chomp $line; $line=~s/^\s+|\s+$//g; @s=split(/\s+/,$line); $srs{$s[0]}=$s[1] } close FILE } chomp; @t=split; if(exists $srs{$t[0]}){ $samples{$srs{$t[0]}}=1; $hash{$t[1]}{$srs{$t[0]}}+=$t[2] } if(eof STDIN){ print "Symbol"; for $s (sort keys %samples){ print "\t$s" } print "\n"; for $g (sort keys %hash){ print "$g"; for $s (sort keys %samples){ if(exists $hash{$g}{$s}){ print "\t$hash{$g}{$s}" }else{ print "\t0" } } print "\n" } }' > ath_sel20240529.nMatrix.txt
```

Points to be noticed:
1. `athaliana_se.tsv.bz2` is the count file downloaded from the DEE2 database
2. `ath_SRR_20240529.txt` is the SRR-SRS mapping (a two-column tab-delimited text file) with collected SRS accessions in above step 3.
3. The output file `ath_sel20240529.nMatrix.txt` (tab-delimited) is the raw count matrix, with columns for samples and rows for genes.

### Duplicate removal

Some samples (SRS) would be repeatedly submitted to the NCBI SRA database. The following steps were applied for removing duplications.

```
wdlin@comp04:SOMEWHERE/ath$ ../scripts/duplicateDetect.pl ath_sel20240529.nMatrix.txt > ath_sel20240529.nMatrix.dupReport

wdlin@comp04:SOMEWHERE/ath$ head ath_sel20240529.nMatrix.dupReport
Reading matrix
Compute hash
Compare
Report
DUP: ERS1647356 ERS1827231
DUP: ERS1647357 ERS1827232
DUP: ERS1647358 ERS1827235
DUP: ERS1647359 ERS1827236
DUP: SRS1042458 SRS2817876
DUP: SRS1121919 SRS2218890

wdlin@comp04:SOMEWHERE/ath$ head -1 ath_sel20240529.nMatrix.txt | perl -ne 'chomp; s/^\s+|\s+$//g; @t=split; print "$_\n" for @t' | perl -ne 'chomp; if($.==1){ open(FILE,"<ath_sel20240529.nMatrix.dupReport"); while($line=<FILE>){ chomp $line; if($line=~/^DUP/){ @s=split(/\s+/,$line); shift @s; shift @s; for $x (@s){ $duplicate{$x}=1 } }} close FILE; print "$_\tnondup\n" }else{ print "$_\t"; if(exists $duplicate{$_}){ print "0\n" }else{ print "1\n" } }' > ath_sel20240529.dup.txt

wdlin@comp04:SOMEWHERE/ath$ head ath_sel20240529.dup.txt
Symbol  nondup
DRS007600       1
DRS007601       1
DRS007602       1
DRS014211       1
DRS014212       1
DRS030797       1
DRS030798       1
DRS047331       1
DRS047332       1
```

The first command was to use the script `duplicateDetect.pl` (in our `scripts` directory) for identifying duplicate columns (samples). Inside the output file (`ath_sel20240529.nMatrix.dupReport` here), lines started with `DUP:` are for duplicated samples. The perl oneliner was to read the header column from the raw count matrix (`head -1 ath_sel20240529.nMatrix.txt`) and generate a 0-1 matrix (`ath_sel20240529.dup.txt`) based on the duplication report. The 0-1 matrix was for indicating which samples are nonduplicated. For samples reported in the ducplication report, only the first sample from each line was specified as nonduplicated.

The last command for duplication removal was to apply `matrixSelection.pl` (in our `scripts` directory). This script takes at least four parameters:
1. selection matrix: in this case, `ath_sel20240529.dup.txt` is the selection matrix. Note that column headers are treated as selection targets.
2. source matrix: a tab-delimited matrix file, with columns for samples.
3. output prefix: an output filename would be in the form `<output prefix>.<target>`.
4. targets: one or more column headers from the selection matrix could be selected for column selection from the source matrix.
```
wdlin@comp04:SOMEWHERE/ath$ ../scripts/matrixSelection.pl
matrixSelection.pl <selMatrix> <sourceMatrix> <outPrefix> [<selTarget>]+

wdlin@comp04:SOMEWHERE/ath$ ../scripts/matrixSelection.pl ath_sel20240529.dup.txt ath_sel20240529.nMatrix.txt ath_sel20240529.nMatrix.txt nondup
```

In this example, the read count matrix without duplicated samples would be named `ath_sel20240529.nMatrix.txt.nondup`.

### In case no sample classification required

We believe that certain normalization is needed for the co-expression database but not the raw counts. In case that no sample classification requried. The following R commands were adopted for the normalization task using the TMM method (PMID: 20196867).
```
library("edgeR")
x <- read.delim("ath_sel20240529.nMatrix.txt.nondup",row.names="Symbol")
dge <- DGEList(counts=x)
dge <- calcNormFactors(dge)
v <- voom(dge,normalize="none")
write.csv(x=v$E,file="sel20240529.nMatrix.TMM")
quit()
```

You may modify the `write.csv` command or use the following perl oneliner to transfer the output CSV file into a tab-delimited text file. Our java program for co-expression computation accepts only tab-delimited matrix files.

```
wdlin@comp04:SOMEWHERE/ath$ cat sel20240529.nMatrix.TMM | perl -ne 'chomp; @t=split(/,/); $nonFirst=0; for $x (@t){ $x=~s/^"|"$//g; print "\t" if $nonFirst; $nonFirst=1; print "$x" } print "\n"' > sel20240529.nMatrix.TMM.txt
```

### Sample classification part 1, downloading metadata from NCBI

We firstly generate a list of SRS accessions of nonduplicated samples.
```
wdlin@comp04:SOMEWHERE/ath$ head -1 ath_sel20240529.nMatrix.txt.nondup | perl -ne 'chomp; @t=split; shift @t; for $x (@t){ print "$x\n" }' > ath_sel20240529.list
```

The `biosampleRetrieveBySRS.pl` (in our `scripts` directory) was used for retriving metadata from NCBI.
```
wdlin@comp04:SOMEWHERE/ath$ ../scripts/biosampleRetrieveBySRS.pl
biosampleRetrieve.pl <listFile> <outSrsBios> <outXML>
```
Points to be noticed:
1. The [NCBI EDirect utility](https://www.ncbi.nlm.nih.gov/books/NBK179288/) is required for running this script
2. This script will write retrieved metadata XML into `<outXML>` and SRS-BioSample accession pairs into `<outSrsBios>`. You may use the line numbers in `<outSrsBios>` to check numbers of SRS records with successfully retrieved metadata.
3. This script will *append* contents to the two output files, and it will process only SRS accessions not in `<outSrsBios>`. That is, you may simply repeat the same command a few number of times for retrieving metadata for the same list without taking care of the outputs. NOTE: It is possible that the NCBI contains no metadata for some SRS accessions. Just mark those SRS accessions kept being searched for a number of times and check them in the NCBI webpage.
4. This script doesn't support parallel processing. You may apply a command like `split -l 2500 -d ath_sel20240529.list ath_sel20240529.list.` to split the list into smaller lists for parallel processing (surely separate output files for separate input lists). Note that NCBI has some query number restriction per second given an API key. Be sure not to exceed the limitation.
5. Variable `$maxTry` was hard-coded as `3` for the number of re-try an `esearch` command. Modify it if needed.
6. Inside the script, the first two `esearch` commands were used for retrieving the corresponding BioSample accession of an SRS accession. They are our current best practices for retrieving BioSample accessions from SRS accessions. Modify them if needed.
7. The last `esearch` command in the script was to extract metadata of the BioSample accession corresponding to an SRS accession. The reason that we extract metadata form BioSample but not SRA is that the metadata from BioSample is generally more detailed than that from SRA.

Again, it is possible that we might retrieve no metadata for some SRS. So we may apply a similar technique to remove them from the count matrix. (`ath_sel20240529.map` is the (merged) `<outSrsBios>` output file)
```
wdlin@comp04:SOMEWHERE/ath$ cat ath_sel20240529.map | perl -ne 'if($.==1){ print "SRS\tgot\n" } chomp; @t=split; print "$t[0]\t1\n"' > ath_sel20240529.gotMetadata

wdlin@comp04:SOMEWHERE/ath$ ../scripts/matrixSelection.pl ath_sel20240529.gotMetadata ath_sel20240529.nMatrix.txt.nondup ath_sel20240529.nMatrix got
```

### Sample classification part 2, an example of arabidopsis ecotypes

Due to the complexity of human-input metadata, we don't have a completely automatic classification method. Here we present an approximation that used to give us enough number of samples after classification. In this session, we present what we had done on arabidopsis ecotypes.

Suppose that `ath_SRS_20240529.txt` is the metadata XML file that we obtained using the script described in the last session. The following two commands helped us for understanding the diversity inside the metadata.
```
wdlin@comp04:SOMEWHERE/ath$ cat ath_SRS_20240529.txt | perl -ne 'chomp; if(/<Attribute attribute_name="(.+?)"/){ $hash{"$1"}++ } if(eof){ for $k (sort {$hash{$b}<=>$hash{$a}} keys %hash){ print "$k\t$hash{$k}\n" } }' > ath_SRS_20240529.attributes

wdlin@comp04:SOMEWHERE/ath$ head ath_SRS_20240529.attributes
tissue  14761
source_name     13210
genotype        11304
ecotype 10567
age     8852
treatment       7384
geo_loc_name    5283
INSDC status    3126
INSDC first public      3126
INSDC center name       3126

wdlin@comp04:SOMEWHERE/ath$ cat ath_SRS_20240529.attributes | perl -ne 'chomp; @t=split(/\t/); $cmd="cat ath_SRS_20240529.txt | grep \"\\\"$t[0]\\\"\" | uniq | head -30"; print "ATTR: $_\n"; system $cmd' | less
ATTR: tissue    14761
      <Attribute attribute_name="organism part" harmonized_name="tissue" display_name="tissue">flower bud (EFO_0001924)</Attribute>
      <Attribute attribute_name="tissue_type" harmonized_name="tissue" display_name="tissue">seedlings</Attribute>
      <Attribute attribute_name="tissue_type" harmonized_name="tissue" display_name="tissue">root</Attribute>
      <Attribute attribute_name="tissue_type" harmonized_name="tissue" display_name="tissue">seedling</Attribute>
      <Attribute attribute_name="tissue_type" harmonized_name="tissue" display_name="tissue">shoot apex</Attribute>
      <Attribute attribute_name="tissue_type" harmonized_name="tissue" display_name="tissue">Inflorescence meristems and young floral buds</Attribute>
      <Attribute attribute_name="tissue_type" harmonized_name="tissue" display_name="tissue">whole tisues</Attribute>
      <Attribute attribute_name="tissue" harmonized_name="tissue" display_name="tissue">1cm long root tips</Attribute>
      <Attribute attribute_name="tissue_type" harmonized_name="tissue" display_name="tissue">cultured cell line MM2d</Attribute>
(deleted)
```
In the metadata XML file, each BioSample is associated with a number of *attributes* which may have different values. For example, samples may have `tissue` attributes of values `root`, `seedling`, .... The first perl oneliner was to collect all attributes and rank them from the most frequently recorded attribute to the least frequently recorded attribute. In file `ath_SRS_20240529.attributes`, we may find that the `tissue` attribute was ranked first, which *should be* corresponding to tissue information. It was also found that the `ecotype` attribute was ranked fourth and that *should be* corresponding to ecotype information. The last perl oneliner command was to list first few nonredundant records for each attribute using simple linux commands (so might be inaccurate). This helped us for quick browsing possible values of each attributes.

Since the above initial observation suggested us that `ecotype` could be an attribute relate with ecotype information, we applied the following perl oneliner to extract (lower-cased) values of attribute `ecotype`.
```
wdlin@comp04:SOMEWHERE/ath$ cat ath_SRS_20240529.txt | perl -ne 'chomp; if(/<Attribute attribute_name="(.+?)".*?>(.+?)</){ print "$2\n" if $1 eq "ecotype" }' | perl -ne 'chomp; $hash{lc($_)}++; if(eof){ for $k (sort {$hash{$b}<=>$hash{$a}} keys %hash){ print "$k\t$hash{$k}\n" } }' > extraction/ecotype0.xls

wdlin@comp04:SOMEWHERE/ath$ head extraction/ecotype0.xls
col-0   5538
columbia        1906
col-0 (efo_0005148)     692
col0    250
landsberg erecta        223
col-0 (cs70000) 106
bay x sha ril   83
columbia (col-0)        81
columbia (efo_0005147)  75
wassilewskija   69
```
We intended to save the tab-delimited text file with extension `.xls` because it is convenient to do the next *manual* curation step by using Excel. By importing `ecotype0.xls` into Excel, we added one more column named `col0` for identifying those `ecotype` values that should be refering to arabidopsis col-0 ecotype. Here, Excel formulas like FIND can be used for some quick and inaccurate identification. No matter how, a manual confirmation is needed. In the following example, it was shown that ecotype col-0 could be recorded under attribute `ecotype` with values like `col-0`, `columbia`, `col-0 (efo_0005148)`, `col-0 (cs70000)`, ... and many others. Note that the *last* column `col0` contains values of `TRUE` and `FALSE`.

![Excel editing of ecotype0.xls](https://github.com/wdlingit/maccu/blob/main/pic/ecotype0_excel.png)

After the first round of manual curation, we saved the curation table for of col-0 *values* into a tab-delimited text file `ecotype0.txt` (quotes removed, if any). The following perl oneliner was applied to compute (i) attribute counts associated with curated col-0 *values* (in `ecotype0.txt`) and (ii) attribute counts in the metadata file. In so doing, we may discover attributes other than `ecotype` that also store ecotype information (recall that the metadata were human-inputted).
```
wdlin@comp04:SOMEWHERE/ath$ cat ath_SRS_20240529.txt | perl -ne 'if($.==1){ open(FILE,"<extraction/ecotype0.txt"); $line=<FILE>; while($line=<FILE>){ chomp $line; $line=~s/^\s+|\s+$//g; @s=split(/\t/,$line); $hash{$s[0]}=0 if $s[-1] eq "TRUE"; } close FILE } chomp; if(/<Attribute attribute_name="(.+?)".*?>(.+?)</){ $attr=$1; $val=lc($2); $cnt{$attr}++; $match{$attr}++ if exists $hash{$val} } if(eof STDIN){ print "attr\tmatch\ttotal\n"; for $attr (sort keys %match){ $x=0; $x=$match{$attr} if exists $match{$attr}; print "$attr\t$x\t$cnt{$attr}\n" } }' > extraction/ecotype1.xls

wdlin@comp04:SOMEWHERE/ath$ head extraction/ecotype1.xls
attr    match   total
Genotype        9       12
Matrial 22      26
Submitter Id    1       3115
accession       120     292
agent   29      135
background cultivar     9       9
background ecotype      114     215
background strain       18      18
cell line       2       63
```
Again, the output table was saved with extension `.xls` for importing to Excel for manual curation. In Excel, we added one `ratio` column that computes the ratio that an attribute associated with curated col-0 values. We also added a `selection` column that simply check if the ratio is greater than 0.1 or not. Note that the simple check on ratios would be convenient but not accurate. So, again, manual curation is needed. In the following picture, you may find that we execlude `Genotype` even if more than 0.1 of it appearances were associated with col-0 like values. It is also surprising (and actually not surprising) that attribute `accession` was found to contain ecotype information. The curated table of ecotype *attributes* was saved into a tab-delimited text file `ecotype1.txt`.

![Excel editing of ecotype0.xls](https://github.com/wdlingit/maccu/blob/main/pic/ecotype1_excel.png)

**FOR A SHORT SUMMARY**, now we have `ecotype1.txt` contains attributes we considered containing ecotype infromation in the metadata file. Note that the last column in `ecotype1.txt` is containing values of `TRUE` and `FALSE`. Also, in `ecotype0.txt`, we have values that we considered indicating col-0 ecotype, where `TRUE` and `FALSE` are under the col0 column.
```
wdlin@comp04:SOMEWHERE/ath$ head extraction/ecotype1.txt
attr    match   total   ratio   selection
Genotype        9       12      0.7500  FALSE
Matrial 22      26      0.8462  FALSE
Submitter Id    1       3115    0.0003  FALSE
accession       120     292     0.4110  TRUE
agent   29      135     0.2148  FALSE
background cultivar     9       9       1.0000  TRUE
background ecotype      114     215     0.5302  TRUE
background strain       18      18      1.0000  TRUE
cell line       2       63      0.0317  FALSE

wdlin@comp04:SOMEWHERE/ath$ head extraction/ecotype0.txt
value   count   col0
col-0   5538    TRUE
columbia        1906    TRUE
col-0 (efo_0005148)     692     TRUE
col0    250     TRUE
landsberg erecta        223     FALSE
col-0 (cs70000) 106     TRUE
bay x sha ril   83      FALSE
columbia (col-0)        81      TRUE
columbia (efo_0005147)  75      FALSE
```

So it is possible for us to iterate all samples in the metadata file and see if any possible ecotype attribute is assigned with a possible col-0 value for every sample. To do that, we applied the `biosampleClassify.pl` script. Note that it generates a *classification* matrix with the same number of columns as that in the `<valueFile>` file and the same number of rows as the number of SRS accessions in the metadata file. In the following example, it was shown that DRS014211 and DRS014212 are the first two SRS accessions considered not related with col-0. Note that our approach might not be fully accruate, but classified samples would be based on specified attributes and specified values in the metadata file.
```
wdlin@comp04:SOMEWHERE/ath$ ../scripts/biosampleClassify.pl
biosampleClassify.pl <attrFile> <valueFile> <biosampleXML>

wdlin@comp04:SOMEWHERE/ath$ ../scripts/biosampleClassify.pl extraction/ecotype1.txt extraction/ecotype0.txt ath_SRS_20240529.txt > extraction/ath_SRS_20240529.ecotype

wdlin@comp04:SOMEWHERE/ath$ head extraction/ath_SRS_20240529.ecotype
SRS     count   col0
DRS007600       0       1
ERS1174633      0       1
DRS007601       0       1
ERS1174634      0       1
DRS007602       0       1
ERS1174635      0       1
DRS014211       0       0
DRS014212       0       0
DRS073469       0       1
```

Again, we applied the `matrixSelection.pl` script to extract the potion of col-0 samples from a count matrix by taking the classification matrix as the selection matrix.
```
wdlin@comp04:SOMEWHERE/ath$ ../scripts/matrixSelection.pl extraction/ath_SRS_20240529.ecotype ath_sel20240529.nMatrix.txt.nondup extraction/ath_sel20240529.nMatrix col0
```

In our practice, we would do normalization on the count matrix of all collected col-0 samples here (refer above TMM method part for the normalization steps).

### Sample classification part 3, an example of arabidopsis tissues

Suppose that we have an attribute file `tissue1.txt` (like `ecotype1.txt` in above) and a value file (like `ecotype0` in above). We can similarly generate a classification matrix for tissues.
```
wdlin@comp04:SOMEWHERE/ath$ tail extraction/tissue1.txt
strain  8       204     0.0392  FALSE
tag     1       18      0.0556  FALSE
time    5       1453    0.0034  FALSE
tissue  12875   14761   0.8722  TRUE
tissue type     69      95      0.7263  TRUE
tissue/cell type        2       3       0.6667  TRUE
tissue_type     133     146     0.9110  TRUE
tissuie 3       3       1.0000  TRUE
tissus  2       23      0.0870  TRUE
treatment       16      7384    0.0022  FALSE

wdlin@comp04:SOMEWHERE/ath$ head extraction/tissue0.txt
value   count   leaf    rosette root    shoot   flower  seedling        seed    whole   CNT     OR
leaf    1659    TRUE    FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   1       TRUE
seedlings       1175    FALSE   FALSE   FALSE   FALSE   FALSE   TRUE    FALSE   FALSE   1       TRUE
root    1132    FALSE   FALSE   TRUE    FALSE   FALSE   FALSE   FALSE   FALSE   1       TRUE
leaves  1030    TRUE    FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   1       TRUE
seedling        1018    FALSE   FALSE   FALSE   FALSE   FALSE   TRUE    FALSE   FALSE   1       TRUE
shoot   600     FALSE   FALSE   FALSE   TRUE    FALSE   FALSE   FALSE   FALSE   1       TRUE
whole seedling  573     FALSE   FALSE   FALSE   FALSE   FALSE   TRUE    FALSE   FALSE   1       TRUE
whole seedlings 429     FALSE   FALSE   FALSE   FALSE   FALSE   TRUE    FALSE   FALSE   1       TRUE
whole plant     399     FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   TRUE    1       TRUE

wdlin@comp04:SOMEWHERE/ath$ ../scripts/biosampleClassify.pl extraction/tissue1.txt extraction/tissue0.txt ath_SRS_20240529.txt > extraction/ath_SRS_20240529.tissue

wdlin@comp04:SOMEWHERE/ath$ head extraction/ath_SRS_20240529.tissue
SRS     count   leaf    rosette root    shoot   flower  seedling        seed    whole   CNT     OR
DRS007600       0       0       0       0       0       1       0       0       0       0       1
ERS1174633      0       0       0       0       0       1       0       0       0       0       1
DRS007601       0       0       0       0       0       1       0       0       0       0       1
ERS1174634      0       0       0       0       0       1       0       0       0       0       1
DRS007602       0       0       0       0       0       1       0       0       0       0       1
ERS1174635      0       0       0       0       0       1       0       0       0       0       1
DRS014211       0       0       0       0       0       0       0       0       0       0       0
DRS014212       0       0       0       0       0       0       0       0       0       0       0
DRS073469       0       0       0       0       0       0       0       0       0       0       0
```

Given that we have the normalized log-count-per-million matrix of only col-0 samples saved in tab-delimited text file `sel20240529.Col0.TMM`, the following command can be applied for generating portions of tissues extracted from the normalized count matrix.
```
wdlin@comp04:SOMEWHERE/ath$ ../scripts/matrixSelection.pl ath_SRS_20240529.tissue coexDB_202406/ath/sel20240529.Col0.TMM coexDB_202406/ath/sel20240529.Col0.TMM leaf rosette root shoot flower seedling seed whole

wdlin@comp04:SOMEWHERE/ath$ find coexDB_202406/ath/ | perl -ne 'chomp; next if -d "$_"; print "$_\n"' | perl -ne 'chomp; $msg=`head -1 $_`; chomp $msg; @t=split(/\t/,$msg); $cnt=@t; $cnt--; print "$_\t$cnt\n"'
coexDB_202406/ath/sel20240529.TMM       19746
coexDB_202406/ath/sel20240529.Col0.TMM  9060
coexDB_202406/ath/sel20240529.Col0.TMM.flower   242
coexDB_202406/ath/sel20240529.Col0.TMM.leaf     2072
coexDB_202406/ath/sel20240529.Col0.TMM.root     1084
coexDB_202406/ath/sel20240529.Col0.TMM.rosette  587
coexDB_202406/ath/sel20240529.Col0.TMM.seed     323
coexDB_202406/ath/sel20240529.Col0.TMM.seedling 2854
coexDB_202406/ath/sel20240529.Col0.TMM.shoot    548
coexDB_202406/ath/sel20240529.Col0.TMM.whole    527
```
The last command is for numbers of data columns (samples) in the extracted matrixes.

### Sample classification part 4 (optional), iteratively refine attributes & values for selection

In above example, we started from one attrabute, collect *values* of our interests, decide *attributes*, and the use lastly adopted attributes and values for sample classification. Actually the process is flexible. For example, we can use those decided *attributes* to search more *values* for our decision. For example, the following perl oneliner was to use decided *attributes* (in file `dev1.txt`) to collect more *values* for making decision. Just remember to give an attribute file and a value file for generating a selection matrix.

```
wdlin@comp01:SOMEWHERE/dm$ cat dm_sel20240531.txt | perl -ne 'if($.==1){ open(FILE,"<extraction/dev1.txt"); while($line=<FILE>){ $line=~s/^\s+|\s+$//g; @s=split(/\t/,$line); $hash{$s[0]}=1 if $s[-1] eq "TRUE" } close FILE; } chomp; if(/<Attribute attribute_name="(.+?)".*?>(.+?)</){ print "$2\n" if exists $hash{$1} }' | perl -ne 'chomp; $hash{lc($_)}++; if(eof){ for $k (sort {$hash{$b}<=>$hash{$a}} keys %hash){ print "$k\t$hash{$k}\n" } }' > extraction/dev2.xls
```

### Sample classification part 5 (optional), import customized logic into the selection

Some coding approach could be applied because the above attribute-value selection method is simple for general scenario and may not fit some complex cases. For example, for fly, strain K-12 and substrain mg1655 can both be values for attributes *strain* and *substrain*. In this case, a sample might be considered as both K12 and mg1655. This could be true to some people but some other might want to separate those K-12 samples without any substrain info from samples with specific strain/substrain info of mg1655. In this case, we can fix the selection matrix by coding.

```
wdlin@comp01:SOMEWHERE/ec$ head -10  extraction/ec_sel20240531.strain
SRS     number  strain  mg1655  w3110   bw25113 ncm3722 K-12    b       rb001   ar3110  wo153   OR
DRS200394       0       1       1       0       0       0       0       0       0       0       0       1
DRS200395       0       1       1       0       0       0       0       0       0       0       0       1
DRS200396       0       1       1       0       0       0       0       0       0       0       0       1
DRS200397       0       1       1       0       0       0       0       0       0       0       0       1
DRS200398       0       1       1       0       0       0       0       0       0       0       0       1
ERS1139734      0       0       0       0       0       0       0       0       0       0       0       0
ERS1146196      0       0       0       0       0       0       0       0       0       0       0       0
ERS1203249      0       1       1       0       0       0       1       0       0       0       0       1
ERS1203251      0       1       1       0       0       0       1       0       0       0       0       1

wdlin@comp01:SOMEWHERE/ec$ head -10  extraction/ec_sel20240531.strain | perl -ne 'if($.==1){ print ; next } chomp; @t=split(/\t/); ($srs,$num,$strain,$mg1655,$w3110,$bw25113,$ncm3722,$k12,$b,$rb001,$ar3110,$wo153,$or)=@t; $k12=0 if ($mg1655 || $w3110 || $bw25113 || $ncm3722); @t=($srs,$num,$strain,$mg1655,$w3110,$bw25113,$ncm3722,$k12,$b,$rb001,$ar3110,$wo153,$or); print join("\t",@t)."\n"'
SRS     number  strain  mg1655  w3110   bw25113 ncm3722 K-12    b       rb001   ar3110  wo153   OR
DRS200394       0       1       1       0       0       0       0       0       0       0       0       1
DRS200395       0       1       1       0       0       0       0       0       0       0       0       1
DRS200396       0       1       1       0       0       0       0       0       0       0       0       1
DRS200397       0       1       1       0       0       0       0       0       0       0       0       1
DRS200398       0       1       1       0       0       0       0       0       0       0       0       1
ERS1139734      0       0       0       0       0       0       0       0       0       0       0       0
ERS1146196      0       0       0       0       0       0       0       0       0       0       0       0
ERS1203249      0       1       1       0       0       0       0       0       0       0       0       1
ERS1203251      0       1       1       0       0       0       0       0       0       0       0       1
```

In this case, we previously considered mg1655, w3110, bw25113, ncm3722 as K-12 substrains and would like to make samples marked with these substrain not being marked by K-12 (ex: ERS1203249 and ERS1203251). So the perl oneliner was to transfer a column of values into values (`($srs,$num,$strain,$mg1655,$w3110,$bw25113,$ncm3722,$k12,$b,$rb001,$ar3110,$wo153,$or)=@t`) and make some simple logic (`$k12=0 if ($mg1655 || $w3110 || $bw25113 || $ncm3722)`). In so doing, the outputted selection matrix should exclude those samples with specified substrain info from K-12.
