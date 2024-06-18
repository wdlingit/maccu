## Steps to build a co-expression database from DEE2 datasets

Steps in this document were done in Ubuntu 20 servers with 128GB memory. For human and mouse data, some steps may take up to more than 500GB memory.

### Processing the metadata file and aggregate the read counts

The metadata tables made by DEE2 was used for the initial sample qualification. The following steps were done using Excel.

1. The metadata tables provide QC results for SRR's, which rather correspond to technical replicates. SRR's are the basic records for samples in DEE2. To qualify biological replicates, i.e., SRS accessions, which also in the metadata tables, we collected SRS accessions where their corresponding SRR's were all PASS in the QC column.
2. To collect count data associated with SRR's from the DEE2 database, we collected SRR accessions (i) under SRS accessions collected in step1, (ii) with `experiment_library_strategy` of `RNA-Seq`, and (iii) `experiment_library_selection` with `cDNA`, `RANDOM`, `PolyA`, or `Oligo-dT`.
3. Save SRR-SRS mapping (two columns) into a tab-delimited text file.

Download the count file and use the perl oneliner command like the following example to extract and aggregate read counts into biological replicates, i.e., SRS accessions.
```
wdlin@login02:SOMEWHERE/ath$ head ath_SRR_20240529.txt
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
2. `ath_SRR_20240529.txt` is the SRR-SRS mapping (a two-column tab-delimited text file) with collected SRS accessions in above step1.
3. The output file `ath_sel20240529.nMatrix.txt` (tab-delimited) is the raw count matrix, with columns for samples and rows for genes.

### Duplicate removal

Some samples (SRS) would be repeatedly submitted to the NCBI SRA database. The following steps were applied for removing duplications.

```
wdlin@comp01:SOMEWHERE/ath$ ../scripts/duplicateDetect.pl ath_sel20240529.nMatrix.txt > ath_sel20240529.nMatrix.dupReport

wdlin@login02:SOMEWHERE/ath$ head ath_sel20240529.nMatrix.dupReport
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

wdlin@comp01:SOMEWHERE/ath$ head -1 ath_sel20240529.nMatrix.txt | perl -ne 'chomp; s/^\s+|\s+$//g; @t=split; print "$_\n" for @t' | perl -ne 'chomp; if($.==1){ open(FILE,"<ath_sel20240529.nMatrix.dupReport"); while($line=<FILE>){ chomp $line; if($line=~/^DUP/){ @s=split(/\s+/,$line); shift @s; shift @s; for $x (@s){ $duplicate{$x}=1 } }} close FILE; print "$_\tnondup\n" }else{ print "$_\t"; if(exists $duplicate{$_}){ print "0\n" }else{ print "1\n" } }' > ath_sel20240529.dup.txt

wdlin@login02:SOMEWHERE/ath$ head ath_sel20240529.dup.txt
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

The first command was to use the script `duplicateDetect.pl` (in our `scripts` directory) for identifying duplicate columns (samples). Inside the output file (`ath_sel20240529.nMatrix.dupReport` here), lines started with `DUP:` are for duplicated samples. The perl oneliner was to read the header column from the raw count matrix (`head -1 ath_sel20240529.nMatrix.txt`) and generate a 0-1 matrix (`ath_sel20240529.dup.txt`) based on the duplication report. The 0-1 matrix was for indicating which samples are nonduplicated. For samples reported in the ducplication report, only the first samples were specified as nonduplicated.

The last command for duplication removal was to apply `matrixSelection.pl` (in our `scripts` directory). This script takes at least four parameters:
1. selection matrix: in this case, `ath_sel20240529.dup.txt` is the selection matrix. Note that column headers are treated as selection targets.
2. source matrix: a tab-delimited matrix file, with columns for samples.
3. output prefix: a output filename would be in the form `<output prefix>.<target>`.
4. targets: one or more column headers from the selection matrix could be selected for column selection from the source matrix.
```
wdlin@login02:SOMEWHERE/ath$ ../scripts/matrixSelection.pl
matrixSelection.pl <selMatrix> <sourceMatrix> <outPrefix> [<selTarget>]+

wdlin@comp01:SOMEWHERE/ath$ ../scripts/matrixSelection.pl ath_sel20240529.dup.txt ath_sel20240529.nMatrix.txt ath_sel20240529.nMatrix.txt nondup
```

In this example, the read count matrix without duplicated samples would be named `ath_sel20240529.nMatrix.txt.nondup`.

### In case no sample classification required

We believe that certain normalization is needed for the co-expression database but not the row counts. In case that no sample classification requried. The following R commands were adopted for the normalization task using the TMM method (PMID: 20196867).
```
library("edgeR")
x <- read.delim("ath_sel20240529.nMatrix.txt.nondup",row.names="Symbol")
dge <- DGEList(counts=x)
dge <- calcNormFactors(dge)
v <- voom(dge,normalize="none")
write.csv(x=v$E,file="sel20240529.nMatrix.TMM")
quit()
```

You may modify the `write.csv` command or use the following perl oneliner to transfer the output CSV file into a tab-delimited text file.

```
wdlin@comp04:SOMEWHERE/ath$ cat sel20240529.nMatrix.TMM | perl -ne 'chomp; @t=split(/,/); $nonFirst=0; for $x (@t){ $x=~s/^"|"$//g; print "\t" if $nonFirst; $nonFirst=1; print "$x" } print "\n"' > sel20240529.nMatrix.TMM.txt
```

### Simple classification part 1, downloading metadata from SRA


