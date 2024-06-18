## Steps to build co-ex db from DEE2 datasets

### Processing the metadata file and aggregate the read counts

The metadata tables made by DEE2 was used for the initial sample qualification. The following steps were done using Excel.

1. The metadata tables provide QC results for SRR's, which rather correspond to technical replicates. SRR's are the basic records for samples in DEE2. To qualify biological replicates, i.e., SRS accessions, which also in the metadata tables, we collected SRS accessions where their corresponding SRR's were all PASS in the QC column.
2. To collect count data associated with SRR's from the DEE2 database, we collected SRR accessions (i) under SRS accessions collected in step1, (ii) with `experiment_library_strategy` of `RNA-Seq`, and (iii) `experiment_library_selection` with `cDNA`, `RANDOM`, `PolyA`, or `Oligo-dT`.
3. Save SRR-SRS mapping (two columns) into a tab-delimited text file.

Download the count file and use the perl-oneliner command like the following example to extract and aggregate read counts into biological replicates, i.e., SRS accessions.
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
3. The output file `ath_sel20240529.nMatrix.txt` is the raw count matrix, with columns for samples and rows for genes.
