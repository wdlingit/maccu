## Steps to build a co-expression database from DEE2 datasets

This document contains steps for downloading specified SRS metadata and classification of samples. Steps in this document were done in a Ubuntu 20 server with 128GB memory. For human and mouse data, some steps may take up to more than 500GB memory.

**NOTE**: The difference between this version and [the last version](https://github.com/wdlingit/maccu/blob/main/DB_from_DEE2.000000.md) is that ALL manual Excel operations were replaced by perl codes. TRUE and FALSE values used in Excel were replaced by 1's and 0's.

### Processing the metadata file and aggregate the read counts

The [metadata tables made by DEE2](https://dee2.io/metadata/) were used for the initial sample qualification. In this document, we use `athaliana_metadata.tsv` for describing the methods.

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ wget https://dee2.io/metadata/athaliana_metadata.tsv

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head -1 athaliana_metadata.tsv | perl -ne 'chomp; @t=split; for($i=0;$i<@t;$i++){ print "$i\t$t[$i]\n" }'
0       SRR_accession
1       QC_summary
2       SRX_accession
3       SRS_accession
4       SRP_accession
5       GEO_series
6       Experiment_title
```
The metadata tables provide QC results of SRR accessions, which rather correspond to technical replicates. SRR's are the basic records in DEE2. To qualify biological replicates, i.e., SRS accessions, which also in the metadata tables, we collected SRS accessions where their corresponding SRR's were all PASS in the QC column.

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat athaliana_metadata.tsv | perl -ne 'next if $.==1; chomp; @t=split; $srsHash{$t[3]}{$t[0]}=1; if($t[1]=~/pass/i){ $srrPassHash{$t[0]}=1 }else{ $srrPassHash{$t[0]}=0 } if(eof){ ($cnt,$pCnt,$fCnt)=(0,0,0); for $k (keys %srrPassHash){ $cnt++; if($srrPassHash{$k}){ $pCnt++; }else{ $fCnt++; } } print "SRR: $cnt\t$pCnt\t$fCnt\n"; ($cnt,$pCnt,$fCnt)=(0,0,0); for $k1 (keys %srsHash){ $cnt++; $flag=1; for $k2 (keys %{$srsHash{$k1}}){ $flag=0 if $srrPassHash{$k2}==0; } if($flag){ $pCnt++; }else{ $fCnt++; } } print "SRS: $cnt\t$pCnt\t$fCnt\n"; }'
SRR: 119477     50063   69414
SRS: 90960      44495   46465

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat athaliana_metadata.tsv | perl -ne 'next if $.==1; chomp; @t=split; $srsHash{$t[3]}{$t[0]}=1; if($t[1]=~/pass/i){ $srrPassHash{$t[0]}=1 }else{ $srrPassHash{$t[0]}=0 } if(eof){ for $k1 (sort keys %srsHash){ $flag=1; for $k2 (keys %{$srsHash{$k1}}){ $flag=0 if $srrPassHash{$k2}==0; } if($flag){ for $k2 (sort keys %{$srsHash{$k1}}){ print "$k1\t$k2\n" } } } }' > SRS_SRR.allpass
```
In above, the first perl oneliner showed total, pass, and non-pass numbers of SRR's, as well as, total, all-pass, non-all-pass numbers of SRS's. Here we have 44495 SRS's with all-pass SRR's. The second perl oneliner saved SRS-SRR pairs of all-pass SRS's.

To collect SRS's that are belonging to RNAseq samples, we need metadata info in addition to those provided in `athaliana_metadata.tsv`. Script `runinfoRetrieve.pl` (in our `scripts` directory) was used for retriving metadata of SRR's from NCBI.

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/runinfoRetrieve.pl
runinfoRetrieve.pl <srrList> <processed> <outCSV>
```
Points to be noticed:
1. The [NCBI EDirect utility](https://www.ncbi.nlm.nih.gov/books/NBK179288/) is required for running this script
2. This script will write retrieved metadata in CSV format into `<outCSV>` and processed SRR's into `<processed>`. You may use the line numbers in `<processed>` to check numbers of SRR records with successfully retrieved metadata.
3. This script will *append* contents to the two output files, and it will process only SRR accessions not in `<processed>`. That is, you may simply repeat the same command a few times for retrieving metadata for the same SRR list without taking care of the outputs. NOTE: It is possible that the NCBI contains no metadata for some SRR accessions. Just remove those SRR accessions kept being searched for a number of times and check them in the NCBI webpage.
4. This script doesn't support parallel processing. You may apply a command like `split -l 6025 SRS_SRR.allpass.SRR SRS_SRR.allpass.SRR.` to split the list into smaller lists for parallel processing (surely separate output files for separate input lists). Note that NCBI has some query number restriction per second given an API key. Be sure not to exceed the limitation.
5. Variable `$maxTry` was hard-coded as `3` for the number of re-try an `efetch` command.
6. Variable `$chunkSize` was hard-coded as `100` so that every 100 SRR accessions would be queried by one single command. This would largely improve the query efficiency. In our experiences, 6000 SRRs would took only a few minutes.

Assuming that `SRS_SRR.allpass.SRR` is for the SRR list and `SRS_SRR.allpass.SRR.out` is the output CSV file of `runinfoRetrieve.pl`. The next perl oneliner extracts `LibraryStrategy`, `LibrarySource`, `LibrarySelection`, `Sample`, and `BioSample` from the CSV output. Note that it requires the `Text::CSV` perl module.

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat SRS_SRR.allpass.SRR.out | perl -MText::CSV -ne 'if($.==1){ open(FILE,"<SRS_SRR.allpass.SRR"); while($line=<FILE>){ chomp $line; $hash{$line}=1; } close FILE; @attrArr=("LibraryStrategy","LibrarySource","LibrarySelection","Sample","BioSample"); $csv=Text::CSV->new({ binary => 1, auto_diag => 1 }); } chomp; if($csv->parse($_)){ @t=$csv->fields }else{ die "ERROR: $_\n" } if($t[0] eq "Run"){ %idxHash=(); for($i=0;$i<@t;$i++){ $idxHash{$t[$i]}=$i } }elsif(exists $hash{$t[0]}){ print "$t[0]"; for $k (@attrArr){ if(exists $idxHash{$k}){ print ",$t[$idxHash{$k}]" }else{ print "," } } print "\n" }' | sort | uniq > SRS_SRR.allpass.SRR.info

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head SRS_SRR.allpass.SRR.info
DRR008476,RNA-Seq,TRANSCRIPTOMIC,cDNA,DRS007600,SAMD00009103
DRR008477,RNA-Seq,TRANSCRIPTOMIC,cDNA,DRS007601,SAMD00009101
DRR008478,RNA-Seq,TRANSCRIPTOMIC,cDNA,DRS007602,SAMD00009102
DRR016112,RNA-Seq,TRANSCRIPTOMIC,cDNA,DRS014211,SAMD00013248
DRR016113,RNA-Seq,TRANSCRIPTOMIC,cDNA,DRS014211,SAMD00013248
DRR016114,RNA-Seq,TRANSCRIPTOMIC,cDNA,DRS014212,SAMD00013247
DRR016115,RNA-Seq,TRANSCRIPTOMIC,cDNA,DRS014212,SAMD00013247
DRR016116,RNA-Seq,TRANSCRIPTOMIC,cDNA,DRS014212,SAMD00013247
DRR018424,RNA-Seq,TRANSCRIPTOMIC,unspecified,DRS016105,SAMD00015876
DRR021335,RNA-Seq,TRANSCRIPTOMIC,RANDOM,DRS030798,SAMD00018417
```

There are some SRS accessions been submitted with multiple *library strategies* and/or *library selections*. They were excluded for safty. The last three columns are for numbers of different *strategies*, *sources*, and *selections*, respectively.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat bck/SRS_SRR.allpass.SRR.info | perl -ne 'if($.==1){ open(FILE,"<SRS_SRR.allpass"); while($line=<FILE>){ chomp $line; @s=split(/\s+/,$line); $hash{$s[1]}=$s[0]; } close FILE } chomp; @t=split(/,/,$_,-1); unshift @t,$hash{$t[0]}; print join("\t",@t)."\n";' | perl -ne 'chomp; @t=split; $lineHash{$t[1]}=$_; $srrHash{$t[1]}=$t[0]; $hash2{$t[0]}{$t[2]}=1; $hash3{$t[0]}{$t[3]}=1; $hash4{$t[0]}{$t[4]}=1; if(eof){ for $k (sort keys %lineHash){ $cnt2=keys %{$hash2{$srrHash{$k}}}; $cnt3=keys %{$hash3{$srrHash{$k}}}; $cnt4=keys %{$hash4{$srrHash{$k}}}; print "$lineHash{$k}\t$cnt2\t$cnt3\t$cnt4\n"; } }' | perl -ne '@t=split; print if $t[-3]>1 || $t[-2]>1 || $t[-1]>1' | sort | head
DRS235147       DRR221886       OTHER   TRANSCRIPTOMIC  cDNA    DRS235147       SAMD00218963    2       1       1
DRS235147       DRR221902       RNA-Seq TRANSCRIPTOMIC  cDNA    DRS235147       SAMD00218963    2       1       1
DRS235148       DRR221887       OTHER   TRANSCRIPTOMIC  cDNA    DRS235148       SAMD00218964    2       1       1
DRS235148       DRR221903       RNA-Seq TRANSCRIPTOMIC  cDNA    DRS235148       SAMD00218964    2       1       1
DRS235149       DRR221888       OTHER   TRANSCRIPTOMIC  cDNA    DRS235149       SAMD00218965    2       1       1
DRS235149       DRR221904       RNA-Seq TRANSCRIPTOMIC  cDNA    DRS235149       SAMD00218965    2       1       1
DRS235150       DRR221889       OTHER   TRANSCRIPTOMIC  cDNA    DRS235150       SAMD00218966    2       1       1
DRS235150       DRR221905       RNA-Seq TRANSCRIPTOMIC  cDNA    DRS235150       SAMD00218966    2       1       1
DRS235152       DRR221891       OTHER   TRANSCRIPTOMIC  cDNA    DRS235152       SAMD00218968    2       1       1
DRS235152       DRR221907       RNA-Seq TRANSCRIPTOMIC  cDNA    DRS235152       SAMD00218968    2       1       1
```

Finally we collected SRR accessions (i) with `library strategy` of `RNA-Seq`, and (ii) `library selection` be `cDNA`, `RANDOM`, `PolyA`, or `Oligo-dT`.

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat SRS_SRR.allpass.SRR.info | perl -ne 'if($.==1){ open(FILE,"<SRS_SRR.allpass"); while($line=<FILE>){ chomp $line; @s=split(/\s+/,$line); $hash{$s[1]}=$s[0]; } close FILE } chomp; @t=split(/,/,$_,-1); unshift @t,$hash{$t[0]}; print join("\t",@t)."\n";' | perl -ne 'chomp; @t=split; $lineHash{$t[1]}=$_; $srrHash{$t[1]}=$t[0]; $hash2{$t[0]}{$t[2]}=1; $hash3{$t[0]}{$t[3]}=1; $hash4{$t[0]}{$t[4]}=1; if(eof){ for $k (sort keys %lineHash){ $cnt2=keys %{$hash2{$srrHash{$k}}}; $cnt3=keys %{$hash3{$srrHash{$k}}}; $cnt4=keys %{$hash4{$srrHash{$k}}}; print "$lineHash{$k}\n" if $cnt2==1 && $cnt3==1 && $cnt4==1; } }' | perl -ne 'chomp; @t=split; print "$t[1]\t$t[0]\n" if ($t[2] eq "RNA-Seq") && (($t[4] eq "cDNA") || ($t[4] eq "Oligo-dT") || ($t[4] eq "RANDOM") || ($t[4] eq "PolyA"))' > SRS_SRR.selected

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ wc -l SRS_SRR.selected
40363 SRS_SRR.selected
```

Download [the count file](https://dee2.io/mx/) and use the `SRS_aggr.R` script (in our `scripts` directory, requires the `rhdf5` library) to aggregate read counts into biological replicates, i.e., SRS accessions.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head SRS_SRR.selected
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

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/SRS_aggr.R athaliana_se.h5 SRS_SRR.selected sel20261001.nMatrix.txt
```
Points to be noticed:
1. `athaliana_se.h5` is the count file downloaded from the DEE2 database. It is in the HDF5 format.
2. `SRS_SRR.selected` is the SRR-SRS mapping (a two-column tab-delimited text file) with SRS accessions collected in the above step.
3. The output file `sel20261001.nMatrix.txt` (tab-delimited) is the raw count matrix, with columns for samples and rows for genes.

### Duplicate removal

Some samples (SRS) might be repeatedly submitted to the NCBI SRA database. The following steps were applied for removing duplications.

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/duplicateDetect.pl sel20261001.nMatrix.txt > sel20261001.nMatrix.dupReport

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head sel20261001.nMatrix.dupReport
Reading matrix
Compute hash
Compare
Report
DUP: DRS518865  SRS26357459
DUP: DRS518866  SRS26357460
DUP: DRS518867  SRS26357461
DUP: ERS14405101        SRS9822351
DUP: ERS14405102        SRS9822352
DUP: ERS14405103        SRS9822347

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head -1 sel20261001.nMatrix.txt | perl -ne 'chomp; s/^\s+|\s+$//g; @t=split; print "$_\n" for @t' | perl -ne 'chomp; if($.==1){ open(FILE,"<sel20261001.nMatrix.dupReport"); while($line=<FILE>){ chomp $line; if($line=~/^DUP/){ @s=split(/\s+/,$line); shift @s; shift @s; for $x (@s){ $duplicate{$x}=1 } }} close FILE; print "$_\tnondup\n" }else{ print "$_\t"; if(exists $duplicate{$_}){ print "0\n" }else{ print "1\n" } }' > sel20261001.dup.txt

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head sel20261001.dup.txt
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
The first command was to use the script `duplicateDetect.pl` (in our `scripts` directory) to identify duplicate columns (samples). Inside the output file (`sel20261001.nMatrix.dupReport` here), lines started with `DUP:` are for duplicated samples. The perl oneliner was to read the header column from the raw count matrix (`head -1 sel20261001.nMatrix.txt`) and generate a 0-1 matrix (`sel20261001.dup.txt`) based on the duplication report. The 0-1 matrix was for indicating which samples are nonduplicated. For samples reported in the ducplication report, only the first sample from each line was specified as nonduplicated.

The last command for duplication removal was to apply `matrixSelection.pl` (in our `scripts` directory). This script takes at least four parameters:
1. selection matrix: in this case, `sel20261001.dup.txt` is the selection matrix. Note that column headers are treated as selection targets.
2. source matrix: a tab-delimited matrix file, with columns for samples.
3. output prefix: an output filename would be in the form `<output prefix>.<target>`.
4. targets: one or more column headers from the selection matrix could be selected for column selection from the source matrix.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/matrixSelection.pl
matrixSelection.pl <selMatrix> <sourceMatrix> <outPrefix> [<selTarget>]+

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/matrixSelection.pl sel20261001.dup.txt sel20261001.nMatrix.txt sel20261001.nMatrix.txt nondup
```

In this example, the read count matrix without duplicated samples would be named `sel20261001.nMatrix.txt.nondup`.

### In case no sample classification required

We believe that certain normalization is needed for the co-expression database but not the raw counts. In case that no sample classification requried. The R script `TMM.R` (in our `scripts` directory, requires the `edgeR` library) was adopted for the normalization task using the TMM method (PMID: 20196867).

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/TMM.R sel20261001.nMatrix.txt.nondup Symbol sel20261001.nMatrix.TMM
```

Its output is a tab-delimited matrix file.

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat sel20261001.nMatrix.TMM | perl -ne 'chomp; @t=split(/\t/); $cnt=@t; $hash{$cnt}++; if(eof){ for $x (sort keys %hash){ print "$x\t$hash{$x}\n" } }'
37809   32834
```

### Sample classification part 1, downloading metadata from NCBI

We firstly generate a list of SRS accessions of nonduplicated samples. You may take `sel20261001.nMatrix.txt.nondup` as the input if you didn't do the TMM step.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head -1 sel20261001.nMatrix.TMM | perl -ne 'chomp; @t=split; for $x (@t){ print "$x\n" if length($x)>0 }' > sel20261001.SRRs
```

Then generate a hint file of SRS to BioSample mapping by extracting them from `SRS_SRR.allpass.SRR.info`.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat SRS_SRR.allpass.SRR.info | perl -ne 'chomp; @t=split(/,/,$_,-1); $cnt=@t; $hash{$t[4]}=$t[5]; if(eof){ for $k (sort keys %hash){ print "$k\t$hash{$k}\n" } }' > SRS_bios.hints

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head SRS_bios.hints
DRS007600       SAMD00009103
DRS007601       SAMD00009101
DRS007602       SAMD00009102
DRS014211       SAMD00013248
DRS014212       SAMD00013247
DRS016105       SAMD00015876
DRS030797       SAMD00018418
DRS030798       SAMD00018417
DRS047331       SAMD00060395
DRS047332       SAMD00060396
```

The `biosampleRetrieveBySRS.pl` (in our `scripts` directory) was used for retriving metadata of BioSamples from NCBI. The reason of using BioSample metadata is because it is richer than that of SRA records so we better obtain BioSample accessions for SRS records. The script contains two edirect queries for each SRS, the first one is for obtaining corresponding BioSample accessions, and the second one is for getting BioSample metadata. With the help of the hint file, we can save the efforts of getting BioSample accessions.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/biosampleRetrieveBySRS.pl
biosampleRetrieve.pl <listFile> <outSrsBios> <outXML>

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/biosampleRetrieveBySRS.pl -hint SRS_bios.hints sel20261001.SRRs sel20261001.SRRs.bios sel20261001.xml
```
Points to be noticed:
1. The [NCBI EDirect utility](https://www.ncbi.nlm.nih.gov/books/NBK179288/) is required for running this script
2. This script will write retrieved metadata XML into `<outXML>` and processed SRS-BioSample accession pairs into `<outSrsBios>`. You may use the line numbers in `<outSrsBios>` to check numbers of SRS records with successfully retrieved metadata.
3. This script will *append* contents to the two output files, and it will process only SRS accessions not in `<outSrsBios>`. That is, you may simply repeat the same command a few number of times for retrieving metadata for the same list without taking care of the outputs. NOTE: It is possible that the NCBI contains no metadata for some SRS accessions. Just mark those SRS accessions kept being searched for a number of times and check them in the NCBI webpage.
4. This script doesn't support parallel processing. You may apply a command like `split -l 4726 -d sel20261001.SRRs sel20261001.SRRs.` to split the list into smaller lists for parallel processing (surely separate output files for separate input lists). Note that NCBI has some query number restriction per second given an API key. Be sure not to exceed the limitation.
5. Option `-hint` is for providing the hint of BioSample accessions for SRS accessions.
6. Option `-maxtry` (default `3`) is for the number of re-try an `esearch` command. Modify it if needed.
6. Inside the script, the first two `esearch` commands were used for retrieving the corresponding BioSample accession of an SRS accession. They are our current best practices for retrieving BioSample accessions from SRS accessions. Modify them if needed.
7. The last `esearch` command in the script was to extract metadata of the BioSample accession corresponding to an SRS accession.

### Sample classification part 2, an example of arabidopsis ecotypes

Due to the complexity of human-input metadata, we *currently* don't have a completely automatic classification method. Here we present an approximation that used to give enough number of samples after classification. In this session, we present what we had done on arabidopsis ecotypes.

Suppose that `sel20261001.xml` is the metadata XML file that we obtained using the script described in the last session. The following two commands helped us for understanding the diversity inside the metadata.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat sel20261001.xml | perl -ne 'chomp; if(/<Attribute attribute_name="(.+?)"/){ $hash{"$1"}++ } if(eof){ for $k (sort {$hash{$b}<=>$hash{$a}} keys %hash){ print "$k\t$hash{$k}\n" } }' > sel20261001.attributes

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head sel20261001.attributes
tissue  29769
source_name     21917
genotype        21474
ecotype 19706
geo_loc_name    18336
age     16510
treatment       15898
collection_date 12488
dev_stage       7251
isolate 5030
```
In the metadata XML file, each BioSample is associated with a number of *attributes* which may have different values. For example, samples may have `tissue` attributes of values `root`, `seedling`, .... The first perl oneliner was to collect all attributes and rank them from the most frequently recorded attribute to the least frequently recorded attribute. In file `sel20261001.attributes`, we found that `tissue` was ranked first, which *should be* corresponding to tissue information. It was also found that the `ecotype` attribute was ranked fourth and that *should be* corresponding to ecotype information.

Since the above initial observation suggested us that `ecotype` could be an attribute relate with ecotype information, we applied the following perl oneliner to extract (lower-cased) values of attribute `ecotype`.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat sel20261001.xml | perl -ne 'chomp; if(/<Attribute attribute_name="(.+?)".+?display_name="(.+?)".*?>(.+?)</){ print "$3\n" if $1 eq "ecotype" }' | perl -ne 'chomp; $hash{lc($_)}++; if(eof){ print "value\tcount\n"; for $k (sort {$hash{$b}<=>$hash{$a}} keys %hash){ print "$k\t$hash{$k}\n" } }' > extraction/ecotype0.in

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head extraction/ecotype0.in
value   count
col-0   10792
columbia        3924
col0    752
col-0 (efo_0005148)     606
landsberg erecta        364
missing 254
columbia (col-0)        152
not collected   104
col-0 (cs70000) 90
```

The next perl oneliner command was for simple (and dirty) ecotype classification for col-0 and ler. It reads `ecotype0.in` and appends 0-1 columns for each classified ecotypes. Some rules were found incorrect and fixed. (ex: `efo_0005154`)
```
cat extraction/ecotype0.in |
perl -ne '
    chomp;
    if($.==1){ print "$_\tcol0\tler\tOR\n"; next; }
    ($value,$count)=split(/\t/);
    if($value=~/col-0/ ||
       $value=~/col_0/ ||
       $value=~/col0/ ||
       $value=~/col - 0/ ||
       $value=~/columbia background/ ||
       $value=~/^columbia$/ ||
       $value=~/^coloumbia$/ ||
       $value=~/^colombia 0$/ ||
       $value=~/^colombia-0$/ || 
       $value=~/^columbia - 0$/ || 
       $value=~/columia-0/ || 
       $value=~/^columbia0$/ || 
       $value=~/^columbia_0$/ || 
       $value=~/efo_0005147/ || 
       $value=~/col 0 \(cs70000\)/ || 
       $value=~/a.thalianaecotype columbia/ || 
       $value=~/columbia-0 ecotype/ || 
       $value=~/columbia-0 background/ || 
       $value=~/wild-type columbia-0/ || 
       $value=~/columbia-0 \(/ || 
       $value=~/columbia \(/){
        $col0=1 
    }else{ 
        $col0=0 
    } 
    if($value=~/landsberg erecta/ || 
       $value=~/lansberg erecta/ || 
       $value=~/landsberg ecotype/ || 
       $value=~/efo_0005154/ || 
       $value=~/ler-0/){ 
        $ler=1 
    }else{ 
        $ler=0 
    } 
    if($col0 || $ler){ 
        $or=1 
    }else{ 
        $or=0 
    } 
    print "$_\t$col0\t$ler\t$or\n";
' > extraction/ecotype0.out
```

This needs some iterative fix of the rules. The next two perl oneliners can provide some help.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat extraction/ecotype0.out | perl -ne 'chomp; @t=split(/\t/); print "$_\n" if $t[-1]==0' | head
value   count   col0    ler     OR
missing 254     0       0       0
not collected   104     0       0       0
bay x sha ril   83      0       0       0
not applicable  82      0       0       0
wassilewskija   81      0       0       0
tre-1   48      0       0       0
mutant  47      0       0       0
tol-0   47      0       0       0
abd-0   47      0       0       0

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat extraction/ecotype0.out | perl -ne 'chomp; if($.==1){ @header=split(/\t/); next } @t=split(/\t/); $hash{"count"}+=$t[1]; for($i=2;$i<@t;$i++){ $hash{$header[$i]}+=$t[1] if $t[$i] } if(eof){ shift @header; for $x (@header){ print "$x\t$hash{$x}\n" } }'
count   19706
col0    16983
ler     428
OR      17411
```
The first one lists `values` that were not taken into consideration. The second one summarizes numbers of classifications.

As the curation table of col-0 and ler was saved in `ecotype0.out`, the following perl oneliner was applied to compute (i) attribute counts associated with curated col-0 or ler *values* (in `ecotype0.out`) and (ii) attribute counts in the metadata file. In so doing, we may discover attributes other than `ecotype` that also store ecotype information (recall that the metadata were human-inputted).
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat sel20261001.xml | perl -ne 'if($.==1){ open(FILE,"<extraction/ecotype0.out"); $line=<FILE>; while($line=<FILE>){ chomp $line; $line=~s/^\s+|\s+$//g; @s=split(/\t/,$line); $hash{$s[0]}=0 if $s[-1]==1; } close FILE } chomp; if(/<Attribute attribute_name="(.+?)".*?>(.+?)</){ $attr=$1; $val=lc($2); $cnt{$attr}++; $match{$attr}++ if exists $hash{$val} } if(eof STDIN){ print "attr\tmatch\ttotal\tratio\n"; for $attr (sort {$match{$b}<=>$match{$a}} keys %match){ $x=0; $x=$match{$attr} if exists $match{$attr}; print "$attr\t$x\t$cnt{$attr}\t".sprintf("%.2f",($x/$cnt{$attr}))."\n" } }' > extraction/ecotype1.in

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head extraction/ecotype1.in
attr    match   total   ratio
ecotype 17411   19706   0.88
genotype        2218    21474   0.10
ecotype background      406     513     0.79
accession       256     452     0.57
ecotype/background      225     232     0.97
cultivar        205     3685    0.06
background ecotype      166     279     0.59
plant line      114     172     0.66
genetic background      90      377     0.24
```

Again, we use the following perl oneliner for writing rules of curation.
```
cat extraction/ecotype1.in | 
perl -ne '
    chomp; 
    @t=split(/\t/); 
    if($.==1){ 
        print "$_\tselection\n"; 
        next 
    } 
    $sel=0; 
    $sel=1 if $t[-1]>0.5; 
    if(($t[0] eq "cell line") || 
       ($t[0] eq "cell type") || 
       ($t[0] eq "cell_line") || 
       ($t[0] eq "cultivar") || 
       ($t[0] eq "genetic background") || 
       ($t[0] eq "genetic background ecotype") || 
       ($t[0] eq "subspecific genetic lineage name") || 
       ($t[0] eq "strain/background") || 
       ($t[0] eq "strain/ecotype") || 
       ($t[0] eq "strain/line")){ 
        $sel=1 
    } 
    if(($t[0] eq "female parent strain")){ 
        $sel=0 
    }  
    push @t,$sel; 
    print join("\t",@t)."\n"
' > extraction/ecotype1.out
```

**FOR A SHORT SUMMARY**, now we have `ecotype1.out` contains attributes we considered containing ecotype infromation in the metadata file. Note that the last column in `ecotype1.out` is containing values of 1's and 0's. Also, in `ecotype0.out`, we have values that we considered indicating col-0 or ler ecotype.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head extraction/ecotype1.out
attr    match   total   selection
ecotype 17411   19706   0.88    1
genotype        2218    21474   0.10    0
ecotype background      406     513     0.79    1
accession       256     452     0.57    1
ecotype/background      225     232     0.97    1
cultivar        205     3685    0.06    1
background ecotype      166     279     0.59    1
plant line      114     172     0.66    1
genetic background      90      377     0.24    1

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head extraction/ecotype0.out
value   count   col0    ler     OR
col-0   10792   1       0       1
columbia        3924    1       0       1
col0    752     1       0       1
col-0 (efo_0005148)     606     1       0       1
landsberg erecta        364     0       1       1
missing 254     0       0       0
columbia (col-0)        152     1       0       1
not collected   104     0       0       0
col-0 (cs70000) 90      1       0       1
```

So it is possible for us to iterate all samples in the metadata file and see if any possible ecotype attribute is assigned with a possible col-0 or ler value for every sample. To do that, we applied the `biosampleClassify.pl` script. Note that it generates a *classification* matrix with the same number of columns as that in the `<valueFile>` file and the same number of rows as the number of SRS accessions in the metadata file. In the following example, it was shown that DRS014211 and DRS014212 are the first two SRS accessions considered not related with col-0. Note that our approach might not be fully accurate, but classified samples would be based on specified attributes and specified values in the metadata file.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/biosampleClassify.pl
biosampleClassify.pl <attrFile> <valueFile> <biosampleXML>

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/biosampleClassify.pl extraction/ecotype1.out extraction/ecotype0.out sel20261001.xml > test.out

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head test.out
SRS     count   col0    ler     OR
DRS007600       0       1       0       1
ERS1174633      0       1       0       1
DRS007601       0       1       0       1
ERS1174634      0       1       0       1
DRS007602       0       1       0       1
ERS1174635      0       1       0       1
DRS014211       0       0       0       0
DRS014212       0       0       0       0
DRS030797       0       1       0       1
```

Here, it is possible to apply the `matrixSelection.pl` script to extract the potion of col-0 or ler samples from a count matrix by taking the classification matrix as the selection matrix.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/matrixSelection.pl test.out sel20261001.nMatrix.txt.nondup testOut col0 ler
```

In our practice, we would do normalization on the count matrix of all collected col-0 and ler samples here (refer above TMM method part for the normalization step).

### Sample classification part 3, an example of arabidopsis tissues

Suppose that we have an attribute file `tissue1.out` (like `ecotype1.out` in above) and a value file (like `ecotype0.out` in above). We can similarly generate a classification matrix for tissues.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head extraction/tissue1.out
attr    match   total   ratio   selection
tissue  28412   29769   0.95    1
source_name     14104   21917   0.64    1
organism part   2918    3016    0.97    1
dev_stage       2672    7251    0.37    0
sample_type     407     1914    0.21    0
plant structure 323     504     0.64    1
tissue_type     305     324     0.94    1
dev stage       304     546     0.56    0
organ   292     292     1.00    1

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ head extraction/tissue0.out
value   count   leaf    rosette root    shoot   flower  inflorescence   pollen  anther  seedling        hypocotyl       cotyledon       seed    embryo  endosperm       whole       aerial  OR
leaf    3190    1       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       1
root    2939    0       0       1       0       0       0       0       0       0       0       0       0       0       0       0       0       1
seedlings       2314    0       0       0       0       0       0       0       0       1       0       0       0       0       0       0       0       1
seedling        2088    0       0       0       0       0       0       0       0       1       0       0       0       0       0       0       0       1
leaves  1591    1       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       1
whole seedling  1391    0       0       0       0       0       0       0       0       1       0       0       0       0       0       0       0       1
shoot   1283    0       0       0       1       0       0       0       0       0       0       0       0       0       0       0       0       1
whole seedlings 977     0       0       0       0       0       0       0       0       1       0       0       0       0       0       0       0       1
roots   741     0       0       1       0       0       0       0       0       0       0       0       0       0       0       0       0       1

wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ ../scripts/biosampleClassify.pl extraction/tissue1.out extraction/tissue0.out sel20261001.xml > extraction/sel20261001.tissue
```

Given that we have the normalized log-count-per-million matrix of only col-0 samples saved in tab-delimited text file `sel20261001.nMatrix.col0.TMM`, the following command can be applied for generating portions of tissues extracted from the normalized count matrix.
```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath/extraction$ ../../scripts/matrixSelection.pl sel20261001.tissue sel20261001.nMatrix.col0.TMM sel20261001.nMatrix.col0.TMM leaf rosette root shoot flower inflorescence pollen anther seedling hypocotyl cotyledon seed embryo endosperm whole aerial

wdlin@comp04:/RAID2/R418/20261001_coexDB$ find coexDB20261001/ | perl -ne 'chomp; next if -d "$_"; print "$_\n"' | perl -ne 'chomp; $msg=`head -1 $_`; chomp $msg; @t=split(/\t/,$msg); $cnt=@t; $cnt--; print "$_\t$cnt\n"'
coexDB20261001/ath/sel20261001.ath.col0.TMM.whole       1004
coexDB20261001/ath/sel20261001.ath.col0.TMM.flower      991
coexDB20261001/ath/sel20261001.ath.col0.TMM.embryo      66
coexDB20261001/ath/sel20261001.ath.col0.TMM     18539
coexDB20261001/ath/sel20261001.ath.ler.TMM      452
coexDB20261001/ath/sel20261001.ath.col0.TMM.root        2908
coexDB20261001/ath/sel20261001.ath.col0.TMM.leaf        4157
coexDB20261001/ath/sel20261001.ath.col0.TMM.pollen      219
coexDB20261001/ath/sel20261001.ath.col0.TMM.endosperm   25
coexDB20261001/ath/sel20261001.ath.col0.TMM.cotyledon   202
coexDB20261001/ath/sel20261001.ath.col0.TMM.seed        798
coexDB20261001/ath/sel20261001.ath.col0.TMM.aerial      231
coexDB20261001/ath/sel20261001.ath.col0.TMM.seedling    5685
coexDB20261001/ath/sel20261001.ath.col0.TMM.anther      69
coexDB20261001/ath/sel20261001.ath.col0.TMM.inflorescence       323
coexDB20261001/ath/sel20261001.ath.col0.TMM.rosette     1146
coexDB20261001/ath/sel20261001.ath.col0.TMM.shoot       1212
coexDB20261001/ath/sel20261001.ath.TMM  37808
coexDB20261001/ath/sel20261001.ath.col0.TMM.hypocotyl   274
```
The last command is for numbers of data columns (samples) in the extracted matrixes.

### Sample classification part 4 (optional), iteratively refine attributes & values for selection

In above example, we started from one attrabute, collect *values* of our interests, decide *attributes*, and the use lastly adopted attributes and values for sample classification. Actually the process is flexible. For example, we can use those decided *attributes* to search more *values* for our decision. For example, the following perl oneliner was to use decided *attributes* (in file `ecotype1.out`) to collect more *values* for making decision. Just remember to give an attribute file and a value file for generating a selection matrix.

```
wdlin@comp04:/RAID2/R418/20261001_coexDB/ath$ cat sel20261001.xml | perl -ne 'if($.==1){ open(FILE,"<extraction/ecotype1.out"); while($line=<FILE>){ $line=~s/^\s+|\s+$//g; @s=split(/\t/,$line); $hash{$s[0]}=1 if $s[-1]==1 } close FILE; } chomp; if(/<Attribute attribute_name="(.+?)".*?>(.+?)</){ print "$2\n" if exists $hash{$1} }' | perl -ne 'chomp; $hash{lc($_)}++; if(eof){ print "value\tcount\n"; for $k (sort {$hash{$b}<=>$hash{$a}} keys %hash){ print "$k\t$hash{$k}\n" } }' > extraction/ecotype2.in
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
