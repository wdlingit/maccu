#!/usr/bin/env perl

my $usage = "biosampleRetrieve.pl <listFile> <outSrsBios> <outXML>\n";

my $listFilename = shift or die $usage;
my $outSrsBios = shift or die $usage;
my $outXML     = shift or die $usage;
my $skipTo = shift;

my $maxTry = 3;

# read list
open(FILE,"<$listFilename");
my @accArr = ();
while(<FILE>){
    chomp;
    s/^\s+|\s+$//g;
    push @accArr, $_;
}
close FILE;

# read outSrsBios for existing records for skipping them
my %finished;
open(FILE,"<$outSrsBios");
while(<FILE>){
    @t=split;
    $finished{$t[0]}=1;
}
close FILE;

# iterate list
my $startIteration = 0;
if(not defined $skipTo){
    $startIteration = 1;
}
open(FILE1,">>$outSrsBios");
open(FILE2,">>$outXML");
for my $acc (@accArr){
    next if exists $finished{$acc};

    print "RETRIEVE: $acc\n";

    $biosAcc = "";
    for($tryNum=0; $tryNum<$maxTry && length($biosAcc)==0; $tryNum++){ # ATTEMPT 1
        $biosAcc = `esearch -db sra -query $acc | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -if Identifiers -contains $acc -block Id -if \@db -equals BioSample -element Id`;
        chomp $biosAcc;
        $biosAcc=~s/^\s+|\s+$//g;
    }
    for($tryNum=0; $tryNum<$maxTry && length($biosAcc)==0; $tryNum++){ # ATTEMPT 2
        if(length($biosAcc)==0){
            $biosAcc = `esearch -db sra -query $acc | efetch -format docsum | xtract -pattern ExpXml -element Biosample`;
            chomp $biosAcc;
            $biosAcc=~s/^\s+|\s+$//g;

            my @arr = split(/\s+/,$biosAcc);
            my %hash;
            $hash{$_}++ for @arr;
            for my $x (keys %hash){
                $biosAcc = $x if ($hash{$x}/@arr)>0.9;
            }
        }
    }
    
    $msg="";
    for($tryNum=0; $tryNum<$maxTry && length($msg)==0 && length($biosAcc)>0; $tryNum++){
        $msg = `esearch -db biosample -query $biosAcc | efetch -format xml`;
        chomp $msg;
        $msg=~s/^\s+|\s+$//g;
    }

    if(length($msg)>0){
        print FILE2 "$msg\n";
        print FILE1 "$acc\t$biosAcc\n";
    }
}
close FILE1;
close FILE2;
