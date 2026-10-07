#!/usr/bin/env perl

my $usage = "biosampleRetrieveBySRS.pl <listFile> <outSrsBios> <outXML>\n";

my $hintFile = "";
my $maxTry = 3;

# retrieve optional parameter
my @arg_idx=(0..@ARGV-1);
for my $i (0..@ARGV-1) {
    if ($ARGV[$i] eq '-hint') {
        $hintFile=$ARGV[$i+1];
        delete @arg_idx[$i,$i+1];
    }elsif ($ARGV[$i] eq '-maxtry') {
        $maxTry=$ARGV[$i+1];
        delete @arg_idx[$i,$i+1];
    }
}
my @new_arg;
for (@arg_idx) { push(@new_arg,$ARGV[$_]) if (defined($_)); }
@ARGV=@new_arg;

# regular parameters
my $listFilename = shift or die $usage;
my $outSrsBios = shift or die $usage;
my $outXML     = shift or die $usage;

# read list
open(FILE,"<$listFilename");
my @accArr = ();
while(<FILE>){
    chomp;
    s/^\s+|\s+$//g;
    push @accArr, $_;
}
close FILE;

# read hint if specified
my %hintHash = ();
if(length($hintFile)){
    open(FILE,"<$hintFile");
    while(<FILE>){
        chomp;
        my @t=split;
        $hintHash{$t[0]} = $t[1];
    }
    close FILE;
}

# read outSrsBios for existing records for skipping them
my %finished;
open(FILE,"<$outSrsBios");
while(<FILE>){
    chomp;
    my @t=split;
    $finished{$t[0]}=1;
}
close FILE;

# iterate list
open(FILE1,">>$outSrsBios");
open(FILE2,">>$outXML");
for my $acc (@accArr){
    next if exists $finished{$acc};

    print "RETRIEVE: $acc\n";

    $biosAcc = "";
    if(exists $hintHash{$acc}){
        $biosAcc = $hintHash{$acc};
    }
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
