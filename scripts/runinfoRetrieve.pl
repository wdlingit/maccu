#!/usr/bin/env perl

my $usage = "runinfoRetrieve.pl <srrList> <processed> <outCSV>\n";

my $listFilename = shift or die $usage;
my $outSrr = shift or die $usage;
my $outCSV = shift or die $usage;

my $chunkSize = 100;
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

# read processed for existing records for skipping them
my %finished;
open(FILE,"<$outSrr");
while(<FILE>){
    @t=split;
    $finished{$t[0]}=1;
}
close FILE;

# collect acc to be processed
my @queryList = ();
for my $acc (@accArr){
    push @queryList, $acc if not exists $finished{$acc};
}

# iterate list
open(FILE1,">>$outSrr");
open(FILE2,">>$outCSV");

while(@queryList){
    # prepare chunk
    my @chunk = splice(@queryList, 0, $chunkSize);
    my $chunkStr = join(",",@chunk);
    my %chunkHash = ();
    for my $x (@chunk){
        $chunkHash{$x} = 1;
    }
    
    print "RETRIEVE: @chunk\n";
    $msg="";
    for($tryNum=0; $tryNum<$maxTry && length($msg)==0; $tryNum++){
        $msg = `efetch -db sra -id $chunkStr -format runinfo`;
        chomp $msg;
        $msg=~s/^\s+|\s+$//g;
    }
    if(length($msg)>0){
        print FILE2 "$msg\n";

        # output processed acc from retrieved message
        foreach my $line (split(/\n/,$msg)){
            chomp $line;
            my @t=split(/,/,$line);
            print FILE1 "$t[0]\n" if exists $chunkHash{$t[0]};
        }
    }
    sleep 1;
}
close FILE1;
close FILE2;
