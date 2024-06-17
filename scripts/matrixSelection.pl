#!/usr/bin/env perl

my $usage = "matrixSelection.pl <selMatrix> <sourceMatrix> <outPrefix> [<selTarget>]+\n";

my $selMtxFilename = shift or die $usage;
my $srcMtxFilename = shift or die $usage;
my $outPrefix = shift or die $usage;

die $usage if @ARGV<1;

my %target;
for my $x (@ARGV){
    $target{$x}=0;
}

# read selection matrix
open(FILE,"<$selMtxFilename");
# read first line, get index of targets
$line=<FILE>;
chomp $line;
@s=split(/\t/,$line);
$headerCnt = @s;
for(my $i=0;$i<@s;$i++){
    if(exists $target{$s[$i]}){
        $idxTargetHash{$i}=$s[$i];
        $target{$s[$i]}=1;
    }
}
for my $x (sort keys %target){
    if($target{$x}==0){
        die "non-existing target '$x' in $selMtxFilename\n";
    }
}
# read rest lines, collect row names for targets
my %sampleTargetHash = ();
while($line=<FILE>){
    chomp $line;
    @s=split(/\t/,$line);
    die "token count inconsistent in $selMtxFilename at row '$s[0]...'\n" if @s != $headerCnt;
    for(my $i=1;$i<@s;$i++){
        if(exists $idxTargetHash{$i} && $s[$i]){
            $sampleTargetHash{$s[0]}{$idxTargetHash{$i}}=1;
        }
    }
}
close FILE;

# prepare file handlers
for my $x (keys %target){
    open($outFILE{$x},">$outPrefix.$x");
}

# read source matrix
open(FILE,"<$srcMtxFilename");
# read first line
$line=<FILE>;
chomp $line;
@s=split(/\t/,$line);
$headerCnt = @s;
for(my $i=0;$i<@s;$i++){
    if($i==0){ # print to all files
        for my $x (keys %target){
            print {$outFILE{$x}} "$s[$i]";
        }
    }else{
        if(exists $sampleTargetHash{$s[$i]}){
            $idxSampleHash{$i}=$s[$i];
            for $x (keys %{$sampleTargetHash{$s[$i]}}){
                print {$outFILE{$x}} "\t$s[$i]";
            }
        }
    }
}
for my $x (keys %target){
    print {$outFILE{$x}} "\n";
}
# rest lines
while($line=<FILE>){
    chomp $line;
    @s=split(/\t/,$line);
    die "token count inconsistent in $srcMtxFilename at row '$s[0]...'\n" if @s != $headerCnt;
    for(my $i=0;$i<@s;$i++){
        if($i==0){ # print to all files
            for my $x (keys %target){
                print {$outFILE{$x}} "$s[$i]";
            }
        }else{
            if(exists $idxSampleHash{$i}){
                for $x (keys %{$sampleTargetHash{$idxSampleHash{$i}}}){
                    print {$outFILE{$x}} "\t$s[$i]";
                }
            }
        }
    }
    for my $x (keys %target){
        print {$outFILE{$x}} "\n";
    }
}
close FILE;

# close file handlers
for my $x (keys %target){
    close $outFILE{$x};
}
