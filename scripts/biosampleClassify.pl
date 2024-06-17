#!/usr/bin/env perl

my $usage = "biosampleClassify.pl <attrFile> <valueFile> <biosampleXML>\n";

my $attrFilename = shift or die $usage;
my $valueFilename = shift or die $usage;
my $biosampleFilename = shift or die $usage;

# read attribute file
open(FILE,"<$attrFilename");
$line=<FILE>; # skip first line
while($line=<FILE>){
    chomp $line;
    $line=~s/^\s+|\s+$//g;
    my @s=split(/\t/,$line);
    $attrHash{$s[0]}=0 if $s[-1] eq "TRUE";
} close FILE;

# read value file
open(FILE,"<$valueFilename");
# get target index from first line
$line=<FILE>;
chomp $line;
$line=~s/^\s+|\s+$//g;
my @s=split(/\t/,$line);
print "SRS";
for(my $i=1;$i<@s;$i++){
    $idxTargetHash{$i}=$s[$i];
    push @targets, $s[$i];
    print "\t$s[$i]";
}
print "\n";
# rest lines
while($line=<FILE>){
    chomp $line;
    $line=~s/^\s+|\s+$//g;
    my @s=split(/\t/,$line);
    for(my $i=1;$i<@s;$i++){
        $valHash{lc($s[0])}{$idxTargetHash{$i}}=0 if $s[$i] eq "TRUE";
    }
}
close FILE;

# read biosample xml
open(FILE,"<$biosampleFilename");
while(<FILE>){
    if(/<BioSampleSet>/){
        %flag=();
        %srs =();
    }

    if(/<Id db="SRA">(.+)<\/Id>/){
        $srs{$1}=1;
    }
    
    if(/<Attribute attribute_name="(.+?)".*?>(.+?)</){
        $attr=$1;
        $val=lc($2);
        
        if(exists $attrHash{$attr} && exists $valHash{$val}){
            for $x (keys %{$valHash{$val}}){
                $flag{$x}=1;
            }
        }
    }
    
    if(/<\/BioSampleSet>/){
        for $k (sort keys %srs){
            print "$k";
            for $x (@targets){
                if(exists $flag{$x}){
                    print "\t1";
                }else{
                    print "\t0";
                }
            }
            print "\n";
        }
    }
}
