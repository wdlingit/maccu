#!/usr/bin/env perl

$| = 1;

my $usage = "duplicateDetect.pl <count matrix TSV>\n";

my $matrixFilename = shift or die $usage;

# read matrix file, header (samples)
my %indexHash;      # $indexHash{index} = sample

open(FILE,"<$matrixFilename");
my $header=<FILE>;
chomp $header;
my @t=split(/\s+/,$header);
for($i=1;$i<@t;$i++){
    $indexHash{$i} = $t[$i];
}

# read matrix file, count part
my %valueHoA;       # $valueHoA{sample} = [ counts ]

print "Reading matrix\n";
while(<FILE>){
    @t = split;
    for($i=1;$i<@t;$i++){
        push @{$valueHoA{$indexHash{$i}}}, $t[$i];
    }
}
close FILE;

# count sum as hash
my %hashSample;     # $hashSample{hashcode} = [ samples ]

print "Compute hash\n";
for my $s (keys %valueHoA){
    my $arr = $valueHoA{$s};
    my $sum = 0;
    
    for($i=0;$i<@{$arr};$i++){
        $sum += $$arr[$i];
    }
    
    push @{$hashSample{$sum}}, $s;
}

# actual check
my %equalPair;      # $equalPair{$s1}{$s2} = 1
my %dupSamples;     # $dupSamples{sample} = 1

print "Compare\n";
for my $hash (keys %hashSample){
    my @arr = @{$hashSample{$hash}};
    
    if(@arr>1){ # size at least 2
        for(my $i=0;$i<@arr-1;$i++){
            for(my $j=$i+1;$j<@arr;$j++){
                if(equalCheck($arr[$i],$arr[$j])){
                    $equalPair{$arr[$i]}{$arr[$j]} = 1;
                    $equalPair{$arr[$j]}{$arr[$i]} = 1;
                    $dupSamples{$arr[$i]} = 1;
                    $dupSamples{$arr[$j]} = 1;
                }
            }
        }
    }
}

# report
print "Report\n";
while(keys %dupSamples > 0){
    my @arr = sort keys %dupSamples;
    
    # get first
    my $first = $arr[0];
    # report it and all equal samples
    my @rest = sort keys %{$equalPair{$first}};
    
    print "DUP: $first\t".join("\t",@rest)."\n";
    
    # remove them
    push @rest,$first;
    
    for my $s (@rest){
        delete $dupSamples{$s};
    }
}

sub equalCheck {
    my $sample1 = shift;
    my $sample2 = shift;
    
    my $arr1 = $valueHoA{$sample1};
    my $arr2 = $valueHoA{$sample2};

    my $ans=1;
    for(my $i=0;$i<@{$arr1};$i++){
        if($$arr1[$i] != $$arr2[$i]){
            $ans=0;
            last;
        }
    }
    
    return $ans;
}
