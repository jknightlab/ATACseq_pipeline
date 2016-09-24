#!usr/bin/perl


use strict;

open (IN, $ARGV[0]) or die ("Couldn't open file: $! \n");

my $header=$ARGV[1];

open (HEAD, ">$header.txt") or die ("Couldn't open file: $! \n");

my $bin170orig="bin170orig.txt";
my $binGreaterorig="binGreaterorig.txt";


open (OUT1, ">$bin170orig.sam") or die ("Couldn't open file: $! \n");
open (OUT2, ">$binGreaterorig.sam") or die ("Couldn't open file: $! \n");

while (<IN>) {
    chomp;
    if (/^@/){                                                  # removes the headers of the SAM file
    
    print HEAD "$_\n";
    #print "pushed header line\n";
    next;
    }
    close HEAD;
    
    my ($qname, $flag, $chr, $start, $mapq, $length, $mrnm, $mpos, $isize, $seq, $qual, $tag, $vtype, $value) = split (/\t+/);  # get the columns from SAM
    my $isize_changed=abs($isize);    

    if ( $isize_changed >= 50 and $isize_changed <= 170){
        print OUT1 "$qname\t$flag\t$chr\t$start\t$mapq\t$length\t$mrnm\t$mpos\t$isize\t$seq\t$qual\t$tag\t$vtype\t$value\n";
        }      
    if ( $isize_changed >= 171 and $isize_changed <= 650){
        print OUT2 "$qname\t$flag\t$chr\t$start\t$mapq\t$length\t$mrnm\t$mpos\t$isize\t$seq\t$qual\t$tag\t$vtype\t$value\n";
        }

}


close IN;
close OUT1;
close OUT2;


