#!/usr/bin/perl

use warnings;
use strict;

if (@ARGV < 2)
{
    print "USAGE: $0 file1.bed file2.bed

From file2.bed this script will randomely
select fragments of the same length as the
fragments in file1.bed.

Example:

file1.bed

chr1    3       5
chr1    7       10
chr1    12      18
chr1    22      24

file2.bed

chr1    1       3
chr1    5       7
chr1    10      12
chr1    18      22
chr1    24      40

$0 file1.bed file2.bed

chr1	1	3	non_peak
chr1	5	7	non_peak
chr1	18	21	non_peak
chr1	24	30	non_peak

";
    exit;
}

my @peaks = ();

open (PEAK, $ARGV[0]);
while (<PEAK>)
{
    chomp;
    if ($_ =~ /\d/)
    {
        my @temp = split '\t', $_;
        push @peaks, $temp[2]-$temp[1];
    }
}
close (PEAK);

open (NOPEAK, $ARGV[1]);
while (<NOPEAK>)
{
    chomp;
    if ($_ =~ /\d/)
    {
        my @temp = split '\t', $_;
        my $non_peak_length = $temp[2]-$temp[1];
        my $non_peak_start  = $temp[1];
        foreach my $peak_length (0..$#peaks)
        {
            if ($peaks[$peak_length] <= $non_peak_length)
            {
                print "$temp[0]\t$non_peak_start\t", $non_peak_start+$peaks[$peak_length], "\tnon_peak\n";
                $peaks[$peak_length] = 1000;
                last;
            }
        }
    }
}
close (NOPEAK);

