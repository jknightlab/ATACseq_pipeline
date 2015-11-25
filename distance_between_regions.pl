#!/usr/bin/perl

use warnings;
use strict;


# copyright Irina Pulyakhina
# jknight group, Oxford University
# 28 Sept 2015
 

my $help_message = "\nUSAGE: $0 input.bed\n
where input.bed is the input file in a bed format. The
requirement for the input file is for the first column
to contain the chromosome name; for the second column
to contain the start coordinate of a region; for the
third coordinate to contain the end coordinate of a
region. The script will take each two subsequent regions
on a chromosome and calculate the distance between them.\n\n";

if ((@ARGV < 1) or ($ARGV[0] eq "-h") or ($ARGV[0] eq "--help"))
{
    print $help_message;
    exit;
}

my $chromosome_counter = "";
my $previous_peak_end = 0;
my $next_peak_start = 0;
my @distances = ();

open (BED, $ARGV[0]);
while (<BED>)
{
    chomp;
    my @string = split '\t', $_; # splitting each string using tab as a delimiter

    # $chromosome_counter is defined starting from the second iteration through the file
    # $chromosome_counter is the chromosome name of the previous line
    if ($chromosome_counter)
    {
        # comparing chromosome names of the previous and the current line
        if ($chromosome_counter eq $string[0]) 
        {
            # $previous_peak_end will be defined starting from the second iteration through the file
            if ($previous_peak_end)
            {
                # #string[1] is a potential start coordinate of the next peak
                # checking that it is bigger than the end coordinate of the previous peak
                # this might not be the case when the file is not sorted
                if ($string[1] >= $previous_peak_end)
                {
                    $next_peak_start = $string[1];
                    # calculating the distance between the peaks
                    my $dist = $next_peak_start - $previous_peak_end;
                    # adding the distance to the array with the rest of the distances
                    push (@distances, $dist);
                }
            }
        }
    }
    $chromosome_counter = $string[0];
    $previous_peak_end = $string[2];
}
close (BED);

# printing the distances between the peaks
foreach my $d (sort {$a <=> $b} @distances)
{
    print "$d\n";
}
