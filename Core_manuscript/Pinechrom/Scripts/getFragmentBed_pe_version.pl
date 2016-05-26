#!/usr/bin/env perl
###################################
# Author:  AC
# Email: acortes @ well
# Date: Wed Apr  9 10:48:12 2014
###################################
use strict;
use warnings;
###################################
# Instead of extracting the mapping intervals of each read/mate, 
# the script will extract the biological fragment position and 
# output it in bed format.
###################################

# check parameters
if (@ARGV !=2){
    usage();
    exit(1);
}

# Read parameters
my ($in,$out) = @ARGV;

if(! -e $in){
    print "Input bam '$in' is not exists\n";
    usage();
    exit(1);
}

###### Main ############
my $rand=rand();
my $tmpfile = "$rand.bed";

# print $rand,"\n";

# 1) Will first invoke bamToBed
## use the -bedpe flag to get complete fragments. the bam file needs to be sorted by queryname.
`/apps/well/bedtools/2.24.0/bamToBed -bedpe -i $in >$tmpfile`;

# 2) Read the tmp file to get the fragment positions

open(IN,"$tmpfile") or die $!;
my %hash;

while(<IN>){
    chomp;
    # no chrm Reads
    my($chr1,$start1,$end1,$chr2,$start2,$end2,$name,$score,$strand1,$strand2) = split "\t";

    if( $chr1 ne $chr2 ) {
        die "Mates map to different chromosomes!\n";
    }

    push @{$hash{$name}->{$chr1}->{interval}},$start1,$start2,$end1,$end2;

}
close IN;

# 3) Write out
open(OUT,">$out") or die $!;
foreach my $fragment (keys %hash){
    foreach my $chr (keys %{$hash{$fragment}}){
        my @aa = @{$hash{$fragment}->{$chr}->{interval}};
        @aa = sort {$a <=>$b } @aa;

        my $strand = "+";

        my $frag_size = $aa[$#aa] - $aa[0];

        print OUT "$chr\t$aa[0]\t$aa[$#aa]\t$fragment\t$frag_size\t$strand\n";

        # print OUT "$chr\t$aa[0]\t$aa[$#aa]\t$fragment\t$aa[$#aa]-$aa[0]\t$hash{$fragment}->{$chr}->{strand}\n";
        # print OUT join "\t",($chr,$aa[0],$aa[$#aa],$fragment,$aa[$#aa]-$aa[0],$hash{$fragment}->{$chr}->{strand}."\n");
    }
}
close OUT;

END{
    if(defined ($tmpfile) && -e $tmpfile){
        unlink $tmpfile;
    }
}
sub usage{
    print <<EOF;
Usage: getFragmentBed.pl in.bam out.bed
-----------------------------------------------------------
Instead of extracting the mapping intervals of each read/mate, 
the script will extract the biological fragment position and 
output it in bed format.
-----------------------------------------------------------
in.bam   -- input file in bam format
out.bed  -- output file 

EOF

}
