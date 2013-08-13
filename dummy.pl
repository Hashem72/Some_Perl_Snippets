#!/usr/bin/perl
use strict;
use warnings;

use lib "/nfs/users/nfs_h/hk3/My_Perl_Scripts/modules";
use Ensembl;
use FileTools;
use SeqTools;
use MathsTools;
use PrintingTools;



my $dummy_seq = "ggggggggggaaaaattttt";
my $gc_content = &get_gc_content($dummy_seq);
print
    "gc_content = $gc_content\n";








sub get_gc_content($){
    my $dna_seq              = shift;
    $dna_seq                 = uc($dna_seq);
    my $number_of_known_nts  = 0;
    while( $dna_seq =~ m/[ACGT]/gi){
	$number_of_known_nts++;
    }
    my $gc_count   = 0;
    while( $dna_seq =~ m/[GC]/gi){
	$gc_count++;
    }
    my $gc_fraction = $gc_count/$number_of_known_nts;
    my ($gc_dec) = $gc_fraction =~ /(\S{0,6})/;
    return $gc_dec;
}#get_gc_content#

