package SeqTools;
use strict;
use warnings;

##################### subroutines ##################

sub get_revcomp($){
    my $dna_seq  = shift;
    
    my $revcom_seq = reverse $dna_seq;
    $revcom_seq =~ tr/NACGTnacgt/NTGCAntgca/;
    return $revcom_seq;
		}#GET_REVCOMP#

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

sub get_freq_ncts($){
    my $sequence = shift;

    my %freqs = ();

    my $seq_length  = length($sequence);
    #get counts
    my $count_A     = ($sequence =~ tr/Aa//);
    my $count_C     = ($sequence =~ tr/Cc//);
    my $count_G     = ($sequence =~ tr/Gg//);
    my $count_T     = ($sequence =~ tr/Tt//);
    my $all_nct     = ($sequence =~ tr/ACGTacgt//);
    my $non_nct     = $seq_length - $all_nct;

    #get fractions
    my $frac_A   = sprintf("%.2f", $count_A/$seq_length); 
    my $frac_C   = sprintf("%.2f", $count_C/$seq_length);
    my $frac_G   = sprintf("%.2f", $count_G/$seq_length);
    my $frac_T   = sprintf("%.2f", $count_T/$seq_length);
    my $frac_non = sprintf("%.2f", $non_nct/$seq_length); 
    
    #assign to a hash
    $freqs{A} = $frac_A;
    $freqs{C} = $frac_C;
    $freqs{G} = $frac_G;
    $freqs{T} =  $frac_T;
    $freqs{non_nct} = $frac_non;

    return %freqs;
}#get_freq_ncts#



1;
