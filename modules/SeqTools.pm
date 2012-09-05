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




1;
