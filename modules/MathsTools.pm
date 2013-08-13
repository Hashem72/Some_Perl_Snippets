package MathsTools;

use strict;
use warnings;


########################    subroutines ########################## 

sub get_an_array_of_random_integers($$$){
    my $min       = shift;
    my $max       = shift;
    my $number    = shift;
    
    my @rand_integers = ();
    
    my $range = $max - $min;
    my $size = 0;
    while($size < $number){
	my $one_rand = int(rand($range)) + $min;
	push(@rand_integers,$one_rand);
	$size = @rand_integers;
    }
    return @rand_integers;
		       }#get_an_array_of_random_integers#

sub log2($){
    my $n = shift;

    my $L = (log($n))/(log(2));
    return $L;
}#log2#


sub pfm2pwm($$){
    my $pfm      = shift;
    my $nucliotides_prior_prob = shift;
    
    my $pwm = [];
    
    my $number_of_seqs_used = &get_get_total_number_of_seqs_used_in_freq_matrix($pfm);
    print
	"number_of_seqs_in_used is $number_of_seqs_used \n";
	
    my $number_of_rows = @{$pfm};
    my $number_of_cols = @{$pfm->[0]};
    
    for (my $r = 0; $r < $number_of_rows; $r++) {
	for (my $c = 0; $c < $number_of_cols; $c++) {
	    my $numerator = ($pfm->[$r][$c])+($nucliotides_prior_prob->[$c]);
	    my $denomenator = ($number_of_seqs_used + 1) * ($nucliotides_prior_prob->[$c]);
	    my $ration = $numerator/$denomenator;
	    $pwm->[$r][$c] = &log2($ration);
	}
    }
    return $pwm;
}#pfm2pwm#


sub pssm2transfac_format($$$$){
    my $pssm      = shift;
    my $id        = shift;
    my $species   = shift;
    my $output_file = shift;

    open(OUT, ">$output_file") or die "cannot open file $output_file to write into it!\n";
    print OUT "ID"." ".$id."\n";
    print OUT "BF"." ".$species."\n";
    print OUT "P0\tA\tC\tG\tT\n";
    

    my $number_of_rows = @{$pssm};
    my $number_of_cols = @{$pssm->[0]};
    my $counter        = 0; 
    for(my $r=0; $r < $number_of_rows; $r++){
	if($r<10){
	  $counter = "0".$r;  
	}
	else{
	    $counter = $r;
	}
	print OUT $counter."\t";
	for(my $c =0; $c < $number_of_cols; $c++){
	    print OUT $pssm->[$r][$c]."\t";
	}
	print OUT "\n";
    }
    print OUT "XX\n";
    print OUT "//";
    close(OUT);
}#pssm2transfac_format#



sub get_get_total_number_of_seqs_used_in_freq_matrix($){
    my $freq_matrix   = shift;

    my $sum = 0;
    for (my $c = 0; $c < 4; $c++) {
	$sum = $sum + $freq_matrix->[0][$c];
    }
    return $sum;
}#get_get_total_number_of_seqs_used_in_freq_matrix#

sub get_prob_seq_given_pwm($$$$){
    my $pwm            = shift;
    my $seq           = shift;
    my $start         = shift;
    my $subseq_length = shift;

    
    my $end_point  = $start + $subseq_length ;
    my $seq_length = length($seq);
    if($end_point > $seq_length){
    
	die "end of subsequence greater than end of sequence(tag)!\n";
    }
    my $sum        = 0;
    for (my $p = $start; $p < $end_point; $p++) {
	my $one_base = substr($seq, $p, 1);
	my $one_weight;
	if ($one_base eq "A") {
	    $one_weight = $pwm->[$p][0];
	}
	elsif($one_base eq "C"){
	    $one_weight = $pwm->[$p][1];
	}
	elsif($one_base eq "G"){
	    $one_weight = $pwm->[$p][2];
	}
	elsif($one_base eq "T"){
	    $one_weight = $pwm->[$p][3];
	}
	else{
	    die "Unkown base: $one_base found in sequence : $seq\n";
	}
	$sum = $sum + $one_weight;
    }
    return $sum;
}#get_prob_seq_given_pwm#

1;




