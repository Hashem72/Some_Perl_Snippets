#!/usr/bin/perl
use strict;
use warnings;
srand(time);

use lib "/nfs/users/nfs_h/hk3/My_Perl_Scripts/modules";
use Ensembl;
use FileTools;
use SeqTools;
use MathsTools;
use PrintingTools;



my $data_resource = "UW";
my $Rep          = 1;
my $cell_line    = "Hsmm";
my $chr_22 = &FileTools::get_seq_from_fasta_file("/nfs/th_group/hk3/Human_CHRS/chr_22.fasta", "Human_chromosome_22");



my $DHSs_file = "/nfs/th_group/hk3/".$data_resource."_DNaseI_HS/".$cell_line."_For_Paper_Analysis/".$cell_line."_Encode_hotspots.bed";
my $bindind_scores_in_DHS_output_file = "/nfs/th_group/hk3/".$data_resource."_DNaseI_HS/".$cell_line."_For_Paper_Analysis/".$cell_line."_Binding_Scores_in_ENCODE_DHS.txt";





  
my $pwm_file_real = "/nfs/th_group/hk3/".$data_resource."_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_real_tags_rep".$Rep.".txt";
my  $pwm_file_bg   = "/nfs/th_group/hk3/".$data_resource."_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_shifted_tags_rep".$Rep.".txt";


my $pwm_real       = &FileTools::read_in_a_matrix_from_a_file($pwm_file_real);
 my $pwm_bg        =  &FileTools::read_in_a_matrix_from_a_file($pwm_file_bg);



&get_binding_prob_of_all_features_in_a_bed_file($DHSs_file, $bindind_scores_in_DHS_output_file, $chr_22, $pwm_real, $pwm_bg, 5, 14, 100, 0);

exit;





######################################  subroutines ##############################


sub get_binding_prob_of_all_features_in_a_bed_file($$$$$$$$$){
    my $input_bed_file      = shift;
    my $output_file         = shift;
    my $chr_sequence        = shift;
    my $pwm_R               = shift;
    my $pwm_BG              = shift;
    my $start_pos_in_pwm    = shift; 
    my $window_length       = shift; 
    my $length_threshold    = shift; # sequences shorter than this will be ignored
    my $random_or_real      = shift; # for real 0 for random

    
    if($random_or_real eq 0){
	$output_file  = &make_output_file_name_for_random_seqs($output_file);
    }

    open(IN, "$input_bed_file") or die "cannot open file $input_bed_file to read the features\n";
    open(OUT, ">$output_file") or die "cannot open file $input_bed_file  to write into\n";

   
    while(<IN>){
    	chomp $_;
    	my @fields = split("\t",$_);
    	my $chr    = $fields[0];
    	my $start  = $fields[1];
    	my $end    = $fields[2];
    	my $source = $fields[3];
    	my $score  = $fields[4];
    	my $strand  = $fields[5];

	#ignore  a dhs if it contains Ns
	my $feature_length  = $end - $start;
	my $dhs  = uc(substr($chr_sequence, $start, $feature_length));
	next if ($dhs =~ m/N/);


	#ignore it if it is short
	next if ($feature_length < $length_threshold);

	#try getting a random sequence that doesnt contain Ns.
	my  $one_sequence ;
	if($random_or_real == 0){
	    my $sequence_contains_N = "true";
	    while($sequence_contains_N eq "true"){
		$start =  &get_a_random_locus_in_active_arm_of_chr_22();
		$one_sequence  = uc(substr($chr_sequence, $start, $feature_length));
		if( $one_sequence !~ m/N/ ){
		    $sequence_contains_N = "false";
		}
	    }
	}
	
    	
    	
    	  $one_sequence  = uc(substr($chr_sequence, $start, $feature_length));
    	my @binding_socres_for_one_seq = &get_binding_scores_for_one_seq($one_sequence, $pwm_R, $pwm_BG, $start_pos_in_pwm, $window_length);
    	# print out binding scores
    	foreach my $s(@binding_socres_for_one_seq){
    	    print OUT
    		$s.',';
    	}
    	print OUT
	    "\n";
	
     }
    
     close(IN);
     close(OUT);


}#get_binding_prob_of_all_features_in_a_bed_file#





sub get_binding_scores_for_one_seq($$$$$){
    my $sequence         = shift;
    my $pwm_R             = shift;
    my $pwm_BG            = shift;
    my $start_pos_in_pwm  = shift; 
    my $window_length     = shift; 

    my @binding_scores;

    my $seq_length  = length($sequence);
    my $upper_bound = $seq_length - $window_length;
    $sequence      = uc($sequence);
    if($sequence =~ m/N/){
	print "sequences was containing some Ns and therefore it was ignored\n";
    }
    else{
	for(my $p=0; $p<$upper_bound; $p++){
	    my $one_stretch = uc(substr($sequence, $p, $window_length));
	    my $one_stretch_revcomp =  &SeqTools::get_revcomp($one_stretch);


#get binding_score for pve starnd:
	    my $prob_seq_given_real_pwm_pve = &prob_seq_given_real_pwm($pwm_R, $one_stretch, $start_pos_in_pwm);
	    my $prob_seq_given_bg_pwm_pve   = &prob_seq_given_real_pwm($pwm_BG, $one_stretch, $start_pos_in_pwm);
	    my $one_log_ratio_pve           = $prob_seq_given_real_pwm_pve - $prob_seq_given_bg_pwm_pve;
#get binding_score nve starnd:
	    my $prob_seq_fiven_real_pwm_nve = &prob_seq_given_real_pwm($pwm_R, $one_stretch_revcomp, $start_pos_in_pwm);
	    my $prob_seq_given_bg_pwm_nve   = &prob_seq_given_real_pwm($pwm_BG, $one_stretch_revcomp, $start_pos_in_pwm);
	    my $one_log_ratio_nve           =  $prob_seq_fiven_real_pwm_nve - $prob_seq_given_bg_pwm_nve ;
	    
# get max of scores:
	    my $max_score = &get_max($one_log_ratio_pve, $one_log_ratio_nve);
	    my $rounded_score = sprintf "%.5f",$max_score ;	    
	    push(@binding_scores, $rounded_score)
	    
	}
    }
    return @binding_scores;
				   }#get_binding_scores_for_one_seq#




sub get_max($$){
    my $first_numb = shift;
    my $second_numb = shift;
    if($first_numb >= $second_numb){
	return $first_numb;
    }
    else{  
	return $second_numb;
    }
}#get_max#


sub prob_seq_given_real_pwm($$$$){
    my $pwm     =  shift;
    my $seq     = shift;
    my $start   = shift;
    

    my $seq_length = length($seq);
    if(  ($seq_length+$start)>50  ){
	print
	    "out of boundary problem!\n";
    }
    my $sum  =     0;
    for (my $p = 0; $p < $seq_length; $p++) {
        my $one_base = substr($seq, $p, 1);
	my $one_weight;
	if ($one_base eq "A") {
	    $one_weight = $pwm->[$p+$start][0];
	}
	elsif($one_base eq "C"){
	    $one_weight = $pwm->[$p+$start][1];
	}
	elsif($one_base eq "G"){
	    $one_weight = $pwm->[$p+$start][2];
	}
	elsif($one_base eq "T"){
	    $one_weight = $pwm->[$p+$start][3];
	}
	else{
	    die "Unkown base: $one_base found in sequence : $seq\n";
	}
	$sum = $sum + $one_weight;
    }
    return $sum;
}#prob_seq_given_real_pwm#


sub make_output_file_name_for_random_seqs($){
    my $output_file_name  = shift;

    my @parts = split(/\./, $output_file_name);
    my $file_name_without_extension = $parts[0];
    my $new_output_file  = $file_name_without_extension."_From_Randomly_Pickes_Sequences.txt";
    return $new_output_file;
}#make_output_file_name_for_random_seqs#


sub get_a_random_locus_in_active_arm_of_chr_22(){
    
    my $min    = 16500000;
    my $Length = 50000000;
    my $range  = $Length  - $min;
    my $random_integer = int(rand($range)) + $min; 
    return  $random_integer;
}#get_a_random_locus_in_active_arm_of_chr_22#




