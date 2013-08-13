#!/usr/bin/perl
use strict;
use warnings;

use lib "/nfs/users/nfs_h/hk3/My_Perl_Scripts/modules";
use Ensembl;
use FileTools;
use SeqTools;
use MathsTools;
use PrintingTools;


my $tag_resource = "UW";
my $Rep          = 1;
#my $cell_line    = "H1hesc";
my $cell_line = $ARGV[0];



if  (defined($ARGV[0])){
    print
	"you have asked me to look at ". $cell_line."\n";
	
}
else{
    die
	"you should pass a cell line name!\n";

}

chomp($cell_line);
my $start_of_stretch = 5; # this is the start of the stretch of the seq you would like to be considered 
my $length_of_stretch = 15;  # the length of the stretch you want to be considered

my $start_offset;
my $end_offset;
my $tags_bed_file;
my $output_file_for_motif_scores;
my $pwm_file_real;
my $pwm_file_bg ;

if($tag_resource eq "UW"){
    $start_offset =  10;
    $end_offset   = 4;


    $tags_bed_file = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."/wgEncodeUwDnase".$cell_line."AlnRep1_chr22.bed";
     $output_file_for_motif_scores  = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."_For_Paper_Analysis/motif_socres_wgEncodeUwDnase".$cell_line."AlnRep1_chr22_s".$start_of_stretch."_l".$length_of_stretch.".txt";
    
     $pwm_file_real = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_real_tags_rep".$Rep.".txt";
     $pwm_file_bg   = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_shifted_tags_rep".$Rep.".txt";
}
elsif($tag_resource eq "Duke"){
    $start_offset =  15;
    $end_offset   = 15;
    $tags_bed_file = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."/wgEncodeOpenChromDnase".$cell_line."AlnRep1_chr22.bed";
    $output_file_for_motif_scores  = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."_For_Paper_Analysis/wgEncodeOpenChromDnase".$cell_line."AlnRep1_chr22_s".$start_of_stretch."_l".$length_of_stretch.".txt";
    $pwm_file_real = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_real_tags.txt";
    $pwm_file_bg   = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_shifted_tags.txt";

}
else{
    die "Unknown tag resource\n";
}

 my $chr_22 = &FileTools::get_seq_from_fasta_file("/nfs/th_group/hk3/Human_CHRS/chr_22.fasta", "Human_chromosome_22");

 my $pwm_real      = &FileTools::read_in_a_matrix_from_a_file($pwm_file_real);
 my $pwm_bg        =  &FileTools::read_in_a_matrix_from_a_file($pwm_file_bg);


&get_binding_scores($chr_22, $pwm_real, $pwm_bg, $output_file_for_motif_scores, 5, 1, 15, 16500000);

exit;

############################################# subroutines ###############################################


sub get_binding_scores($$$$$$$$){
    my $chr_seq              = shift;
    my $pwm_R                = shift;
    my $pwm_BG               = shift;
    my $output_file          = shift;
    my $start_pos_in_pwm     = shift; # my pws is of length 50, howver the only informative part of it is [5,20] ie stretch of 15 sitting at pos 5, this start_pos_in_pwm is to give this position
    my $step_size            = shift || 1; # step_size for  scanning the chr. if not given take it 1bp
    my $window_length        = shift;
    my $scan_start_pos       = shift || 0; # the first 16m of chr22 is the dead arm so there is no need to scan. this is to tell from where start to scan. if not give it is zero

    my $chr_length      =   length($chr_seq);
    my $score_for_pve   =   0; 
    my $score_for_nve   =   0;
    my $max_score       =   0;
    my  $dummy_counter  =   0;
    open(OUTPUT, ">$output_file") or die "cannot open file $output_file to write in!\n ";

    for(my $p =$scan_start_pos; $p< $chr_length; $p = $p+ $step_size){
	# if( ($p % 1000000) == 0){
	#     print $p;
	#     }
	$dummy_counter++;
	my $one_stretch = uc(substr($chr_seq, $p, $window_length));
	my $one_stretch_revcomp =  &SeqTools::get_revcomp($one_stretch);
	if($one_stretch =~ m/N/){
	    my $max_score  = -1000;
	    print OUTPUT
	    $max_score."\n";

	    next;
	}
	
#get binding_score for pve starnd:
	my $prob_seq_given_real_pwm_pve = &prob_seq_given_real_pwm($pwm_R, $one_stretch, $start_pos_in_pwm);
	my $prob_seq_given_bg_pwm_pve   = &prob_seq_given_real_pwm($pwm_BG, $one_stretch, $start_pos_in_pwm);
	my $one_log_ratio_pve           = $prob_seq_given_real_pwm_pve - $prob_seq_given_bg_pwm_pve;
#get binding_score nve starnd:
	my $prob_seq_fiven_real_pwm_nve = &prob_seq_given_real_pwm($pwm_R, $one_stretch_revcomp, $start_pos_in_pwm);
        my $prob_seq_given_bg_pwm_nve   = &prob_seq_given_real_pwm($pwm_BG, $one_stretch_revcomp, $start_pos_in_pwm);
	my $one_log_ratio_nve           =  $prob_seq_fiven_real_pwm_nve - $prob_seq_given_bg_pwm_nve ;

 # get max of scores:
	$max_score = &get_max($one_log_ratio_pve, $one_log_ratio_nve);
	my $rounded_score = sprintf "%.5f",$max_score ;
	print OUTPUT
	    $rounded_score."\n";
	    
    }
    close(OUTPUT);
}#get_binding_socres#

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

#####################################################################################################








