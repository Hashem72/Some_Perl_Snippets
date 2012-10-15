#!/usr/bin/perl
use strict;
use warnings;

use lib "/nfs/users/nfs_h/hk3/My_Perl_Scripts/modules";
use Ensembl;
use FileTools;
use SeqTools;
use MathsTools;
use PrintingTools;
my $data_centre = "UW";
my $cell_line = "H1hesc";

my $prior_probs     = [0.24, 0.26, 0.26, 0.24];
my $tag_lenghts;
my $bed_file;
my $pwm_real_tags_file;
my $pwm_shifte_tags_file;
if ($data_centre eq "UW") {
     $tag_lenghts     = 36;  
     $bed_file        = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."/wgEncodeUwDnase".$cell_line."AlnRep1_chr22.bed";
     $pwm_real_tags_file  = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_real_tags.txt";
     $pwm_shifte_tags_file  = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_shifted_tags.txt";
}
elsif($data_centre eq "Duke"){
    
    $tag_lenghts     = 20;
    $bed_file        = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."/wgEncodeOpenChromDnase".$cell_line."AlnRep1_chr22.bed";
    $pwm_real_tags_file  = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_real_tags.txt";
    $pwm_shifte_tags_file  = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."_For_Paper_Analysis/pwm_shifted_tags.txt";
}
else{
    die "Unknown data centre!\n";
}





my %unique_features = &FileTools::GET_UNIQUE_FEATURES_FROM_A_BED_FILE_V2($bed_file);
my $unique_features = keys %unique_features;
print
    "found ". $unique_features . " unique tags\n";

my $Human_chr_22             = &Ensembl::GET_CHROMOSOME_SEQUENCE("Human","22","false");
my $freq_matrix_real_tags    = &get_position_frequency_matrix($bed_file, 10, 4, 0, $Human_chr_22, $tag_lenghts);
my $pwm_real_tags            = &MathsTools::pfm2pwm($freq_matrix_real_tags, $prior_probs);


my $freq_matrix_shifted_tags =  &get_position_frequency_matrix($bed_file, 10, 4, 40, $Human_chr_22, $tag_lenghts);

my $pwm_shifted_tags         = &MathsTools::pfm2pwm($freq_matrix_shifted_tags, $prior_probs );



#&MathsTools::get_get_total_number_of_seqs_used_in_freq_matrix($freq_matrix_real_tags);
#&PrintingTools::print_out_a_matrix($pwm_real_tags );
&PrintingTools::print_out_a_matrix_into_a_file($pwm_real_tags, $pwm_real_tags_file, ",");
print
    "here is the freq  matrix for shifted tags:\n";
#&MathsTools::get_get_total_number_of_seqs_used_in_freq_matrix($freq_matrix_shifted_tags);
#&PrintingTools::print_out_a_matrix($pwm_shifted_tags );
&PrintingTools::print_out_a_matrix_into_a_file($pwm_shifted_tags, $pwm_shifte_tags_file, ",");

# my $number_of_rows  = @{$freq_matrix_real_tags};
# my $number_of_cols  = @{$freq_matrix_real_tags[0]};
# print
#     "rows are : $number_of_rows\n and number of collumns is $number_of_cols\n";

exit;



#####################################  SUBROUTINES ###################


sub get_position_frequency_matrix($$$$$){
    my $input_bed_file  = shift;
    my $start_offset    = shift;
    my $end_offset      = shift;
    my $bg_shift        = shift; #only for back ground calculations otherwise must be zero
    my $chr_seq         = shift;
    my $tag_length      = shift;


    
    my %unique_tags = &FileTools::GET_UNIQUE_FEATURES_FROM_A_BED_FILE_V2($input_bed_file);
    my $number_seqs_with_ns = 0;
    my $length = $tag_length + $start_offset + $end_offset;
   #initialize composzition matrix and assign zero to each of its cells
    my $com_matrix = [];
    for (my $r = 0; $r < $length; $r++) {
	for (my $c = 0; $c < 4; $c++) {
	    $com_matrix->[$r][$c] = 0;
	}
    }




    foreach my $key(keys %unique_tags){
	my $one_feat = $unique_tags{$key};
	my $chr    = $one_feat->{Chr};
	my $start  = $one_feat->{Start};
	my $end  = $one_feat->{End};
	my $source_tag = $one_feat->{Source};
	my $strand = $one_feat->{Strand};
	
	
	
	if($chr =~ m/^chr/i){
	    $chr =~ s/chr//;
	}
	my $new_start;
	my $new_end;
	
	if($strand eq "+"){
	    $new_start = $start - $start_offset + $bg_shift;
	    $new_end   = $end + $end_offset  + $bg_shift;
	}
	elsif($strand eq "-"){
	    $new_start = $start - $end_offset - $bg_shift;
	    $new_end   = $end   + $start_offset- $bg_shift;
	}
	else{
	    die "Unknown strand found at $chr\t $start\t $end \t $source_tag\t $strand\n";
	}
	my $seq_stretch = uc(substr($chr_seq,$new_start,$length));
	if($seq_stretch =~ m/N/){
	    $number_seqs_with_ns++;
	    next;
	}
	my $seq_length = length($seq_stretch);
	if($seq_length ne $length){
	    die "length are not matching!\n";
	}
	if ($strand eq "-") {
	    $seq_stretch = &SeqTools::get_revcomp($seq_stretch);
	}

	for (my $p = 0; $p < $length; $p++) {
	    my $one_base = substr($seq_stretch,$p,1);
	    if ($one_base eq "A") {
		$com_matrix->[$p][0] = $com_matrix->[$p][0] +1;
	    }
	    elsif ($one_base eq "C") {
		$com_matrix->[$p][1] = $com_matrix->[$p][1] +1;
	    }
	    elsif ($one_base eq "G") {
		$com_matrix->[$p][2] = $com_matrix->[$p][2] +1;
	    }
	    elsif ($one_base eq "T") {
		$com_matrix->[$p][3] = $com_matrix->[$p][3] +1;
	    }
	    else{
		die "Expected nuceotides but got $one_base!\n";
	    }
	}#p#
	
    }#foreach#
    print 
	"Found $number_seqs_with_ns of seqs containing at least one N!\n";
    return $com_matrix;
}#get_position_frequency_matrix#




