#!/usr/bin/perl
use strict;
use warnings;

use lib "/nfs/users/nfs_h/hk3/My_Perl_Scripts/modules";
use Ensembl;
use FileTools;
use SeqTools;
use MathsTools;
use PrintingTools;
use K_Mers;

my $k = 8;
my$repeat_masked = "true";
#&Ensembl::get_chromosme_sequence_and_write_into_a_fasta_file("Human", "22", "true", "/nfs/th_group/hk3/Human_CHRS/chr_22_RepeatMasked.fasta");

my $chr_file = '';
my $rp_extention = '';
if($repeat_masked eq "true"){
    $chr_file = &FileTools::get_seq_from_fasta_file("/nfs/th_group/hk3/Human_CHRS/chr_22_RepeatMasked.fasta","Human_chromosome_22");
    $rp_extention = "RepeatMasked";
}
else{
    $chr_file = &FileTools::get_seq_from_fasta_file("/nfs/th_group/hk3/Human_CHRS/chr_22.fasta", "Human_chromosome_22");
    $rp_extention = "UnRepeatMasked"
}
my @shifting_lengths_for_real_tags            = (0);
my @rand_integers = &MathsTools::get_an_array_of_random_integers(10,500,50);


my @shifting_lengths_for_background_estimation = @rand_integers; #(10,20,30,40,50,60,70,80,90,100);
my $shifting_lengths_for_background_estimation_ref = \@shifting_lengths_for_background_estimation;

my $tags_bed_file = "/nfs/th_group/hk3/Duke_DNaseI_HS/K562/wgEncodeOpenChromDnaseK562AlnRep1_chr22.bed";
my $k_mers_counts_file_real = "/nfs/th_group/hk3/Duke_DNaseI_HS/K562_For_Paper_Analysis/".$rp_extention."_".$k."_mers_counts_real.txt";


my $k_mers_counts_file_bg = "/nfs/th_group/hk3/Duke_DNaseI_HS/K562_For_Paper_Analysis/".$rp_extention."_".$k."_mers_counts_bg.txt";

my @k_mers = &K_Mers::make_arrays_of_k_mers($k);
my $k_mers_ref = \@k_mers;

my %hash_of_k_mers = &K_Mers::array_to_hash_with_zero_values($k_mers_ref);
my $hash_of_k_mers_ref = \%hash_of_k_mers;

&K_Mers::count_number_of_occurrances_of_k_mers_in_tags($tags_bed_file, $k_mers_counts_file_real, 4, -15, \@shifting_lengths_for_real_tags, $hash_of_k_mers_ref, $chr_file);

&K_Mers::count_number_of_occurrances_of_k_mers_in_tags($tags_bed_file, $k_mers_counts_file_bg, 4, -15, $shifting_lengths_for_background_estimation_ref, $hash_of_k_mers_ref, $chr_file);

exit;


########################################################





   
