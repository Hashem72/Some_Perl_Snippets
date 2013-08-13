#!/usr/bin/perl
use strict;
use warnings;

use lib "/nfs/users/nfs_h/hk3/My_Perl_Scripts/modules";
use Ensembl;
use FileTools;
use SeqTools;
use MathsTools;
use PrintingTools;


my $tag_resource          = "Duke";
my $cell_line             = "Gm12878";
my $Rep                   = "rep1";
my $five_prime_offset     = 10;  
my $three_prime_offset    = 10;
my $tag_length;

my $tags_bed_file;
my $dinucleo_output_file;

if($tag_resource eq "UW"){
    $tag_length            = 36;

#for randomly Generated tags use the following two lines:
    #$tags_bed_file          = "/nfs/th_group/hk3/UW_DNaseI_HS/Gm12878/Randomly_Generated_Tags.bed";
    #$dinucleo_output_file  = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."_For_Paper_Analysis/dinucleotides_frequency_from_Gm12878_randomly_generated_tags.txt";
    $tags_bed_file         = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."/wgEncodeUwDnase".$cell_line."AlnRep1_chr22.bed";
    $dinucleo_output_file  = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."_For_Paper_Analysis/wgEncodeUwDnase".$cell_line."Aln".$Rep."_chr22_dinucleotides_frequency.txt";   
}
elsif($tag_resource eq "Duke"){
    $tag_length             = 20;
    #for randomly Generated tags use the following two lines:
    $tags_bed_file          =  "/nfs/th_group/hk3/Duke_DNaseI_HS/Gm12878/Randomly_Generated_Tags.bed";
    $dinucleo_output_file   =  "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."_For_Paper_Analysis/dinucleotides_frequency_from_Gm12878_randomly_generated_tags.txt";
   # $tags_bed_file          = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."/wgEncodeOpenChromDnase".$cell_line."AlnRep1_chr22.bed";
    #$dinucleo_output_file   = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."_For_Paper_Analysis/wgEncodeOpenChromDnase".$cell_line."Aln".$Rep."_chr22_dinucleotides_frequency.txt";
}
else{
    die "Unknown tag resource\n";
}


my $chr_22 = &FileTools::get_seq_from_fasta_file("/nfs/th_group/hk3/Human_CHRS/chr_22.fasta", "Human_chromosome_22");
print
    "chromosome sequence extracted!\n";



     
########################

my @dinucleotides = ("AA", "AC", "AG", "AT",
		     "CA", "CC", "CG", "CT", 
		     "GA", "GC", "GG", "GT",
		     "TA", "TC", "TG", "TT");
my $stetch_length = $tag_length + $five_prime_offset + $three_prime_offset;

#initialize an array of hashes
my @dinucleos ;
for (my $p = 0; $p<$stetch_length-1; $p++){
    foreach my $d(@dinucleotides) {
	$dinucleos[$p]{$d}=0;
    }

}


#get_unique tags:
my %unique_features = &FileTools::GET_UNIQUE_FEATURES_FROM_A_BED_FILE_V2($tags_bed_file );
my $number_of_unique_tags = keys %unique_features;
print
    "Found ". $number_of_unique_tags. " unique tags!\n";
foreach my $key(keys %unique_features){
	my $one_feat  = $unique_features{$key};
	my $chro      = $one_feat->{Chr};
	my $start     = $one_feat->{Start};
	my $end       = $one_feat->{End};
	my $source    = $one_feat->{Source};
	my $strand    = $one_feat->{Strand};
	my $score     = $one_feat->{Score};
	if ($strand ne "+" && $strand ne "-" ){
	print
	    ">>>>".$strand ."<<<<\n";
	print $chro."\t".$start."\t".$end."\t".$source."\t".$score."\t".$strand."\n";
	
	}

   
    #apply offsets
    my $modified_start = $start - $five_prime_offset;
    my $modified_end   = $end   + $three_prime_offset;
    my $subsequence_length = $end + $three_prime_offset - ($start- $five_prime_offset);
    
    # extract this subsequence from  the whole chromosome
	

    my $seq_stretch = uc(substr($chr_22, $modified_start, $subsequence_length));
	next if $seq_stretch =~ m/N/;
	
if ($strand eq "-") {
    $seq_stretch = uc(&SeqTools::get_revcomp($seq_stretch));
}

    for (my $p = 0; $p < $subsequence_length-1; $p++) {
	
	my $one_dinucleotide = substr($seq_stretch, $p, 2);
	$dinucleos[$p]{$one_dinucleotide} = $dinucleos[$p]{$one_dinucleotide} + 1;
    }
}


open DINCLEOTIDES, ">$dinucleo_output_file" or die "cannot open file $dinucleo_output_file to write into!\n";

for (my $x = 0; $x<$stetch_length-1 ; $x++){
    foreach my $d(@dinucleotides) {
	print  DINCLEOTIDES
	    $dinucleos[$x]{$d}."\t";
    }
    print DINCLEOTIDES
	"\n";
}

close(DINCLEOTIDES);
