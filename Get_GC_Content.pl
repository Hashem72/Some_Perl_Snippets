#!/usr/bin/perl
use strict;
use warnings;

use lib "/nfs/users/nfs_h/hk3/My_Perl_Scripts/modules";
use Ensembl;
use FileTools;
use SeqTools;
use MathsTools;
use PrintingTools;







my $data_centre       = "Duke";
my $cell_line         = "Gm12878";
my $replicate         = 1;
my $five_prime_offset;
my $three_prime_offset;
my $bed_file;
my $output_file;


#############################


if ($data_centre eq "UW") {
    $five_prime_offset     = 10;
    $three_prime_offset    = 10;
     $bed_file              = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."/wgEncodeUwDnase".$cell_line."AlnRep".$replicate."_chr22.bed";
     $output_file    = "/nfs/th_group/hk3/UW_DNaseI_HS/".$cell_line."_For_Paper_Analysis/GC_Content_rep".$replicate.".txt";

}
elsif($data_centre eq "Duke"){
    $five_prime_offset      = 10;
    $three_prime_offset     = 10;
    $bed_file               = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."/wgEncodeOpenChromDnase".$cell_line."AlnRep".$replicate."_chr22.bed";
    $output_file     = "/nfs/th_group/hk3/Duke_DNaseI_HS/".$cell_line."_For_Paper_Analysis/GC_Content".$replicate.".txt";
 }
else{
    die "Unknown data centre!\n";
}




##########################








my $Human_chr_22             = &Ensembl::GET_CHROMOSOME_SEQUENCE("Human","22","false");
print "fetched the chromosome sequence!\n";
&get_gc_content_v2($bed_file, $output_file, $Human_chr_22, $five_prime_offset ,  $five_prime_offset );

exit;


################## subroutines ###################################


sub get_gc_content_v2($$$$$$){
    my $inp_bed_file         = shift;
    my $otp_bed_file         = shift;
    my $chromosome_sequence  = shift;
    my $five_prime_offset    = shift;
    my $three_prime_offset   = shift;
    
    #my $off_set_length       = shift;
   # my $shifting_length      = shift; # if you want to move reads around and get an estimaiton of background
   # my $stretch_length       = shift; # length of stretch you want to get its gc content, the max must be equal to tag length 
    open(OUTPUT, ">$otp_bed_file") or die "cannot open file $otp_bed_file to write in!\n";
    print OUTPUT
	"chr\tstart\tscore\tgc_content\tstrand\n";

    #get unique features
    my %unique_tags = &FileTools::GET_UNIQUE_FEATURES_FROM_A_BED_FILE_V2($inp_bed_file);
    foreach my $key(keys %unique_tags){
	my $one_feat        = $unique_tags{$key};
	my $chr             = $one_feat->{Chr};
	my $start           = $one_feat->{Start};
	my $end             = $one_feat->{End};
	my $source          = $one_feat->{Source};
	my $score           = $one_feat->{Score};
	my ($rounded_score) = $score =~ /(\S{0,6})/;
	my $strand          = $one_feat->{Strand};

	my $new_start       = $start - $five_prime_offset;
	my $new_end         = $end   + $three_prime_offset;
        my $stretch_length  = $new_end - $new_start ;

	my $seq_stretch = uc(substr($chromosome_sequence, $new_start, $stretch_length ));

	if ($strand eq "-"){
	    $seq_stretch = &SeqTools::get_revcomp($seq_stretch);
	}
       #if sequences contained N, ignore it
	next if ($seq_stretch =~ m/N/);
	
	my $gc_content =  &SeqTools::get_gc_content($seq_stretch);
	#print
	 #   $gc_content."\n";
	print OUTPUT
	    $gc_content. "\n";
	
    }
    close(OUTPUT);
}#get_gc_content_v2#










sub get_gc_content($$$$$$){
    my $inp_bed_file         = shift;
    my $otp_bed_file         = shift;
    my $chromosome_sequence  = shift;
    my $off_set_length       = shift;
    my $shifting_length      = shift; # if you want to move reads around and get an estimaiton of background
    my $stretch_length       = shift; # length of stretch you want to get its gc content, the max must be equal to tag length 
    open(OUTPUT, ">$otp_bed_file") or die "cannot open file $otp_bed_file to write in!\n";
    print OUTPUT
	"chr\tstart\tscore\tgc_content\tstrand\n";

    #get unique features
    my %unique_tags = &FileTools::GET_UNIQUE_FEATURES_FROM_A_BED_FILE_V2($inp_bed_file);
    foreach my $key(keys %unique_tags){
	my $one_feat = $unique_tags{$key};
	my $chr = $one_feat->{Chr};
	my $start = $one_feat->{Start};
	my $end   = $one_feat->{End};
	my $source  = $one_feat->{Source};
	my $score   = $one_feat->{Score};
	my ($rounded_score) = $score =~ /(\S{0,6})/;
	my $strand  = $one_feat->{Strand};

	my $new_start;
	my $new_end;
	my $seq_stretch;

	if($strand eq "+"){
	    $new_start = $start - $off_set_length + $shifting_length;
	    $seq_stretch = uc(substr($chromosome_sequence, $new_start, $stretch_length ));
	}
	elsif ($strand eq "-"){
	    $new_end     = $end + $off_set_length - $shifting_length;
	    $new_start   = $new_end - $stretch_length;
	    $seq_stretch = uc(substr($chromosome_sequence, $new_start, $stretch_length ));
	    $seq_stretch = &SeqTools::get_revcomp($seq_stretch);
	}
	else{
	    print "Unknown strand: $chr\t$start\tend\tscore\tstrand\n";
	}

        #if the whole sequence was made of Ns, ingore it:
	my $number_of_known_nts  = 0;
	while( $seq_stretch =~ m/[ACGT]/gi){
	    $number_of_known_nts++;
	}
	next if $number_of_known_nts eq 0;
	
	

	
	my $gc_content =  &SeqTools::get_gc_content($seq_stretch);
	print OUTPUT
	    $chr."\t".$start."\t".$rounded_score."\t".$gc_content. "\t". $strand."\n";
	
    }
    close(OUTPUT);
}#get_gc_content#


