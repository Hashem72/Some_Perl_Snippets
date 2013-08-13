#!/usr/bin/perl
use strict;
use warnings;

use lib "/nfs/users/nfs_h/hk3/My_Perl_Scripts/modules";
use Ensembl;
use FileTools;
use SeqTools;
use MathsTools;
use PrintingTools;

my $data_source = "UW";
my $cell_line = "Gm12878";


my $chr_22 = &FileTools::get_seq_from_fasta_file("/nfs/th_group/hk3/Human_CHRS/chr_22.fasta", "Human_chromosome_22");
my $input_bed_file = "/nfs/th_group/hk3/".$data_source."_DNaseI_HS/".$cell_line."_For_Paper_Analysis/".$cell_line."_Encode_hotspots.bed";
my $output_txt_file = "/nfs/th_group/hk3/".$data_source."_DNaseI_HS/".$cell_line."_For_Paper_Analysis/".$cell_line."_Freq_of_ncts_in_ENCODE_DHS.txt";

&print_out_freq_of_ncts_in_DHS($input_bed_file, $output_txt_file, $chr_22);

print "Frequencies of ncts were printed into : $output_txt_file\n";

exit;


#################################### subroutines ###########


sub print_out_freq_of_ncts_in_DHS($$$){
    my $input_bed_file  = shift;
    my $output_file     = shift;
    my $chr_seq         = shift;
    
    open(IN, "$input_bed_file") or die "cannot open file $input_bed_file to read the data!\n";
    open(OUT, ">$output_file") or die "cannot open file $output_file to write the data!\n";
    print OUT
	"A\tC\tG\tT\tnon_ACGT\n";
    
    while (<IN>) {
	chomp $_;
	my @fields = split("\t", $_);
	my $chr    = $fields[0];
	my $start  = $fields[1];
	my $end    = $fields[2];
	my $source = $fields[3];
	my $score  = $fields[4];
	my $strand  = $fields[5];
	if ( !(  defined($chr) && defined($start) && defined($end) && defined($source) && defined($score) && defined($strand)   )) {
	    die "required six columns but it didnt get it at: \n $_";
	}
	
	my $feature_length  = $end - $start;
	my $one_sequence  = uc(substr($chr_seq, $start, $feature_length));
	my %freqs = &SeqTools::get_freq_ncts($one_sequence);
	print OUT
	    $freqs{A}."\t".$freqs{C}."\t".$freqs{G}."\t".$freqs{T}."\t".$freqs{non_nct}."\n";
	    }

    close(IN);
    close(OUT);
}#print_out_freq_of_ncts_in_DHS#
