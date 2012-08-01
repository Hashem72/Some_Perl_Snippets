#!/usr/bin/perl
use strict;
use warnings;

use lib "/nfs/users/nfs_h/hk3/src/lib/perl5";
use lib '/nfs/users/nfs_h/hk3/src/bioperl-live';
use lib '/nfs/users/nfs_h/hk3/src/ensembl/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-compara/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-variation/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-functgenomics/modules';
use Bio::EnsEMBL::Registry;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Perl;
use Bio::FeatureIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;

######################### Note: ############################################################################
#  this is to convert wig files (ouput form hotspot algorithms to gff files. wig files are bed like files) #
#  they include chr start and end                                                                          #
#                                                                                                          #
############################################################################################################


use IO::File;
use Getopt::Long;

use Bio::EnsEMBL::Registry;



my ($input_file, $output_file, $length_threshold, $data_source, $data_type , $Chr);
my $help= '';
if(!GetOptions('input_file|i=s'       =>\$input_file,
	       'output_file|o=s'      => \$output_file,
	       'length_threshold|l=s' => \$length_threshold,
	       'data_source|s=s'      => \$data_source,
	       'data_type|t=s'        => \$data_type,
               'chr|c=s'              => \$Chr,
               'help|h'               => \$help)
 || !(defined($input_file) && defined($output_file) && defined($length_threshold) && defined($data_type) && defined($Chr)   ) || $help )
{
    print <<END_USAGE;
Usage:
$0  --input_file=input_file --output_file=output_file --length_threshold=length_threshold --data_source=data_source --data_type=data_type --Chr=chr
$0  --help
--input_file       /-i  Name of input file which is in Wig (or bed three collumns:chr start and end).
--output_file      /-o  Name of output file which will be as a gff(chr data source start end startn data type ).
--length_threshold /-l  Sequences shorter than this in wig file will not be considered.Put zero if you want each sequence to be considered.
--data_source      /-s  Data source, call it what ever you like.
--data_type        /-t  Data type: whatever you like.
--chr             /-c  chromosome name: example: chr22
--help             /-h  To see this text.
Example usage:
$0 -i inputfile.wig -o outputfile.gff -l 200 -s a_text_here -t anohter_text_here -c chr22
END_USAGE
exit(1);
}#end of if#










	#DHS in ENCODE duke and UW hotspot (http://encodewiki.ucsc.edu/EncodeDCC/index.php/Locations_of_ENCODE_Data#Post_Jan_2011_Freeze_data) is in wig files. this is to convert a wig to a giff
	
	my $wig_input_file      = $input_file;
	my $gff_output_file      = $output_file;
	my $length              = $length_threshold;

	
	
	open(WIG, "$wig_input_file") or die "Cannot open file $wig_input_file to read the data!\n";
	my @lines = <WIG>;
	close(WIG);
	
	my $gff_out    = new Bio::Tools::GFF(-file => ">$gff_output_file", -gff_version =>1);
	my $counter = 0;
	foreach my $line(@lines){
		chomp($line);
		next if ($line =~ m/track/);
		next if ($line !~ m/$Chr/);
		my @line_spec = split("\t", $line);
		my $chr = $line_spec[0];
		if($Chr ne $chr){
			print
				"Chr = ".  $Chr . "   chr = ". $chr . "\n";
			#die "check your chromosome name because what you have given me is not matching with what I have found in the file!\n";
		}
		my $start = $line_spec[1];
		my $end  = $line_spec[2];
		my $score = $line_spec[3];
	
	
		my $one_length = $end- $start;
		if($one_length >=   $length){
			$counter++;
			my $feature = new Bio::SeqFeature::Generic(-start => $start,
 			-end         => $end,
 			-seq_id      => $chr,
 			-source_tag  => $data_source."_".$counter,
 			-primary     => $data_type,
 			-score       => $score
 		);
 
			$gff_out->write_feature($feature);
			
		}
		
	}
	print
			
		"There were ". $counter . " features longer than or equal to ". $length .  "\n";	
