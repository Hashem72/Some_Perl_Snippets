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
#  given a gff file including information about some sequences, this it to contact Ensembl and fetch the   #
#  sequences with reuired flanking lengths.                                                                #
############################################################################################################


use IO::File;
use Getopt::Long;

use Bio::EnsEMBL::Registry;




# make connection to Ensembl Registry
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);





my ($species, $input,$output,$flanking_length,$repeatmasking);
my $help='';

if( !GetOptions('species|s=s'         => \$species,
                'input|i=s'           => \$input,
                'output|o=s'          => \$output,
                'flanking_length|f=s' => \$flanking_length,
                'repeatmasking|r=s'   => \$repeatmasking,
                'help|h'              => \$help) 
  || !(defined($species) && defined($input) && defined($output) && defined($flanking_length) && defined($repeatmasking)  ) || $help )
{
    print <<END_USAGE;
Usage:
$0 --species=species --input=input --output=output --flanking_length=flanking_length --repeatmasking=repeatmasking
$0 --help
--species  /-s name of species must be known to Ensembl for example Human
--input   /-i input gff file name
--ouput  /-o ouput fasta file name
--flanking_length  \-f flanking length you wish to get sequences from both sides
--repeatmasking   \-r whether to repeatmask or not true or false is accepted.
--help            \-h to print this message.
Example usage:
$0 -s Human -i input_file.gff -o output_file.fasta -l 100 -r false
END_USAGE
exit(1);
}



    


	# to read a feature from a gff file and then fetch their seqs from Ensemb
	
my $flank              = $flanking_length;
my $gff_input_file     = $input;
my $fasta_output_file  = $output;
my $repeatmask        = $repeatmasking;


my $gff_in = Bio::Tools::GFF->new(-file => $gff_input_file,  -gff_version => 1);
my $seqout  = Bio::SeqIO->new('-format'=> 'Fasta', -file=> ">$fasta_output_file");
#get the sliece adaptor:
my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');

my $feature;
while($feature = $gff_in->next_feature()){
		my $start = $feature->start();
		my $end   = $feature->end();
		my $chr   = $feature->seq_id();
		my $primery    = $feature->primary_tag();
		$chr     =~ s/chr//;
		my $source_tag = $feature->source_tag();
		my $new_start = $start-$flank;
		my $new_end    = $end+$flank;
		my $seq_id    = $source_tag."_".$primery;  #$chr."_".$source_tag."_".$new_start."_".$new_end;
		my $seq_slice = $slice_adaptor->fetch_by_region('chromosome',$chr,$new_start,$new_end);
		my $seq = uc($seq_slice->seq());
		my $seq_length = length($seq);
		if($repeatmask eq 'true'){
		    my $masked_seq_slice = $seq_slice->get_repeatmasked_seq();
		    $seq   = uc($masked_seq_slice->seq());
		    $seq_length = length($seq);
		} 
		if($seq =~ m/N/){
		    print
			" a sequence ignored because it was including Ns\n";	
		}
		else{
		    my $seq_obj  = Bio::Seq->new(-seq => $seq, -id =>$seq_id);
			$seqout->write_seq($seq_obj);
		    print
		    "fetched seq $seq_id with length $seq_length.\n";

		}
}
$gff_in->close();
$seqout->close();



exit;
