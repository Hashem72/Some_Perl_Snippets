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
# this is to get sequences of features from the Ensembl. The user can choose to substitute sequence        #
#  their reverse complement if the sequence is in negative starnd                                          #
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





my ($species,$chr,  $input,$output,$shift_length,$repeatmasking, $revCom);
my $help='';

if( !GetOptions('species|s=s'         => \$species,
                'input|i=s'           => \$input,
		'chr|c=s'             => \$chr,
                'output|o=s'          => \$output,
                'shift_length|l=s' => \$shift_length,
                'repeatmasking|r=s'   => \$repeatmasking,
		'revCom|v=s'          =>  \$revCom,
                'help|h'              => \$help) 
  || !(defined($species) &&  defined($chr)   && defined($input) && defined($output) && defined($shift_length) && defined($repeatmasking) && defined($revCom) ) || $help )
{
    print <<END_USAGE;
Usage:
$0 --species=species --chr=chr --input=input --output=output --revCom=revCom --shift_length=shift_length --repeatmasking=repeatmasking
$0 --help
--species  /-s name of species must be known to Ensembl for example Human
--chr     /-c chromosome of interest eg 22 and not chr22
--input   /-i input as bed file, IT should contain at least six columns: chr_name, start and end, name, score and starnd. 
--ouput  /-o ouput fasta file name
--revCom /-v if true, then it will get the  reverse complement of sequences with negative strand, otherwise just the sequence.
--shift_length  \-l (start , end) of a feature  in pve strand will be replaced with (start+shift_length, end+flanking_lenght) and in nve strand will be (start-shift_length, end-shift_length) in other words pve features will be shifted to right and negative features will be shifted to left. 
--repeatmasking   \-r whether to repeatmask or not true or false is accepted.
--help            \-h to print this message.
Example usage:
$0 -s Human -c chr22 -i input_file.bed -o output_file.fasta -l 10 -r false -v true
END_USAGE
exit(1);
}



#my $chr_sequence = &GET_CHROMOSOME_SEQUENCE("Human", "chr22", "false");#($species, $chr,$repeatmasking );
#print 
 #   "chr22 fetched\n";

&get_feature_seqs_and_write_them_into_a_file($input, $output, $shift_length,$species, $chr, $repeatmasking, $revCom);


sub get_feature_seqs_and_write_them_into_a_file($$$$$$$){
    my $input_file_name      = shift;
    my $out_put_file_name    = shift;
    my $shifting_length      = shift;
    my $organism             = shift;
    my $chromosome           = shift;
    my $repeat_mask          = shift;
    my $rev_comp             = shift;             
    
    
    
    my $chr_sequence = &GET_CHROMOSOME_SEQUENCE($organism, $chromosome, $repeat_mask  );
    
    my $seqout  = Bio::SeqIO->new('-format'=> 'Fasta', -file=> ">$out_put_file_name");
    
    open(INPUT, $input_file_name) or die "coudn't open file $input_file_name to read the features!";
    while(<INPUT>){
	my $line = $_;
	chomp($line);
	my @fields = split("\t", $line);
	my $chr    = $fields[0];
	my $start  = $fields[1];
	my $end    = $fields[2];
	my $strand = $fields[5];
	
	my $seq_id = $start."_".$end."_".$strand;
	if( $strand eq '+'){
	    $start = $start + $shifting_length;
	   $end   = $end   + $shifting_length;
	}
	elsif ($strand eq "-"){
	    $start = $start - $shifting_length;
	    $end   = $end  - $shifting_length;
	}
	else{
	    die "uknown strand: = $strand detected\n";
	}
	
	my $stretch_length = $end - $start;
	my $stretch_seq    = uc(substr($chr_sequence, $start, $stretch_length));
	if($rev_comp eq 'true' && $strand eq '-'){
	    $stretch_seq = &GET_REVCOMP($stretch_seq);
	}
	my $seq_object = Bio::Seq->new(-seq => $stretch_seq, -id => $seq_id);
	$seqout->write_seq($seq_object);
    }
    close(INPUT);
    $seqout->close();
    
						}#get_feature_seqs_and_write_them_into_a_file#


sub GET_CHROMOSOME_SEQUENCE($$$){
    my $spec      = shift;
    my $chromosome   = shift;
    my $repeat_mask  = shift;
    
    $chromosome =~ s/chr//;
#get slice adaptor
    my $sa      = $registry->get_adaptor($spec,'Core','slice');
    my $slice = $sa->fetch_by_region('chromosome', $chromosome);
    my $sequence;
    if($repeat_mask eq 'true'){
	my $slice_masked  = $slice->get_repeatmasked_seq();
	print STDERR
	    "Repeatmasking, please be patient!\n";
	$sequence = uc($slice_masked->seq());
    }
    else{
	$sequence = uc($slice->seq());
    }
    return $sequence;   
			    }# GET_CHROMOSOME_SEQUENCE#

sub GET_REVCOMP($){
    my $dna_seq  = shift;
    
    my $revcom_seq = reverse $dna_seq;
    $revcom_seq =~ tr/NACGTnacgt/NTGCAntgca/;
    return $revcom_seq;
		}#GET_REVCOMP#






