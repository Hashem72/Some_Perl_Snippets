package Ensembl;
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

use IO::File;
use Getopt::Long;






# make connection to Ensembl Registry
my $registry = 'Bio::EnsEMBL::Registry';


#load registry
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => '5306'
    );



######################################################### subroutines ####################

sub GET_CHROMOSOME_SEQUENCE($$$){
    my $species      = shift;
    my $chromosome   = shift;
    my $repeat_mask  = shift;
    
    
#get slice adaptor
    my $sa      = $registry->get_adaptor($species,'Core','slice');
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
sub get_chromosme_sequence_and_write_into_a_fasta_file($$$$){
    my $species          = shift;
    my $chro             = shift;
    my $repeat_mask      = shift;
    my $fasta_file_name  = shift;

    my $chr_seq = &GET_CHROMOSOME_SEQUENCE($species, $chro, $repeat_mask);
    my $seq_id  = $species."_chromosome_".$chro;
    my $seqout = Bio::SeqIO->new('-format'=> 'Fasta', -file => ">$fasta_file_name");
    my $seq_obj = Bio::Seq->new(-seq => $chr_seq, -id => $seq_id);
    $seqout->write_seq($seq_obj);
    my $L = length($chr_seq);
    print
	"a seq with length $L written into $fasta_file_name\n";
}#get_chromosme_sequence_and_write_into_a_fasta_file#

1;
