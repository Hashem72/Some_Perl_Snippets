#!/usr/bin/perl
use strict;
use warnings;
use lib '/nfs/users/nfs_h/hk3/src/bioperl-live';
use lib '/nfs/users/nfs_h/hk3/src/ensembl/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-compara/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-variation/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-functgenomics/modules';
use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use Bio::AlignIO;



# make connection to Ensembl Registry
my $registry = 'Bio::EnsEMBL::Registry';
#my $Chr17_seq = &getChromosomeSequence('Human','true','17');
#&pickRandomSubregionsofChromosome($Chr17_seq, 10);
&pickSomeRandomSubregionsofChromosome('Human',50, 'false',10,'22');
exit;


sub getChromosomeSequence($$;$){
	my $Species       = shift;
	my $RepeatMasked  = shift;
	my ($Chromosome)        = (@_ > 0) ? shift @_: '1'; #chromosome '1' = default if no chromosome specified
	
	#load registry
	$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
	);
	
	#get slice adpator:
	my $slice_adaptor          = $registry->get_adaptor($Species,'Core', 'Slice');
	my $slice                  = $slice_adaptor->fetch_by_region('chromosome', $Chromosome);
	my $slicemasked            = $slice->get_repeatmasked_seq();
	#changes type to Bio:EnsEMBL::RepeatMaskedSlice	
	my $sequence;
	if($RepeatMasked eq 'true'){
		print STDERR
			"RepeatMasking, please be patient!\n";
		$sequence = uc($slicemasked->seq())
	}
	else{
		$sequence = uc($slice->seq());
	}
		return $sequence;
}#getChromosomeSequence#

sub pickRandomSubregionsofChromosome($$){
	my $Chro_Sequence        = shift;
	my $Length_of_Sequence  = shift;

	my $Sequence;
	my $SuccessfullSelection = 'false';
	my $Length_of_Chromosome = length($Chro_Sequence);
	while($SuccessfullSelection ne 'true'){
		my $Allowed_Start_Position = $Length_of_Chromosome-$Length_of_Sequence-1;
		my $Start_Position         = int(rand($Allowed_Start_Position));
			$Sequence              = substr($Chro_Sequence,$Start_Position,$Length_of_Sequence);
			$SuccessfullSelection  = 'true';
			my $i                  = 0;
			while(  ($SuccessfullSelection eq 'true') && ($i< $Length_of_Sequence)){
				if(substr($Sequence,$i,1) eq 'N'){
					$SuccessfullSelection = 'false'
				}
				$i++;
			}
	}
	return $Sequence;
}#pickRandomSubregionsofChromosome#

sub pickSomeRandomSubregionsofChromosome($$$$;$){
	my $Species              = shift;
	my $Number_of_Sequences  = shift;
	my $RepeatMasked         = shift;
	my $Length_of_Sequences  = shift;	
	my ($Chromosome)        = (@_ > 0) ? shift @_: '1'; #chromosome '1' = default if no chromosome specified
	
	
	my @Vector_of_Sequences_and_Ids ;
	my $Chro_Sequence         = &getChromosomeSequence($Species,$RepeatMasked,$Chromosome );
	for(my $i=0; $i< $Number_of_Sequences; $i++){
		print STDERR 
			"getting random sequence number: ", $i, "\n";
		my $one_seq = &pickRandomSubregionsofChromosome($Chro_Sequence,$Length_of_Sequences);
		my $one_id  = "seq_".$i;
		my $one_pair = {
			SEQ     => $one_seq,
			SEQ_ID  =>$one_id 
		};
		push(@Vector_of_Sequences_and_Ids,$one_pair)
	}
	my $OutPutFile = '/nfs/users/nfs_h/hk3/RandomlyPickedSequences_from_chr_22.fasta';
	my $seqout = Bio::SeqIO->new(-file => ">$OutPutFile", -format => 'Fasta');
		foreach my $pair(@Vector_of_Sequences_and_Ids){
			my $id = $pair->{SEQ_ID};
			my $seq = $pair->{SEQ};
			my 	$bioseq = Bio::Seq->new(
			-seq => $seq,
			-id  => $id
			);
			$seqout ->write_seq($bioseq);
		}	
}#pickSomeRandomSubregionsofChromosome#
