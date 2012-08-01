#Please Note: this snippet is very inefficient for getting large number of random sequences. It spends much of its run time communicating with the server.If you plan to get large number of random sequences, please use GET_RANDOM_GENOMIC_SUBSEQUENCE_2.pl.


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
&get_some_random_sequences_from_the_genome('Human',10,1000000,'false','17');
#&get_random_piece_of_genomic_sequence('Human',5,'false','17');
exit;


sub get_random_piece_of_genomic_sequence($$$;$){
	my $Species             = shift; 
	my $Sequence_Length     = shift;
	my $Repeat_Masked       = shift;
	my ($Chromosome)        = (@_ > 0) ? shift @_: '1'; #chromosome '1' = default if no chromosome specified
	
	#load registry
	$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
	);
	
	#get slice adpator:
	my $slice_adaptor          = $registry->get_adaptor($Species,'Core', 'Slice');
	my $slice                  = $slice_adaptor->fetch_by_region('chromosome', $Chromosome,1,100);
	my $Chromosome_length      = $slice->seq_region_length();
	my $successfullselection   = 'false';
	my $sequence;
	while($successfullselection ne 'true'){
		my $allowed_start_position = $Chromosome_length - $Sequence_Length -1 ;
		my $start_position  = int(rand($allowed_start_position));
			print
				"Randomly selected starting position is: ", $start_position, "\n";
		my $end_position    = $start_position+ $Sequence_Length-1;
		my $slice           = $slice_adaptor->fetch_by_region('chromosome', $Chromosome,$start_position,$end_position);
		my $slicemasked     = $slice->get_repeatmasked_seq();
			#changes type to Bio:EnsEMBL::RepeatMaskedSlice
			if( $Repeat_Masked eq 'true'){
				$sequence= uc($slicemasked->seq())
			}
			else{
				$sequence = uc($slice->seq());
			}
			$successfullselection = 'true';
			my $i= 0;
			while(($successfullselection eq 'true') && ($i < $Sequence_Length ) ){
				if(substr($sequence,$i,1) eq 'N'){
					$successfullselection = 'false';		
				}
				$i++;
			}
	}
	return $sequence;
	}#get_random_piece_of_genomic_sequence#
	
	sub get_some_random_sequences_from_the_genome($$$$;$){
		my $Species                = shift;
		my $Length_of_Sequences    = shift;
		my $Number_of_Sequences    = shift;
		my $Repeat_Masked          = shift;
		my ($Chromosome)        = (@_ > 0) ? shift @_: '1'; #chromosome '1' = default if no chromosome specified 
		
		my @Vector_of_piars_of_ids_and_seqs;
		for(my $i= 0; $i< $Number_of_Sequences; $i++){
			print
				"getting random sequence number: ", $i, "\n";
			my $one_seq = &get_random_piece_of_genomic_sequence($Species,$Length_of_Sequences,$Repeat_Masked,$Chromosome);
			my $one_id  = "Seq_".$i;
			my $one_piar = {
				SEQ    => $one_seq,
				SEQ_ID => $one_id
			};
			push(@Vector_of_piars_of_ids_and_seqs,$one_piar)
		}
		my $OutPutFile = '/lustre/scratch103/sanger/hk3/Five_Vertebrate/Randomly_Picked_Sequences/RandomlyPickedSequences_for_cebpa.fasta';
		my $seqout = Bio::SeqIO->new(-file => ">$OutPutFile", -format => 'Fasta');
		foreach my $pair(@Vector_of_piars_of_ids_and_seqs){
			my $id = $pair->{SEQ_ID};
			my $seq = $pair->{SEQ};
			my 	$bioseq = Bio::Seq->new(
			-seq => $seq,
			-id  => $id
			);
			$seqout ->write_seq($bioseq);
		}
	}#get_some_random_sequences_from_the_genome#

   
