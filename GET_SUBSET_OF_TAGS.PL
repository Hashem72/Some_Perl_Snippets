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
use IO::File;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

# make connection to Ensembl Registry
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);




##################################################################################################
# This is to get a subset of tags for motif analysis. It can be either a random or               #
# outliers either from positive (top socred B tags) or negative side(top socred A tag)           #
#Note: for B tags use 'ascending' for A tags use 'descending' else for randomly picked tags use  #
#random                                                                                          #
##################################################################################################  

&GET_SEQUENCES("Human","/nfs/th_group/hk3/UW_DNaseI_HS/H1hesc/file_A.bed",500, 'random',36,"/nfs/th_group/hk3/UW_DNaseI_HS/H1hesc/Motif_Discovery/H1hesc_file_A_random_3.fasta", 'true');
#my %bed_feats_file_B =&GET_SUBSET_OF_TAGS("/nfs/th_group/hk3/UW_DNaseI_HS/H1hesc/file_A.bed",50,'descending');



exit;

sub GET_SEQUENCES($$$$$$$){
    my $species                           = shift;
    my $input_bed_file                    = shift;
    my $required_number_of_tags           = shift;
    my $ascending_or_descenging_or_random = shift;
    my $offset                            = shift;
    my $output_fasta_file                 = shift;
    my $repeat_mask                       = shift;
    
    my $seqout = Bio::SeqIO->new('-format'=> 'Fasta', -file=> ">$output_fasta_file");
    my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');
    my %subset_of_tags = &GET_SUBSET_OF_TAGS($input_bed_file, $required_number_of_tags, $ascending_or_descenging_or_random);
    my @scores = keys %subset_of_tags;
    my $number_of_ignored_seq = 0;
    foreach my $score(@scores){
	my $one_bed_feature = $subset_of_tags{$score};
	my $chr = $one_bed_feature->{Chr};
	$chr =~ s/chr//;
	my $start = $one_bed_feature->{Start};
	my $end   = $one_bed_feature->{End};
	my $source = $one_bed_feature->{Source};
	my $score = $one_bed_feature->{Score};
	my $strand = $one_bed_feature->{Strand};
	my $seq;
	my $seq_start;
	my $seq_end;
	my $seq_id = $species."_".$chr."_".$source."_".$strand;
	if($strand eq '+'){
	    $seq_start = $start-$offset;
	     $seq_end    = $start+$offset;
	     
	     my $seq_slice = $slice_adaptor->fetch_by_region('chromosome',$chr,$seq_start,$seq_end);
	    $seq = uc($seq_slice->seq());
	     if($repeat_mask eq  'true'){
		my $masked_seq_slice = $seq_slice->get_repeatmasked_seq();
		$seq = uc($masked_seq_slice->seq());
	    }
	}
	elsif($strand eq '-'){
	    $seq_start = $end-$offset;
	    $seq_end   = $end+$offset;
	    my $seq_slice = $slice_adaptor->fetch_by_region('chromosome',$chr,$seq_start,$seq_end);
	    $seq = uc($seq_slice->seq());
	    if($repeat_mask eq 'true'){
		my $masked_seq_slice = $seq_slice->get_repeatmasked_seq();
		$seq = uc($masked_seq_slice->seq());
		$seq = &GET_REVCOMP($seq);
	    }
	}
	else{
	    print
		"Unknown strand!: at moment only + and - are accepted!\n";
	}
#ignore the seq if contains N othewise print it into the fasta file:
	if($seq =~ m/N/){
	    $number_of_ignored_seq++;
	}
	else{
	    my $seq_object = Bio::Seq->new(-seq => $seq, -id => $seq_id);
	    $seqout->write_seq($seq_object);
	}
	
    }
    print
	$number_of_ignored_seq. " sequences ignored because they were containing some Ns!\n";
    $seqout->close();
		  }#GET_SEQUENCES#


sub GET_SUBSET_OF_TAGS($$$){
    my $input_bed_file                    = shift;
    my $number_required_top_scores        = shift;
    my $ascending_or_descenging_or_random = shift;

    my %subset_of_tags;
   

    my %bed_features = &READ_A_TAG_FILE_INTO_A_HASH($input_bed_file);
    if($ascending_or_descenging_or_random eq 'random'){
	my @scores         = keys %bed_features;
	my $number_of_keys = @scores;
	my @random_indices = &GET_SOME_UNIQUE_RANDOM_INTEGERS($number_of_keys,$number_required_top_scores);
	foreach my $r( @random_indices){
	    my $one_key = $scores[$r];
	    my $one_bed_obj = $bed_features{$one_key};
	    $subset_of_tags{$one_key} = $one_bed_obj; 
	    
	}
	return %subset_of_tags;
    }
    else{
	my @sorted_scores;
	if($ascending_or_descenging_or_random eq "ascending"){
	    @sorted_scores = sort { $a <=> $b } keys %bed_features; 
	}
	elsif($ascending_or_descenging_or_random eq "descending"){
	    @sorted_scores = sort { $b <=> $a } keys %bed_features;
	}
	else{
	    print
		"Found some thing which is not ascending or descending or random!\n";
	}

	my $number_of_unique_scores = @sorted_scores;
	if($number_required_top_scores > $number_of_unique_scores){
	    print
		"\nNote: Number of top scoredd motifs you asked is greater than number of uniquely scored tags so all of them will be considered\n";
	    $number_required_top_scores = $number_of_unique_scores
	}
	my $upper_bound = $number_of_unique_scores-$number_required_top_scores;
	for my $index($upper_bound..$number_of_unique_scores-1){
	    my $score = $sorted_scores[$index];
	    print 
		"index = ". $index . " score = ". $score."\n"; 
	    my $one_bed_object = $bed_features{$score};
	    $subset_of_tags{$score} = $one_bed_object;
	}
	return %subset_of_tags;
	
    }
		       }#GET_SUBSET_OF_TAGS#



sub GET_SOME_UNIQUE_RANDOM_INTEGERS($$){
    my $range            = shift;
    my $requred_number   = shift;
    
    my @unique_random_integers;
    my $number_generated_unique_integers = 0;
    while( $number_generated_unique_integers < $requred_number){
	my $one_random = int(rand($range));
	push(@unique_random_integers,$one_random);
	my %seen;
	@unique_random_integers = grep { ! $seen{$_}++ } @unique_random_integers;
	$number_generated_unique_integers = @unique_random_integers;
    }
    my $number = @unique_random_integers;
    return @unique_random_integers;

}#GET_SOME_UNIQUE_RANDOM_INTEGERS#


sub READ_A_TAG_FILE_INTO_A_HASH($){
    #read a taf file into a hash: scores of tags are keys and lines are values
    my $input_file_name = shift;

    my %result;
    open(IN,$input_file_name) or die "cannot open file $input_file_name to read the data\n";
    while(<IN>){
	chomp;
	my $line = $_;
	my @line_spec = split("\t", $line);
	my $chr    = $line_spec[0];
	my $start  = $line_spec[1];
	my $end    = $line_spec[2];
	my $source = $line_spec[3];
	my $score  = $line_spec[4];
	my $strand = $line_spec[5];
	my $one_bed_line_object = {
	    Chr    => $chr,
	    Start  => $start,
	    End    =>$end,
	    Source => $source,
	    Score  => $score,
	    Strand => $strand
	};
	$result{$score}= $one_bed_line_object;
    }
    close(IN);
    return %result;
}#READ_A_TAG_FILE_INTO_A_HASH#


sub GET_REVCOMP($){
    my $dna_seq  = shift;
    
    my $revcom_seq = reverse $dna_seq;
    $revcom_seq =~ tr/NACGTnacgt/NTGCAntgca/;
    return $revcom_seq;
		}#GET_REVCOMP#

