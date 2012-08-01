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


&TESTING("/nfs/th_group/hk3/UW_DNaseI_HS/K562/testing.bed", "/nfs/th_group/hk3/UW_DNaseI_HS/K562/dummy_testin.bed");
#&REMOVE_DUPLICATE_FEATURES_FROM_A_BED_FILE("/nfs/th_group/hk3/UW_DNaseI_HS/K562/wgEncodeUwDnaseK562RawDataRep1_chr22_only.bed", "/nfs/th_group/hk3/UW_DNaseI_HS/K562/wgEncodeUwDnaseK562AlnRep1_chr22_Unique_features_2.bed");
exit;
###############  subroutines #################################################

sub TESTING($$){
    my $input     = shift;
    my $output    = shift;


    my $in  =  Bio::FeatureIO->new(-format => 'bed', -file => $input);
    my $out =  Bio::FeatureIO->new(-format => 'bed', -file => ">$output");
    my $feat;
    while($feat = $in->next_feature()){
	
	
	my $start = $feat->start();
	$feat->start($start-1);
	my $end   = $feat->end();
	$end   = $feat->end($end-1);
	my $source   = $feat->source_tag();
	my $score = $feat->score();
	my $Chr   = $feat->seq_id();
	my $starnd = $feat->strand();
	$feat->strand(-1);
print
    "starnd ==== ". $feat->strand()."====\t";
	print
	    $Chr ."\t". $start ."\t".$end ."\t". $source  ."\t". $score . "\t". $starnd . "\n";
	
	$out->write_feature($feat);
    }


$in->close();
$out->close();
}

sub REMOVE_DUPLICATE_FEATURES_FROM_A_BED_FILE(){
#note: this is based on starting position of a feature ie if starting position is duplicated, then the feature is removed.
    my $input_bed_file         = shift;
    my $output_bed_file        = shift;




 my $in = Bio::FeatureIO->new(-format => 'bed', -file => $input_bed_file);
    my $feat;

    my %unique_features;
    my $lineCounter =0;
    while(  $feat = $in->next_feature()){
	my $start = $feat->start();
	$unique_features{$start} = $feat;
    }
    $in->close();

    my $out = Bio::FeatureIO->new (-format=>'bed' , -file => ">$output_bed_file");
    foreach my $start_pos(keys %unique_features){
	my $one_feat = $unique_features{$start_pos};
	$out->write_feature($one_feat);
    }
    
  $out->close();  
    
}#REMOVE_DUPLICATE_FEATURES_FROM_A_BED_FILE#






