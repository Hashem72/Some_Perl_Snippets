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
use Bio::EnsEMBL::Utils::Exception qw(throw);
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
    -user => 'anonymous'
    );

my $species = "Human";

# get slice adaptor:
my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');
throw("Registry configuration file has no data for connecting to <$species>") if (!$slice_adaptor);

#get conservation_score adaptor:
#my $conservation_score_adaptor = $registry->get_adaptor("Multi", "compara", "ConservationScore");
my $cs_adaptor = $registry->get_adaptor("Multi", 'compara', 'ConservationScore');


thow("No conservation score") if(!defined($cs_adaptor));

#get method_link_species_set adaptor
my $mlss_adaptor = $registry->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");

#get method_link_species_set object for gerp conservation scores for mammals
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");

throw("Unable to find method_link_species_set") if (!defined($mlss));




&GET_CONSERVATION_SCORES_OF_SEQS("/nfs/th_group/hk3/UW_DNaseI_HS/H1hesc/group_B_hotspots.gff", 25,100);

####################################   subroutines ##############################################################


sub GET_CONSERVATION_SCORES_OF_SEQS($$$){
    my $gff_input_file         = shift;
    my $flanking_length        = shift;
    my $core_length            = shift;
    
    #get rid of extension
    my @file_spec = split(/\./, $gff_input_file);
    my $ext_upper_case = uc($file_spec[1]);
    if(  $ext_upper_case ne "GFF" ) {
	die "expecting to get a gff file but the extension ($file_spec[1]) doesnt look lile a gff file!\n";
    }
    my $five_prime_side_output_file   = $file_spec[0]."_Five_Prime_Conservation_Scores.txt";
    my $thtee_prime_side_output_file  = $file_spec[0]."_Three_Prime_Conservation_Scores.txt";
    my $core_output_file              =  $file_spec[0]."_Core_Conservation_Scores.txt";
    
    open(FOUT,">$five_prime_side_output_file") or die "Cannot open file $five_prime_side_output_file for writing!\n";
    open(TOUT, ">$thtee_prime_side_output_file") or die "Cannot open file $thtee_prime_side_output_file for writing!\n";
    open (COUT, ">$core_output_file") or die "Cannot open file $core_output_file for writing!\n";
    
    my $gff_in = Bio::Tools::GFF->new(-file =>$gff_input_file, -version =>1 );
    my $feature;
    while($feature = $gff_in->next_feature()){
	#original start and end
	my $start = $feature->start();
	my $end   = $feature->end();
	my $chr   = $feature->seq_id();
	$chr      =~ s/chr//;
	print
	    $chr."\t".$start."\t".$end."\n";

	#five prime
	my $f_start = $start -  $flanking_length;
	my $f_end   = $start   +  $flanking_length;
	#get slice
	my $f_slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $f_start,$f_end);
	throw("No Slice can be created with coordinates $chr:$f_start-$f_end") if (!$f_slice);
	my $f_display_size = $f_slice->end - $f_slice->start + 1; 
	my $f_scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $f_slice, $f_display_size);
	print
	    "FivePrime: scores are not defined for this slice" if (!defined($f_scores));
	
	foreach my $f_score (@$f_scores) {
	    if (defined $f_score->diff_score) {
		#printf("position %d observed %.4f expected %.4f difference %.4f\n",  $score->position, $score->observed_score, $score->expected_score, $score->diff_score);
		print FOUT
		    $f_score->observed_score.",";
	    }
	}
	print FOUT
	    "\n";

	
	
#thre prime 
	my $t_start = $end- $flanking_length;
	my $t_end   = $end + $flanking_length;
	#get slice
	my $t_slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $t_start,$t_end);
	throw("No Slice can be created with coordinates $chr:$t_start-$t_end") if (!$t_slice);
	my $t_display_size = $t_slice->end - $t_slice->start +1;
	my $t_scores       = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $t_slice, $t_display_size);
	print
	    "ThreePrime: scores are not defined for this slice" if (!defined($t_scores));
	

	foreach my $t_score (@$t_scores) {
	    if (defined $t_score->diff_score) {
	
		print TOUT
		    $t_score->observed_score.",";
	    }
	}
	print TOUT
	    "\n";

###core
	my $core_left_gap =  int(($end-$start - $core_length)/2);
	my $c_start = $start + $core_left_gap;
	my $c_end  =  $c_start + $core_length;
	#get slice
	my $c_slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $c_start,$c_end);
	throw("No Slice can be created with coordinates $chr:$c_start-$c_end") if (!$c_slice);
	my $c_display_size = $c_slice->end - $c_slice->start + 1; 
	my $c_scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $c_slice, $c_display_size);
	print
	    "Core: scores are not defined for this slice" if (!defined($c_scores));
	

	foreach my $c_score (@$c_scores) {
	    if (defined $c_score->diff_score) {
	
		print COUT
		    $c_score->observed_score.",";
	    }
	}
	print COUT
	    "\n";	

    }
    close(FOUT);
    close(TOUT);
    close(COUT);
    $gff_in->close();
    
				    }#GET_CONSERVATION_SCORES_OF_SEQS#

sub GET_CONSERVATION_SOCRES_OF_SEQUENCES($$$$){
	my $gff_input_file       = shift;
	my $path_to_data         = shift;
	my $flank                = shift;
	my $start_or_end         = shift;
	
	

	
	my $gff_input_file_full_name = $path_to_data."/".$gff_input_file;
	my @file_specs               = split(/\./,$gff_input_file );
	my $output_file_full_name    = $path_to_data."/Conservation_Scores_".$start_or_end."_".$file_specs[0].".txt";
	
	open(OUT,">$output_file_full_name") or die "Cannot open file $output_file_full_name to write in!\n";
	
	my $gff_in = Bio::Tools::GFF->new(-file =>$gff_input_file_full_name, -version =>1 );
	
	my $feature;
	while($feature = $gff_in->next_feature()){
		my $start = $feature->start();
		my $end   = $feature->end();
		my $center;
		if($start_or_end eq 'start'){
			$center = $start;
		}
		elsif($start_or_end eq 'end'){
			$center = $end;
		}
		else{
			print
				"Unknown value for start_or_end!\n";
				exit;
		}
		my $new_start = $center- $flank;
		my $new_end   = $center+ $flank;
		my $slice_length = $new_end-$new_start;
		my $chr   = $feature->seq_id();
		$chr     =~ s/chr//;
		my $source_tag = $feature->source_tag();
		#get slice
		my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $new_start,$new_end);
		throw("No Slice can be created with coordinates $chr:$new_start-$new_end") if (!$slice);
		#get conservation adaptor
		my $cs_adaptor = $registry->get_adaptor("Multi", 'compara', 'ConservationScore');
		my $display_size = $slice->end - $slice->start + 1; 
		my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $display_size);
		print
			"scores are not defined for this slice" if (!defined($scores));
		print
			"looking at conservation of sequence : Chr = ". $chr . " start= ". $slice->start. " end= ".$slice->end. "\n";
		
		foreach my $score (@$scores) {
   			 if (defined $score->diff_score) {
				#printf("position %d observed %.4f expected %.4f difference %.4f\n",  $score->position, $score->observed_score, $score->expected_score, $score->diff_score);
				print OUT
					$score->observed_score.",";
    			}
		}
		print OUT
			"\n";
	}
	close(OUT);
	$gff_in->close();
	   
}#GET_CONSERVATION_SOCRES_OF_SEQUENCES#
