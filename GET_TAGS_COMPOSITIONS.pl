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






# make connection to Ensembl Registry
my $registry = 'Bio::EnsEMBL::Registry';


#load registry
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
    );

########################### Note #######################################################
# this is to obtain the tags compositions (frequency of nucleotides at each position)  #
# it requires tag files at bed format (6 columuns), an output file name and offSet     #
# which is used as start-offSet for the start of a tag.                                #
########################################################################################

my ($species,$chr,$tag_length, $tags_input_file, $output_file_name, $offset, $repeatmask);
my $help = '';
if(!GetOptions('species|s=s'             => \$species,
	       'chr|c=s'                 => \$chr,
	       'tag_length|l=s'          => \$tag_length,
	        'tags_input_file|i=s'     => \$tags_input_file,
               'output_file_name|o=s'    => \$output_file_name,
                'offset|f=s'             => \$offset,
	       'repeatmaski|r=s'         => \$repeatmask,
                'help|h'                 => \$help) 
 || !(defined($tags_input_file) && defined($output_file_name) && defined($offset) ) || $help)
{
    print <<END_USAGE;
Usage:
$0  --species=species --chr=chr --tag_length=tag_length --tags_input_file=tags_input_file --output_file_name=output_file_name --offset=offset --repeatmask=repeatmask
$0 --help
--species             /-s species: must be recognized by ensembl for instance Human.
--chr                /-c chromosome: example: 22 and not chr22.
--tag_length         /-l length of tags.
--tags_input_file     /-i A bed file that contains tags.
--output_file_name   /-o  Output file name. this will be a matrix-like with 4 rows and offset+lengthOfTags columns
--offset            /-f  offset. for getting exactly tag compostions, put offset=0, otherwise you will get composition of sequences from start-ffset to end.
--repeatmask        /-r wether to repeatmask or not: true or false is accepted.
--help               /-h to print this text.
Example usage:
perl $0 -s species -c chr -l tag_length -i tags_input_file -o output_file_name -f 0 -r false
END_USAGE
exit(1);
}#if end#



my $chr_seq = &GET_CHROMOSOME_SEQUENCE($species,$chr,$repeatmask);
print 
   "chromosome sequence fetched\n";

&GET_COMPOSITION_OF_TAGS($chr_seq,$tags_input_file,$output_file_name,$tag_length,$offset);
exit;


########################### subroutines ###########


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


sub GET_UNIQUE_FEATURES_FROM_A_BED_FILE($){
    my $input_bed_file     = shift;


    open(IN ,$input_bed_file ) or die "Cannot open file $input_bed_file for reading data\n";
    my @lines = <IN>;
    close(IN);
    my %unique_features;
    foreach my $line( @lines){
	chomp($line);
	my @line_spec = split("\t", $line);
	my $chr    = $line_spec[0];
	my $start  = $line_spec[1];
	my $end    = $line_spec[2];
	my $source = $line_spec[3];
	my $score  = $line_spec[4];
	my $strand = $line_spec[5];
	my $one_object = {
	    Chr => $chr,
	    Start => $start,
	    End   =>$end,
	    Source => $source,
	    Score  => $score,
	    Strand => $strand
	};
	$unique_features{$start}= $one_object;
    }

     return %unique_features;
}#GET_UNIQUE_FEATURES_FROM_A_BED_FILE#





sub GET_COMPOSITION_OF_TAGS($$$$$){
    my $chro_seq = shift;
    my $input_file = shift;
    my $output     = shift;
    my $tag_length = shift;
    my $start_offset = shift;
    

    my %unique_tags = &GET_UNIQUE_FEATURES_FROM_A_BED_FILE($input_file);

    my $new_length = $tag_length + $start_offset; 
    #declare and assign zero to the matrix:
    my @com_matrix;
    for(my $r=0; $r< $new_length;$r++ ){
	for(my $c=0; $c<4; $c++){
	    $com_matrix[$r][$c]=0;
	}
    }

    foreach my$key(keys %unique_tags){
	my $one_feat = $unique_tags{$key};
	my $chr         = $one_feat->{Chr};
	$chr =~ s/chr//;
	my $start       = $one_feat->{Start};
	my $end         = $one_feat->{End};
	my $source_tag  = $one_feat->{Source};
	my $score       = $one_feat->{Score};
	my $strand      = $one_feat->{Strand};
	
	
	my $new_start;
        # if it is on positive strand, get back the same as offset length
	if($strand eq "+"){
	    $new_start = $start-$start_offset;
	}
	#if it is on negative starnad, start form real start but get more seq from three prime end
	elsif($strand eq "-"){
	    $new_start = $start;
	}
	else{
	    die "unkwon stradn at: $start \t $end \t $source_tag \t $score \t $strand \n";
	}
	
#	my $new_start = $start-$start_offset;
	my $seq_stretch = uc(substr($chr_seq,$new_start,$new_length));
	next if( $seq_stretch =~ m/N/);
	my $seq_length = length($seq_stretch);
	if($seq_length ne ($new_length)){
	    die "seq length  $seq_length is not matching with $new_length +1\n";
	}
	if($strand eq "-"){
	    $seq_stretch = &GET_REVCOMP($seq_stretch);
	}
	

	for(my $p=0; $p<$new_length;$p++){
	    my $one_base = substr($seq_stretch,$p,1);
	    if($one_base eq "A"){
		$com_matrix[$p][0]= $com_matrix[$p][0]+1;
	    }
	    elsif($one_base eq "C"){
		$com_matrix[$p][1]= $com_matrix[$p][1]+1;
	    }
	    elsif($one_base eq "G"){
		$com_matrix[$p][2]= $com_matrix[$p][2]+1;
	    }
	    elsif($one_base eq "T"){
		$com_matrix[$p][3]= $com_matrix[$p][3]+1;
	    }
	    else{
		die "expecting nucleotide but got $one_base!\n";
	    }
	}#for over $p#	
    }#foreach#

#print results into ouput file:
    open(OUT,">$output") or die "cannot open file to write in!\n";
    for(my $r=0;$r<$new_length; $r++){
	for(my $c=0; $c<3;$c++){
	    print OUT
		$com_matrix[$r][$c].",";
	}
	print OUT
	    $com_matrix[$r][3];
	print OUT
	    "\n";
    }


}#GET_COMPOSITION_OF_TAGS#




sub GET_REVCOMP($){
    my $dna_seq  = shift;
    
    my $revcom_seq = reverse $dna_seq;
    $revcom_seq =~ tr/ACGTacgt/TGCAtgca/;
    return $revcom_seq;
		}#GET_REVCOMP#

