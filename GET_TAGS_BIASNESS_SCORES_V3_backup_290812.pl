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

########################### Note ###########################################################################################################################
# very similar to V2 BUT:                                                                                                                                  #
# We know that  some of the nucleotides just before starting(5prime side) positions of pve starnads are biased too.                                        #
#(and also some nucleotides just after the 3prime side of tags in negative strands). Now we think these nucleotides will                                  #
# help us to discriminate biased tags from non-biased better.                                                                                              #
# therefore for each given tag if it is in pve strand it will move to start-offset, end-offset, if it is in nve strand it will move                        #
#to start+offset, end+offset and then these new fake tags will be scored (compared to background model). note that finally we will write our original tag   #
# coordiantes into output files and not these fake coordiantes.                                                                                            #
############################################################################################################################################################

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
--output_file_name   /-o  Output file name. This is a bed file indeed it is input file just the scores are substituted with new biasness scores.
--offset            /-f  offset. for getting exactly tag compostions, put offset=0, otherwise you will get composition of sequences from start-ffset to end.
--repeatmask        /-r wether to repeatmask or not: true or false is accepted.
--help               /-h to print this text.
Example usage:
perl $0 -s species -c chr -l tag_length -i tags_input_file -o output_file_name -f 0 -r false
END_USAGE
exit(1);
}#if end#


my $prior_probability                    = 0.25;
my $chr_seq                              = &GET_CHROMOSOME_SEQUENCE($species,$chr,$repeatmask);
print 
  "chromosome sequence fetched\n";
my @matrix_of_compositions_of_real_tags  = &GET_PWM($chr_seq,$tags_input_file,$output_file_name,36,$offset,$prior_probability,"real",0);

my @matrix_of_compositions_of_moved_tags = &GET_PWM($chr_seq,$tags_input_file,$output_file_name,36,0,$prior_probability,"background",100);
&GET_BIASNESS_SCORE_FOR_TAGS_V2($chr_seq,$tags_input_file,$output_file_name,\@matrix_of_compositions_of_real_tags,\@matrix_of_compositions_of_moved_tags,36,$offset);
exit;




################################################################# subroutines ##############################################


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



sub GET_REVCOMP($){
    my $dna_seq  = shift;
    
    my $revcom_seq = reverse $dna_seq;
    $revcom_seq =~ tr/NACGTnacgt/NTGCAntgca/;
    return $revcom_seq;
		}#GET_REVCOMP#




sub GET_PWM($$$$$$$$){
    my $chro_seq           = shift;
    my $input_file         = shift;
    my $output             = shift;
    my $tag_length         = shift;
    my $offset             = shift;
    my $prior_prob         = shift;
    my $model               = shift;
    my $transforming_length = shift;

   

    $transforming_length = 0 unless defined($transforming_length);
    my %unique_tags = &GET_UNIQUE_FEATURES_FROM_A_BED_FILE_V2($input_file);
    my $number_seqs_with_ns =0;
   
    #declare and assign zero to the matrix:
    my @com_matrix;
    for(my $r=0; $r< $tag_length; $r++ ){
	for(my $c=0; $c<4; $c++){
	    $com_matrix[$r][$c]=0;
	}
    }
    my $number_of_unique_tags = keys  %unique_tags;
    print 
	"There were ". $number_of_unique_tags." unque tags.\n";
    foreach my$key(keys %unique_tags){
	my $one_feat = $unique_tags{$key};
	my $chr         = $one_feat->{Chr};
	$chr =~ s/chr//;
	my $start       = $one_feat->{Start};
	my $end         = $one_feat->{End};
	my $source_tag  = $one_feat->{Source};
	my $score       = $one_feat->{Score};
	my $strand      = $one_feat->{Strand};
	my $modified_start;
	my $modified_end;
	if($strand eq "+"){
	    $modified_start = $start - $offset;
	    $modified_end   = $end   - $offset;
	    $modified_start = $modified_start +   $transforming_length;
	    $modified_end   = $modified_end  +    $transforming_length;
	}
	elsif($strand eq "-"){
	    $modified_start  = $start + $offset;
	    $modified_end    = $end   + $offset;
	    $modified_start = $modified_start -  $transforming_length;
	    $modified_end   = $modified_end   -  $transforming_length;
	}
	else{
	    die "Unknown strand!\n";
	}
	my $seq_stretch = uc(substr($chr_seq,$modified_start,$tag_length));
	if( $seq_stretch =~ m/N/){
	    $number_seqs_with_ns++;
	    next;
	}
	my $seq_length = length($seq_stretch);
	if($seq_length ne ($tag_length)){
	    die "seq length  $seq_length is not matching with $tag_length +1\n";
	}
	if($strand eq "-"){
	    $seq_stretch = &GET_REVCOMP($seq_stretch);
	}
	for(my $p=0; $p<$tag_length;$p++){
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
    print 
	"there were $number_seqs_with_ns  of seqs containing at least one N!\n";
#make position wieght matrix
    for(my $r=0;$r< $tag_length; $r++){
	for(my $c=0; $c<4;$c++){
	    #note: ln { (n_ij+p_i)/(N+1)*p_i } is almost equal to ln( (f_ij)/p_i ) where p_i is prior probability and f_ij = n_ij/N. this is in fact a scaling up, to get
#rid of gettin ln of zero.
	    $com_matrix[$r][$c] = &log2(($com_matrix[$r][$c] +$prior_prob  ) / ( ($number_of_unique_tags+1)*$prior_prob )  );
	   #theoritically we shouldn have scores bigger than 2 if prior_prob is taken equal to 0.25
	}

    }
   
    if($model eq "real"){
	$output =~ s/.bed/_Real_PWM.txt/;
    }
    elsif($model eq "background"){
$output =~ s/.bed/_BG_PWM.txt/;
    }
    else{
	die " Only real and background is acceptable for model parameter!\n";
    }
 my $pwm_testing_file = $output;
    
    open(OUT,">$pwm_testing_file") or die "Cannot open $pwm_testing_file to write in!\n";
    for(my $r=0;$r< $tag_length; $r++){
	for(my $c=0; $c<3;$c++){
	    print OUT
		 $com_matrix[$r][$c].",";
	}
	print OUT
	    $com_matrix[$r][3]."\n";
    }
    close(OUT);
    return @com_matrix;
}#GET_PWM#

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
	if( !( ($chr) && ($start) && ($end) && ($score) && ($source) && ($strand)  )  ){
	    print
		$line."\n";
	    die "required six columns to be defined but in the above mentioned line that wasnt the case\n";
	}
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

sub GET_UNIQUE_FEATURES_FROM_A_BED_FILE_V2($){
    #Note: the uniqueness is defined based on chr, start, strand and source, ie if two features have same start and strand then they are considered equal. and this is the
    #difference between this subroutine and its version of 1 where the uniqueness was only based on start position of the feature.

    my $input_bed_file  = shift;

    open(IN, $input_bed_file) or die "Cannot open $input_bed_file $input_bed_file to read the data\n";
    my @lines = <IN>;
    close(IN);
    my %unique_features;
    foreach my $line(@lines){
	chomp ($line);
	my @line_spec = split("\t", $line);
	my $chr    = $line_spec[0];
	my $start  = $line_spec[1];
	my $end    = $line_spec[2];
	my $source = $line_spec[3];
	my $score  = $line_spec[4];
	my $strand = $line_spec[5];
	
	
	my $one_key = $chr.$strand.$start.$source;
	if( !( ($chr) && ($start) && ($end) && ($score) && ($source) && ($strand)  )  ){
	    print
		$line."\n";
	    die "required six columns to be defined but in the above mentioned line that wasnt the case\n";
	}
	my $one_object = {
	    Chr     => $chr,
	    Start   => $start,
	    End     => $end,
	    Source  => $source,
	    Score   => $score,
	    Strand  => $strand
	};
	$unique_features{$one_key}= $one_object;
    }
return %unique_features;
}#GET_UNIQUE_FEATURES_FROM_A_BED_FILE_V2#

sub GET_BIASNESS_SCORE_FOR_TAGS_V2($$$$$$$){
    my $chr_seq          = shift;
    my $input_file       = shift;
    my $output_file      = shift;
    my $pwm_ref          = shift;
    my $bg_pwm_ref       = shift;
    my $seqLength        = shift;
    my $offset           = shift;

    open(INPUT, $input_file) or die "Cannot open file $input_file to read the data\n";
    open(OUTPUT, ">$output_file") or die "Cannot open file for $output_file writing data into it!\n ";
    my @lines = <INPUT>;
    close(INPUT);

    my $number_seqs_with_ns =0;
    
    #get probability of a sequence given weight matrix model and write them into this temp file for debugging:
    foreach my $line(@lines){
	chomp($line);
	my @line_spec = split("\t", $line);
	my $chr       = $line_spec[0];
	my $start     =  $line_spec[1];
	my $end       =  $line_spec[2];
	my $source    =  $line_spec[3];
	my $score     =  $line_spec[4];
	my $strand    =  $line_spec[5];
	if( !( ($chr) && ($start) && ($end) && ($score) && ($source) && ($strand)  )  ){
	    print
		$line."\n";
	    die "required six columns to be defined but in the above mentioned line that wasnt the case\n";
	}

	#modify the start and end based on offset
	my $modified_start;
	my $modified_end;
	if($strand eq "+"){
	    $modified_start =  $start - $offset;
	    $modified_end   =  $end   - $offset;
	}
	elsif($strand eq "-"){
	    $modified_start = $start + $offset;
	    $modified_end   = $end   + $offset;
	}
	else{
	    die "Unknown strand!\n";
	}
	my $seq_length = $modified_end - $modified_start;
	my $one_seq =  uc(substr($chr_seq,$modified_start,$seq_length));
	if( $one_seq =~ m/N/){
	    $number_seqs_with_ns++;
	    next;
	}
	

	#orient the sequence if it is in negative strand 
	if($strand eq "-"){
	    $one_seq = &GET_REVCOMP($one_seq);
	}
	my $one_log_ratio =  &GET_LOG_RATIO($pwm_ref,$bg_pwm_ref, $one_seq,$seqLength);
	my $one_pms       =  &GET_PROBABILITY_OF_SEQUENCE_GIVEN_PWM($pwm_ref, $one_seq,$seqLength);
	my $one_pbs       =  &GET_PROBABILITY_OF_SEQUENCE_GIVEN_PWM($bg_pwm_ref, $one_seq,$seqLength);
        #print out into file BUT keep the original coordinates
	print OUTPUT
	    $chr."\t".$start."\t".$end."\t".$source."\t".$one_log_ratio."\t".$strand."\n";
    }
    print 
	"$number_seqs_with_ns of sequences ignored because they were containing Ns!\n";
    close(OUTPUT);
}#GET_BIASNESS_SCORE_FOR_TAGS_V2#




sub log2($){
    my $n = shift;

    my $L = (log($n))/(log(2));
    return $L;
}#log2#

sub GET_LOG_RATIO($$$$){
    my $pwm_ref             = shift;
    my $background_pwm_ref  = shift;
    my $seq                 = shift;
    my $seq_length          = shift;

   # $seq_length = length($seq) unless defined($seq_length);
   

    my $prob_seq_give_pwm      = &GET_PROBABILITY_OF_SEQUENCE_GIVEN_PWM($pwm_ref, $seq, $seq_length);
    my $porb_seq_give_bg_model = &GET_PROBABILITY_OF_SEQUENCE_GIVEN_PWM($background_pwm_ref, $seq, $seq_length);
    my $log_ratio = ($prob_seq_give_pwm)-$porb_seq_give_bg_model;
    
    return $log_ratio;
		  }#GET_LOG_RATIO#


sub GET_PROBABILITY_OF_SEQUENCE_GIVEN_PWM ($$$){
    my $pwm_ref       = shift;
    my $seq           = shift;
    my $seq_length    = shift;
    
    #dereference
    my @pwm        = @$pwm_ref;
   
   # $seq_length = length($seq) unless defined($seq_length);
    my $sum = 0;
    $seq     = uc($seq);
    for(my $r=0; $r<  $seq_length; $r++ ){
	my $one_base = substr($seq,$r,1);
	my $one_value;
	if ($one_base eq "A"){
	    $one_value = $pwm[$r][0];
	}
	elsif ($one_base eq "C"){
	    $one_value = $pwm[$r][1];
	}
	elsif ($one_base eq "G"){
	    $one_value = $pwm[$r][2];
	}
	elsif ($one_base eq "T"){
	    $one_value = $pwm[$r][3];
	}
	else {
	    print
		"r= ".$r . "seq = " . $seq."\n";
	    die "was expecting A or C or G or T but got $one_base\n ";
	}
	$sum = $sum + $one_value
    }
    return $sum;
}#GET_PROBABILITY_OF_SEQUENCE_GIVEN_PWM#





