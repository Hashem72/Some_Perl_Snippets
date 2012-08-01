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
# Given a bed file (six columns: chr start end sourcename score strand ) this script   #
# will calculate the biasness score which log2(p(S|M)/p(S|B)) and then will substitute #
# old score (which is a constant score equal to 25 in DHS tags) with this new biasness #
#score. the output is again a bed file.                                                #
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
--output_file_name   /-o  Output is a bed file very similar to the input file. just scores are substituted with new biasness scores.
--offset            /-f  offset. for getting exactly tag compostions, put offset=0, otherwise you will get composition of sequences from start-ffset to end.
--repeatmask        /-r wether to repeatmask or not: true or false is accepted.
--help               /-h to print this text.
Example usage:
perl $0 -s species -c chr -l tag_length -i tags_input_file -o output_file_name -f 0 -r false
END_USAGE
exit(1);
}#if end#


my $prior_probability   = 0.25;
my $chr_seq = &GET_CHROMOSOME_SEQUENCE($species,$chr,$repeatmask);
print 
  "chromosome sequence fetched\n";

my @matrix_of_compositions_of_tags = &GET_PWM($chr_seq,$tags_input_file,$output_file_name,36,$offset,$prior_probability);

#get an average of weights for A, C, G and T that occur from position 10 onwards
my @averages = (0,0,0,0);
my $start =20;
my $end   = 36;
for(my $c=0; $c<4; $c++){
    my $sum_over_one_col = 0;
    for(my $r=$start; $r<$end; $r++){
	$sum_over_one_col = $sum_over_one_col + $matrix_of_compositions_of_tags[$r][$c];
    }
    
#nomlize over length
    $sum_over_one_col = $sum_over_one_col/($end-$start);
    $averages[$c]= $sum_over_one_col;
    print
	"c = ". $c ."prior is equal to ".  $averages[$c] . "\n";
}
&GET_BIASNESS_SCORE_FOR_TAGS($chr_seq,$tags_input_file,$output_file_name,\@matrix_of_compositions_of_tags,$averages[0],$averages[1],$averages[2],$averages[3],10);
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



sub GET_BIASNESS_SCORE_FOR_TAGS($$$$$$$$$){
    my $chr_seq          = shift;
    my $input_file       = shift;
    my $output_file      = shift;
    my $pwm_ref          = shift;
    my $prob_a           = shift;
    my $prob_c           = shift;
    my $prob_g           = shift;
    my $prob_t           = shift;
    my $seqLength       = shift;

    open(INPUT, $input_file) or die "Cannot open file $input_file to read the data\n";
    open(OUTPUT, ">$output_file") or die "Cannot open file for $output_file writing data into it!\n ";
    my @lines = <INPUT>;
    close(INPUT);
    
    #get probability of a sequence given weight matrix model and write them into this temp file for debugging:
    #probability of seq given model;
    my $pms_file = "/nfs/th_group/hk3/UW_DNaseI_HS/K562/pms_scores.txt";
    #prob of model given background
    my $pbs_file = "/nfs/th_group/hk3/UW_DNaseI_HS/K562/pbs_scores.txt";

    open(PMS,">$pms_file") or die "Cannot open file $pms_file to write in!\n";
    open(PBS, ">$pbs_file") or die "Cannot open file $pbs_file to write in!\n";
    
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
	my $seq_length = $end - $start;
	my $one_seq =  uc(substr($chr_seq,$start,$seq_length));
	#orient the sequence if it is in negative strand 
	if($strand eq "-"){
	    $one_seq = &GET_REVCOMP($one_seq);
	}
	
	my $one_log_ratio = &GET_LOG_RATIO($pwm_ref, $one_seq, $prob_a,$prob_c,$prob_g,$prob_t,$seqLength);
	my $one_pms      = &GET_PROBABILITY_OF_SEQUENCE_GIVEN_PWM($pwm_ref, $one_seq,$seqLength);
	#my $one_pbs      = &GET_PROBABILITY_OF_SEQUENCE_GIVEN_BACKGROUND_MODEL($one_seq, $prob_a,$prob_c,$prob_g,$prob_t,$seqLength);
	my $one_pbs = &GET_PROBABILITY_OF_SEQUENCE_GIVEN_BACKGROUND_MODEL_V2($one_seq, $prob_a,$prob_c,$prob_g,$prob_t,$seqLength);
	print PMS
	    $one_pms.",";
	print PBS
	    $one_pbs.",";
	print OUTPUT
	    $chr."\t".$start."\t".$end."\t".$source."\t".$one_log_ratio."\t".$strand."\n";
    }
    close(PMS);
    close(PBS);	
    close(OUTPUT);
}#GET_BIASNESS_SCORE_FOR_TAGS#


sub GET_PWM($$$$$$){
    my $chro_seq = shift;
    my $input_file = shift;
    my $output     = shift;
    my $tag_length = shift;
    my $start_offset = shift;
    my $prior_prob   = shift;
    

    my %unique_tags = &GET_UNIQUE_FEATURES_FROM_A_BED_FILE($input_file);

    my $new_length = $tag_length + $start_offset; 
    #declare and assign zero to the matrix:
    my @com_matrix;
    for(my $r=0; $r< $new_length;$r++ ){
	for(my $c=0; $c<4; $c++){
	    $com_matrix[$r][$c]=0;
	}
    }
    my $number_of_unique_tags = keys  %unique_tags;
    foreach my$key(keys %unique_tags){
	my $one_feat = $unique_tags{$key};
	my $chr         = $one_feat->{Chr};
	$chr =~ s/chr//;
	my $start       = $one_feat->{Start};
	my $end         = $one_feat->{End};
	my $source_tag  = $one_feat->{Source};
	my $score       = $one_feat->{Score};
	my $strand      = $one_feat->{Strand};
	
	

	my $new_start = $start-$start_offset;
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

#make position wieght matrix
    for(my $r=0;$r< $new_length; $r++){
	for(my $c=0; $c<4;$c++){
	    #note: ln { (n_ij+p_i)/(N+1)*p_i } is almost equal to ln( (n_ij)/p_i ) where p_i is prior probability. this is in fact a scaling up, to get
#rid of gettin ln of zero.
	    $com_matrix[$r][$c] = &log2(($com_matrix[$r][$c] +$prior_prob  ) / ( ($number_of_unique_tags+1)*$prior_prob )  );
	   #theoritically we shouldn have scores bigger than 2 if prior_prob is taken equal to 0.25
	}

    }
#    my $pwm_testing_file = "/nfs/th_group/hk3/UW_DNaseI_HS/K562/pwm_testing_file.txt";
#    open(OUT,">$pwm_testing_file") or die "Cannot open pwm_testing_file to write in!\n";
#    for(my $r=0;$r< $new_length; $r++){
#	for(my $c=0; $c<3;$c++){
#	    print OUT
#		 $com_matrix[$r][$c].",";
#	}
#	print OUT
#	    $com_matrix[$r][3]."\n";
#    }
#    close(OUT);
    return @com_matrix;
}#GET_PWM#

sub log2(){
    my $n = shift;
    my $L = (log($n))/(log(2));
    return $L;
}#log2#

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

sub GET_PROBABILITY_OF_SEQUENCE_GIVEN_BACKGROUND_MODEL($$$$$$){
    my $seq        = shift;
    my $pro_of_A   = shift;
    my $pro_of_C   = shift;
    my $pro_of_G   = shift;
    my $pro_of_T   = shift; 
    my $seq_length = shift;
    
    my $sum = 0;
   # $seq_length = length($seq) unless defined($seq_length);
    for(my $p=0;$p<$seq_length; $p++){
	my $one_base = substr($seq,$p,1);
	my $one_value;
	if($one_base eq "A"){
	    $one_value = &log2($pro_of_A);
	}
	elsif ($one_base eq "C"){
	    $one_value = &log2($pro_of_C);
	}
	elsif ($one_base eq "G"){
	    $one_value = &log2($pro_of_G);
	}
	elsif ($one_base eq "T"){
	    $one_value = &log2($pro_of_T);
	}
	else {
	    die "Expecting A or C or G or T but got $one_base\n";
	}
	$sum = $sum+ $one_value;
    }
    return $sum;
    
						       }# GET_PROBABILITY_OF_SEQUENCE_GIVEN_BACKGROUND_MODEL#

sub GET_PROBABILITY_OF_SEQUENCE_GIVEN_BACKGROUND_MODEL_V2($$$$$$){
    my $seq                = shift;
    my $log_of_pro_of_A    = shift;
    my $log_of_pro_of_C    = shift;
    my $log_of_pro_of_G    = shift;
    my $log_of_pro_of_T    = shift; 
    my $seq_length         = shift;

    my $sum    = 0;
    for(my $p=0; $p<$seq_length; $p++){
	my $one_base = substr($seq,$p,1);
	my $one_value;
	if($one_base eq "A"){
	    $one_value = $log_of_pro_of_A;
	}
	elsif($one_base eq "C"){
	    $one_value = $log_of_pro_of_C; 
	}
	elsif($one_base eq "G"){
	    $one_value = $log_of_pro_of_G;
	}
	elsif($one_base eq "T"){
	    $one_value = $log_of_pro_of_T;
	}
	else{
	    die "Expecting A or C or G or T but got $one_base!\n";
	}
	$sum = $sum + $one_value;
    }
    
    return $sum;
}# GET_PROBABILITY_OF_SEQUENCE_GIVEN_BACKGROUND_MODEL_V2#

sub GET_LOG_RATIO($$$$$$$){
    my $pwm_ref    = shift;
    my $seq        = shift;
    my $prob_a     = shift;
    my $prob_c     = shift;
    my $prob_g     = shift;
    my $prob_t     = shift;
    my $seq_length = shift;

   # $seq_length = length($seq) unless defined($seq_length);
   

    my $prob_seq_give_pwm = &GET_PROBABILITY_OF_SEQUENCE_GIVEN_PWM($pwm_ref, $seq, $seq_length);
   # my $porb_seq_give_bg_model = &GET_PROBABILITY_OF_SEQUENCE_GIVEN_BACKGROUND_MODEL($seq, $prob_a,$prob_c,$prob_g,$prob_t, $seq_length);
    my $porb_seq_give_bg_model = &GET_PROBABILITY_OF_SEQUENCE_GIVEN_BACKGROUND_MODEL_V2($seq, $prob_a,$prob_c,$prob_g,$prob_t, $seq_length); 
    my $log_ratio = ($prob_seq_give_pwm)-$porb_seq_give_bg_model;
    
    return $log_ratio;
		  }#GET_LOG_RATIO#

sub GET_REVCOMP($){
    my $dna_seq  = shift;
    
    my $revcom_seq = reverse $dna_seq;
    $revcom_seq =~ tr/ACGTacgt/TGCAtgca/;
    return $revcom_seq;
		}#GET_REVCOMP#





