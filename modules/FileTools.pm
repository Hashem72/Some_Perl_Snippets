package FileTools;
use strict;
use warnings;


############################ subroutine #################


sub get_seq_from_fasta_file($$){
    my $fasta_file_name   = shift;
    my $seq_id            = shift;

    my $seqin = Bio::SeqIO->new('-format'=> 'fasta', -file => "$fasta_file_name");
    my $sequence;
    while(my $seqobj = $seqin->next_seq()){
	my $sequence_id = $seqobj->id();
	if($sequence_id eq $seq_id ){
	     $sequence = $seqobj->seq();
	    last;
	}
    }
    return $sequence;
}#get_seq_from_fasta_file#

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
	if( !( defined($chr) && defined($start) && defined($end) && defined($score) && defined($source) && defined($strand)  )  ){
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

sub read_in_a_matrix_from_a_file($){
    my $input_file_name = shift;

    my @matrix = ();
    my $matrix_ref = \@matrix;
    
    open(FILE, "<$input_file_name") or die "cannot open file $input_file_name to read the data!\n";
    while (<FILE>) {
	chomp $_;
	my @row = split(",");
	push(@matrix, [@row]);
    }
    close(FILE);
    return $matrix_ref;
    
				 }#read_in_a_matrix_from_a_file#


1;
