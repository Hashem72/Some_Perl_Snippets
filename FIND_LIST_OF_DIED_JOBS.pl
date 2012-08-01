#!/usr/bin/perl


my $file_containing_list_of_successful_jobs = "/lustre/scratch103/sanger/hk3/Analysis_Human_Pecan/COMP_22/list_of_comp_files";
my $file_cotaining_list_of_all_jobs         = "/lustre/scratch103/sanger/hk3/Human_Data/Pecan_Alignments/22_List_of_Alignment_Files.txt";

my @failed_jobs = &GET_LIST_OF_FAILED_JOBS($file_containing_list_of_successful_jobs,$file_cotaining_list_of_all_jobs);
&WRITE_LIST_INTO_A_FILE(\@failed_jobs, "/lustre/scratch103/sanger/hk3/Human_Data/Pecan_Alignments/22_List_of_Alignment_Files_failed_jobs.txt");
exit;

###################################################

sub GET_LIST_OF_FAILED_JOBS($$){
    my $done_jobs = shift;
    my $all_jobs  = shift;
    

    my @list_of_successful_jobs = &GET_LIST_FROM_A_FILE($done_jobs);
    my @list_of_all_jobs        = &GET_LIST_FROM_A_FILE($all_jobs);
    my @failed_jobs ;

    my %done_jobs = map{$_ =>1}@list_of_successful_jobs;
    my $counter = 0;
    foreach(@list_of_all_jobs){
	if($done_jobs{$_}){
	    $counter++;
	    print
		$counter. "   ".$_ ."\n";
	}
	else{
	    push(@failed_jobs,$_);
	    print
	    "failed jobs ************ ". $_."\n";
	}
    }
    return @failed_jobs;

}#GET_LIST_OF_FAILED_JOBS#

sub GET_LIST_FROM_A_FILE($){
    my $inputFile = shift;
    my @list;
    open(FILE, $inputFile) or die "Cannot open file $inputFile!";
    while(my $line = <FILE>){
	chomp($line);
	push( @list, $line);
    }
    close(IN);
    return @list;
}#GET_LIST_FROM_A_FILE#


sub WRITE_LIST_INTO_A_FILE($$){
    my $array_of_failed_jobs_ref  = shift;
    my $outputFile                = shift;
    
    my @array_of_failed_jobs = @$array_of_failed_jobs_ref;
    open(OUTPUT, ">$outputFile") or die "Cannot open file $outputFile";
    foreach(@array_of_failed_jobs){
	print OUTPUT
	    $_."\n";
    }
    close(OUTPUT);
}#WRITE_LIST_INTO_A_FILE#
