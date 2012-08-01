#!/usr/bin/perl





use lib "/nfs/users/nfs_h/hk3/src/lib/perl5";
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Perl;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;


#my $gff_output_file = "dummy.gff";

my $dir = "/nfs/th_group/hk3/DHS_Human/UW_fdr0.01";
my $gff_output_file = "/lustre/scratch103/sanger/hk3/Analysis_Human_Pecan/DNaseI/Chr_22_Uni_Washington_DNaseI.gff";
#&get_length_of_features_in_a_gff_file("/nfs/th_group/hk3/DHS_Human/Union_of_Chr_22_Uni_Duke_DNaseI.gff", "/nfs/th_group/hk3/DHS_Human/Length_Union_of_DNaseI_Features_on_Chr_22_Duke.txt");
#&get_DNaseI($dir, $gff_output_file, "chr22", "UW");
&get_dnaseI_specific_to_a_chro_from_a_bed_file("/nfs/th_group/hk3/DHS_Human/multi-tissue.master.ntypes.simple.hg19.bed", "/nfs/th_group/hk3/DHS_Human/Chr22_combined_uw_and_duke_multi_tissue_master_ntypes_simple.gff", "chr22");
exit;


my $feature = new Bio::SeqFeature::Generic(-start => 10, -end=>100,-strand => 1,  
					       -primary => 'primary',
					       -source_tag =>'source_tag',
					       -seq_id => 'seqId',
					       -score => 567,
					       -tag => {
						   new      =>1,
						   author   => 'hashem',
						   silyytag => 'sillytag'
					       }
    );

my $feature_2 = new Bio::SeqFeature::Generic(-start => 10, -end=>100,-strand => 1,  
					     -primary => 'primary',
					     -source_tag =>'source_tag',
					     -seq_id => 'seqId',
					     -score => 567,
					     -tag => {
						 new      =>1,
						 author   => 'hashem',
						 silyytag => 'sillytag'
					     }
    );







      my $gffout = new Bio::Tools::GFF(-file => ">$gff_output_file" ,
-gff_version => 1);


$gffout->write_feature($feature);
$gffout->write_feature($feature_2);
exit;



########################################### subroutines ########################################



sub get_dnaseI_specific_to_a_chro_from_a_bed_file($$$){
    my $bed_file_name          = shift;
    my $gff_output_file_name   = shift;
    my $chr                    = shift;
    
    my $seq_id         = $chr;
    if($chr =~ m/^chr/){
	$seq_id =~ s/chr//;
    }
    open(BED, $bed_file_name ) or die "Cannot open bed file $bed_file_name !";
    my $gff_output = new Bio::Tools::GFF(-file => ">$gff_output_file_name" , -gff_version => 1);
    while(my $line = <BED>){
	if($line =~ m/^$chr/){
	    chomp($line);
	    my @specification = split("\t",$line);
	    
	    my $one_gff_feature = new Bio::SeqFeature::Generic(-seq_id   => $seq_id,
							       -score    => $specification[4],
							       -start    => $specification[1],
							       -end      => $specification[2],
							       -source_tag=> "DNaseI",
							       -primary    =>  "combined_Duke_and_UW" ,
							       
							       
		);
	    $gff_output ->write_feature($one_gff_feature);
	}
	
    }
    close(BED);
						  }#get_dnaseI_specific_to_a_chro_from_a_bed_file#

sub get_DNaseI ($$$$){
    my $path_to_DNaseI_files = shift;
    my $output_gff_file      = shift;
    my $chr                  = shift;
    my $primary              = shift;
    
    opendir(DIR, $path_to_DNaseI_files) || die "Can't open directory $path_to_DNaseI_files $!";
    my $gffout =  new Bio::Tools::GFF(-file => ">$output_gff_file" , -gff_version => 1);
    my $seq_id = $chr;
    if($chr =~ m/^chr/){
		    $seq_id =~ s/chr//;
		}

    
    
    while(my $file = readdir(DIR)){
	next if ($file =~ m/^\./);
	my $full_file_name = $path_to_DNaseI_files."/".$file;
	open(FILE, $full_file_name) or die "Cannot open  $full_file_name";
	while(my $line = <FILE>){
	    if($line =~ m/^$chr/){
		chomp($line);
		my @specifications = split("\t", $line);
		my $start = $specifications[1];
		my $end   = $specifications[2];
		if($start > $end ){
		    $start = $specifications[2];
		    $end   = $specifications[1];
		    
		}
		
		my $feature = new Bio::SeqFeature::Generic(-start      => $start,
							   -end        => $end,
							   -seq_id     => $seq_id,
							   -score      => 0,
							   -strand     => 1,
							   -source_tag => 'DNaseI',
							   -primary    =>  $primary ,
							  # -tag        => {
							   #    curator => 'hashem',
							   #    testign => 'testing',
							   #}
		    );
		$gffout->write_feature($feature);
	    }
	}
    }
    close(FILE);
    closedir(DIR);
}#get_DNaseI#


sub get_length_of_features_in_a_gff_file($$){
    my $a_gff_input_file     = shift;
    my $output_file          = shift;
    open(GFF, $a_gff_input_file) or die "Cannot open $a_gff_input_file!";
    open(OPF, ">$output_file" ) or die "Cannot open $output_file for writtin outputs!";
   
    
    while(my $line = <GFF>){
	
	chomp($line);
	my @spec = split("\t", $line);

	my $one_feature_start = $spec[3];
	my $one_feature_end   = $spec[4];
	my $one_feature_length= $one_feature_end - $one_feature_start;
	if($one_feature_length <  0){
	    
	    print
		"Found unusuall features! start is less than or equal to end!\n";
	}
	print OPF
	    $one_feature_length.",";
	    
    }
    close(GFF);
    close(OPF);
}#get_length_of_features_in_a_gff_file#
