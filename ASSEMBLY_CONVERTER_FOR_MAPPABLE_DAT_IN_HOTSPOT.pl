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

######################### Note: ############################################################################
# this will accept inputs files in bed format with three collumns: chr start end  and will write the   #
# converted data into a bed fromat again.                                                              #
#      #
#                                                            #
###########################################################################################################


use IO::File;
use Getopt::Long;

use Bio::EnsEMBL::Registry;

my $strand = 1;
my $cs_name = "chromosome";

my ( $filename, $species,$old_assembly, $outputfile );
my $help = '';

if ( !GetOptions( 'file|f=s'        => \$filename,
                  'species|s=s'     => \$species,
		  'old_assembly|o=s' => \$old_assembly,
		  'outputfile|r=s'  => \$outputfile,
                  'help|h!'         => \$help )
     || !( defined($filename) && defined($species) && defined($old_assembly) && defined($outputfile))
     || $help )
{
  print <<END_USAGE;

Usage:
  $0 --species=species --file=filename --odlAssembly=oldAssembly --outputFile=outputFile

  $0 --help


    --species / -s  Name of species.  Mappings are currently only
                    available for mouse and human.

    --file / -f     Name of file containing a list of features to map to
                    the most recent assembly.  Only bed format with three comlumns is currently accepted:
                    chr tab stat tab end

     --old_assembly / -o Name of old assembly: for example NCBI36.
    --outputFile / -r Name of outputFile(full path)
    --help    / -h  To see this text.

Example usage:

  $0 -s mouse -f slices.txt

END_USAGE

  exit(1);
} ## end if ( !GetOptions( 'file|f=s'...))




my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db( '-host' => 'ensembldb.ensembl.org',
                                  '-user' => 'anonymous' );

my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );

my $in = IO::File->new($filename);


if ( !defined($in) ) {
  die( sprintf( "Could not open file '%s' for reading", $filename ) );
}

my $out= IO::File->new(">$outputfile");
if(!defined($out)){
die( sprintf( "Could not open file '%s' for writing", $outputfile ))
}
while(my  $line= $in->getline()){
    chomp($line);
 #Skip commment lines
  next if($line =~ m/^#/);
  # Skip lines containing only whitespace.
  next if ($line =~ /^\s*$/);
    my @line_spec = split("\t", $line);
    my $old_chr = $line_spec[0];
    $old_chr =~ s/chr//;
    my $old_start = $line_spec[1];
    my $old_end   = $line_spec[2];

#get a slice from the old region:
    my $old_slice = $slice_adaptor->fetch_by_region($cs_name,$old_chr,$old_start,$old_end,$strand, $old_assembly);
  # Complete possibly missing info.
   $cs_name ||= $old_slice->coord_system_name();
  $old_chr ||= $old_slice->seq_region_name();
  $old_start   ||= $old_slice->start();
  $old_end     ||= $old_slice->end();
  $strand  ||= $old_slice->strand();
  $old_assembly ||= $old_slice->coord_system()->version();

#printf( "# %s\n", $old_slice->name() );

    foreach my $segment(@{ $old_slice->project('chromosome') } ){
	my $new_start =  $old_start + $segment->from_start() - 1;
	my $new_end  =   $old_start + $segment->from_end() - 1;
	my @slice_spec = split(":",$segment->to_Slice()->name());
	my $length = @slice_spec;

	my $one_line = $slice_spec[2]. "\t". $slice_spec[3]. "\t".$slice_spec[4]."\n"; 
	print $out
	    $one_line;
    

    }
}
$in->close();
$out->close();
exit;

