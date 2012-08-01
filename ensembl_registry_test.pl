#!/usr/bin/perl

use strict;
use warnings;
use lib '/nfs/users/nfs_h/hk3/src/bioperl-live';
use lib '/nfs/users/nfs_h/hk3/src/ensembl/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-compara/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-variation/modules';
use lib '/nfs/users/nfs_h/hk3/src/ensembl-functgenomics/modules';
use Bio::EnsEMBL::Registry;

my $registry = "Bio::EnsEMBL::Registry";


$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);


# list of all Ensembl databases installed on a given database host: 
#my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

#foreach my $db_adaptor (@db_adaptors) {
 #   my $db_connection = $db_adaptor->dbc();

 #   printf(
  #      "species/group\t%s/%s\ndatabase\t%s\nhost:port\t%s:%s\n\n",
   #     $db_adaptor->species(),   $db_adaptor->group(),
    #    $db_connection->dbname(), $db_connection->host(),
     #   $db_connection->port()
   # );
#}

#exit;



## Load the databases into the registry and print their names
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
   );

my $slice_adaptor     = $registry->get_adaptor('Human','Core', 'Slice');
my $slice             = $slice_adaptor->fetch_by_region('chromosome', '22');
my $Atts              = $slice->get_all_Attributes();
foreach my $Att(@{$Atts}){
    print
	"Att= ", $Att->name(), "\n";
}#foreach#


exit;
