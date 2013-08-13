#!/bin/perl -w

use Bio::SeqIO;

$in = Bio::SeqIO->new(-file=> "test.fastq", -format =>"fastq"

                          );

$out = Bio::SeqIO->new(-file=>">test", -format =>"fasta");

while ( my $seq = $in->next_seq() ) {$out->write_seq($seq); }
