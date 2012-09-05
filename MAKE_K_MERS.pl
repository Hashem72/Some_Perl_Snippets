#!/usr/bin/perl
use strict;
use warnings;


my $seq = "accgtaccgtacc";
my $word = "accgtacc";
my  $word_length = length($word);

my $count = 0;
while($seq =~ /$word/g){
    pos $seq -= $word_length-1;
    $count++;
}

print
    "number of occurences is: " .$count. "\n";
my @eight_mers = &make_arrays_of_k_mers(8);

my $number_of_eight_mers = @eight_mers;
print
    "there are ". $number_of_eight_mers. " eight mers\n";
exit;


########################################################

sub make_arrays_of_k_mers($){
    my $k = shift;
my @bases = ('A','C','G','T');
my @words = @bases;

for my $i (1..$k-1){
    undef my @newwords;
    foreach my $w (@words){
         foreach my $b (@bases){
             push (@newwords,$w.$b);
         }
    }
    undef @words;
    @words = @newwords;
 }

    
    return @words;
}#make_arrays_of_k_mers#






   
