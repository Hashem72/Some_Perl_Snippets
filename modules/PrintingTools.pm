package PrintingTools;
use strict;
use warnings;
############################ subroutines #######################
sub print_out_a_matrix($){
    my $matrix = shift;
    
    my $numbero_of_rows = @{$matrix};
    my $number_of_cols =  @{$matrix->[0]};
    for (my $r = 0; $r < $numbero_of_rows; $r++) {
	for (my $c = 0; $c < $number_of_cols; $c++) {
	    printf("%.5f ", $matrix->[$r][$c]);
	}
	print
	    "\n";
    }
}#print_out_a_matirx#

sub print_out_a_matrix_into_a_file($$$){
    my $matrix      = shift;
    my $output_file = shift;
    my $seperator   = shift;

    my $nrows = @{$matrix};
    my $ncols = @{$matrix->[0]};
    open(OUT, ">$output_file") or die "cannot open file $output_file to write in!\n";
    for (my $r = 0; $r < $nrows; $r++) {
	for (my $c = 0; $c < $ncols; $c++) {
	    print OUT
		$matrix->[$r][$c];
	    unless($c eq $ncols-1){
		print OUT 
		    $seperator;
	    }
	}
	print OUT "\n";
    }
    close(OUT);
				   }#print_out_a_matrix_into_a_file#

sub print_out_freq_in_each_col($){
    my $matrix      = shift;
    
    my $number_of_rows = @{$matrix};
    my $number_of_cols = @{$matrix->[0]};
    
    for (my $c = 0; $c < $number_of_cols; $c++) {
	my $sum = 0;
	for (my $r = 0; $r < $number_of_rows; $r++) {
	    $sum = $sum + $matrix->[$r][$c];
	}
	print
	    "col = $c\t sum_of_freqs = $sum\n";
    }
 }#print_out_freq_in_each_col#

sub print_out_freq_in_each_row($){
    my $matrix    = shift;

    my $number_of_rows = @{$matrix};
    my $number_of_cols = @{$matrix->[0]};
    for (my $r = 0; $r < $number_of_rows; $r++) {
	my $sum = 0;
	for (my $c = 0; $c < $number_of_cols; $c++) {
	    $sum = $sum + $matrix->[$r][$c];
	}
	print
	    "row = $r\t sum_of_freqs = $sum\n";
    }

}# print_out_freq_in_each_row#

1;

