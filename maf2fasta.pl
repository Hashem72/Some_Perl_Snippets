#########################################################
### Hashem Koohy 12-10-2011
### this is meant to read a maf file and write blocks of alinments from that maf file into fasta files (one fasta file for each alignment block) ####
### Please note: this parser will ignore ancestral sequences
### Please note: if you have an emf file, then first use emf2maf.pl to convert it into a maf file and then use this convertor in order to write sequences into a fasta file
### in this verision, we will obviously loose some information for instance we dont care about alignment scores, strands and number of gaps in sequences. But if required it is easy to implement extracting these information too.
### Alignment blocks with duplicated sequences are ignored here!
###########################################################



#!/usr/bin/perl





use lib "/nfs/users/nfs_h/hk3/src/lib/perl5";
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Perl;

use constant false => 0;
use constant true  => 1;


my $VERSION = "1.0";

if (!@ARGV) {
  print STDERR qq"
  This is maf2fasta v$VERSION - MAF to FASTA file converter.

  Use: perl maf2fasta.pl mafInputFile.maf /path/to/fasta/directory/ ChrName

  Fasta Files will named as Chr_start_end.fasta where start and end corresponds to genomic coordinates of the core sequences (mainly first seq) in genome before alignment.

";
exit;
}



 
 my $inputMafFile       = $ARGV[0];
 my $pathForFastaFiles  = $ARGV[1];
 my $Chr                = $ARGV[2];

my @monkies = ('pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii','macaca_mulatta');
 





# species in 12-way alignment
&convert_maf_to_fast($inputMafFile,$pathForFastaFiles,$Chr, \@monkies, false);
my @species = ('homo_sapiens', 'pan_troglodytes','gorilla_gorilla', 'pongo_abelii','macaca_mulatta', 'callithrix_jacchus','mus_musculus', 'rattus_norvegicus','bos_taurus','sus_scrofa','equus_caballus','canis_familiaris' );
exit;

#####################################  subroutines ###############

sub convert_maf_to_fast($$$$){
    my $input_maf_file                   = shift;
    my $path_to_fasta_files              = shift;
    my $chr                              = shift;
    my $list_of_monkies_to_be_filter_ref = shift;
    my $filter_monkies                   = shift;
    

    my @species_to_filter = @$list_of_monkies_to_be_filter_ref;
    my $alignio = Bio::AlignIO->new(-file => $input_maf_file , , -format => 'maf');
    my $numAli =0;
    my $number_of_not_duplicated_bocks =0;
    while(my $aln = $alignio->next_aln()){
    my @seq_Ids ;
    $numAli++;
    
    
    my $firstSeqId = $aln->get_seq_by_pos(1)->display_id();
    my $isHumanSeq = ($firstSeqId =~ m/"hg16"/);
    if($firstSeqId){
	my $fistSeqStart = $aln->get_seq_by_pos(1)->start()-1;
	my $firstSeqEnd  = $aln->get_seq_by_pos(1)->end();
	my $outputFileName = $path_to_fasta_files.$chr."_".$fistSeqStart."_".$firstSeqEnd.".fasta";
	my $seqout = Bio::SeqIO->new('-foramat'=> 'Fasta', -file => ">$outputFileName");
	
	
	foreach my $seq($aln->each_seq()){
	    
	    $seq    = &set_seqId($seq);
	    my $seqId = $seq->display_id();
	    if($filter_monkies){#if asked for, filter monkies
		if(grep{$_ eq $seqId}   @species_to_filter ){
		    next;
		}
	    }
	    my $isAncestor  = ($seqId =~ m/ancestral_sequences/); 
	    my $strand     = $seq->strand();
	    unless($isAncestor){#ignore ancestor sequences
		if(grep{$_ eq $seqId}  @seq_Ids){#ingore duplicated sequence (if for instance we have two or more  segments of human, only the first segment is written into the fasta file.)
		    # print
		#	"seq is: ", $seqId,"\n";
		    next;
		}
		else{
		    $seqout ->write_seq($seq);
		    push(@seq_Ids,$seqId);
		    #    print
		    #"the seq added to array is", $seqId,"\n";
		}
		
	    }
	    
	}

	&delete_fasta_files_with_only_one_seqs($outputFileName);
    }
    
}#while#
    
    print 
	"total number of alignment blocks : ", $numAli, "\n"; 
			}#convert_maf_to_fast#

 
sub delete_fasta_files_with_only_one_seqs($){#after filtering duplications some fasta files are left with only one sequence, this is to delete them.
    my $file = shift;
    
    my $seqin = Bio::SeqIO->new('-format' => 'Fasta', -file => "$file");
    my $numb_of_seqs =0;
    while(my $seq_obj = $seqin->next_seq()){
	$numb_of_seqs++;
    }
    if($numb_of_seqs < 2){
	unlink($file);
	print
	    "file $file deleted because included only one sequence\n";
}
print 
    "overall $numb_of_seqs deleted!\n";
}#delete_fasta_files_with_only_one_seqs#

 
 sub get_approptiate_seq_id($){
 	my $seq = shift;
 	
 	
 	my $seqId = $seq->display_id();
 	my @id_split = split( /\./, $seqId);
 	
 	
 	my $newId    = $id_split[0];
 	my $finalId  = undef;

 	if($newId =~ m/homo_sapiens/){
 		$finalId = "Hsap";
 	}
 	elsif($newId =~ m/pan_troglodytes/){
 		$finalId = "Ptro";
 	}
 	elsif($newId =~ m/gorilla_gorilla/){
 		$finalId = "Ggor";
 	}
 	elsif($newId =~ m/pongo_abelii/){
 		$finalId = "Pabe";
 	}
 	elsif($newId =~ m/bos_taurus/){
 		$finalId = "Btau";	
 	}
 	elsif($newId =~ m/equus_caballus/){
 		$finalId = "Ecab";
 	}
 	elsif($newId =~ m/macaca_mulatta/){
 		$finalId = "Mmul";
 	}
 	elsif($newId =~ m/callithrix_jacchus/){
 		$finalId = "Cjac";
 	}
 	elsif($newId =~ m/equus_caballus/){
 		$finalId = "Ecab";
 	}
 	elsif($newId =~ m/mus_musculus/){
 		$finalId = "Mmus"
 	}
 	elsif($newId =~ m/rattus_norvegicus/){
 		$finalId = "Rnor";
 	}
 	elsif($newId =~ m/sus_scrofa/){
 		$finalId = "Sscr";
 	}
 	elsif($newId =~ m/ancestral_sequences/){
 		$finalId = "Aseq";
 	}
 	elsif($newId =~ m/canis_familiaris/){
 		$finalId = "Cfam";
 	}
 	else{
 		die "found unknown species: $newId \n";
 	}
 	$seq->display_id($finalId);
 	return $seq;
 }#get_approptiate_seq_id#
 
 sub set_seqId($){
 	my $seq = shift;
 	
 	my $seqId = $seq->display_id();
 	my @id_split = split( /\./, $seqId);
 	
 	
 	my $newId    = $id_split[0];
 	$seq->display_id($newId);
 	return $seq;
 }#set_seqId#
 
 sub is_there_duplicated_sequences($){
 	my $alignment_block  = shift;
 	
 	my $isDuplicated = false;
 	
 	my $nhomo        = 0;
 	my $npan         = 0;
 	my $ngorilla     = 0;
 	my $npongo       = 0;
 	my $nmacaca      = 0;
 	my $callithrix   = 0;
 	my $nmus         = 0;
 	my $nrattus      = 0;
 	my $nbos         = 0;
 	my $nsus         = 0;
 	my $nequus       = 0;
 	my $ncanis       = 0;

	my $number_of_duplicated_blocks =0;
 	foreach my $seq($alignment_block->each_seq()){ 
 		my $id = $seq->display_id();
 		if($id =~ m/^homo_sapiens/){
 			$nhomo++;
 		}
 		if($id =~ m/^pan_troglodytes/ ){
 			$npan++;
 		}
 		if($id =~ m/^gorilla_gorilla/){
 			$ngorilla++;
 		}
 		if($id =~ m/^pongo_abelii/){
 			$npongo++;
 		}
 		if($id =~ m/^macaca_mulatta/){
 			$nmacaca++;
 		}
 		if($id =~ m/^callithrix_jacchus/){
 			$callithrix++;
 		}
 		if($id =~ m/^macaca_mulatta/){
 			$nmus++;
 		}
 		if($id =~ m/^rattus_norvegicus/){
 			$nrattus++;
 		}
 		if($id =~ m/^bos_taurus/){
 			$nbos++;
 		}
 		if($id =~ m/^sus_scrofa/){
 			$nsus++;
 		}
 		if($id =~ m/^equus_caballus/){
 			$nequus++;
 		}
 		if($id =~ m/^canis_familiaris/){
 			$ncanis++;
 		}
 	}
 	
 	if( $nhomo>1 || $npan>1 || $ngorilla>1 || $npongo>1  || $nmacaca>1 || $callithrix>1 ||  $nmus>1 || $nrattus>1   ||  $nbos>1 || $nsus>1 ||  $nequus>1 ||  $ncanis>1){
 		$isDuplicated = true;
 		
 	}
 	return $isDuplicated;
 }#is_duplicated#
