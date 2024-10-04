#!/usr/bin/perl -w

use strict;
use FileHandle;

unless($#ARGV == 2){
    die "uasge: $0 [5'-UTR (fasta file)] [CDS (fasta file)] [run_raccess_contra (binary file)]\n";
}

my $Raccess_bin = $ARGV[2];

my $utr = &readOneSeq($ARGV[0]);
my $cds = &readOneSeq($ARGV[1]);
#print("$utr\n");
#print("$cds\n");
my ($L_utr, $L_cds) = (length $utr, length $cds);
# check start codon and sequence length
{
    if($cds !~ /^A[TU]G/){
	die "The start codon must be ATG.\n";
    }
    if($L_utr < 20){
	die "The length UTR must be at least 20 nt.\n";
    }
    if($L_cds < 90){
	die "The length CDS must be at least 90 nt.\n";
    }
}

# calculate emopec score
my ($emopec_score, $cand_sd) = &calcEmopecScore($utr);
#print "$emopec_score\n";

# calculate accC which is the accessibility calculated by the Contrafold model
my $accC = 0;
{
    my $flseq = $utr . $cds;
    $accC = &calcAccC($flseq);
    
}

my $final_score = 9.0488809 + 0.2473070 * $emopec_score + 0.8546904 * $accC;

my $level = &getExpLevel($final_score);

print "accC\t\t: $accC\n";
print "emopec score\t: $emopec_score\n";
print "predicted SD\t: $cand_sd\n";
print "RBSeval score\t: $final_score\n";
print "Exp. level\t: $level\n";

sub getExpLevel{
    my $sc = $_[0];
    
    if($sc >= 9.79371871522648){
	return "Very high (Top 10 percentile)";
    }
    elsif($sc >= 9.22676640010229){
	return "High (Top 30 percentile)";
    }
    elsif($sc >= 8.67717639210721){
	return "Middle (Top 50 percentile)";
    }
    elsif($sc >= 8.13572274167719){
	return "Low (Bottom 50 percentile)";
    }
    else{
	return "Very low (Bottom 30 percentile)";
    }
}


sub calcAccC{
    my $seq = $_[0];
    
    my $tmp_seq_file = "./tmp$$.fa";
    my $tmp_out_file = "./tmp$$.out";
    
    my $ofh = new FileHandle(">$tmp_seq_file") || die;
    print $ofh ">tmp\n";
    print $ofh "$seq\n";
    $ofh->close();
    
    # check and run raccess
    if($Raccess_bin !~ /run_raccess_contrafold$/){
	die "Raccess binary file name must be run_raccess_contrafold.\n";
    }
    if(!-e $Raccess_bin){
	die "Raccess binary file was not found.\n";
    }
    
    #print STDERR "Run raccess program.\n";
    `$ARGV[2] -seqfile=$tmp_seq_file -outfile=$tmp_out_file -access_len=35`;
    
    my %ACC;
    {
	my $fh = new FileHandle("$tmp_out_file") || die;
	while(<$fh>){
	    chomp;
	    if(/^>/ || /^$/){
		# skip
	    }
	    elsif(/^\d+\t/){
		my ($fm, $inf) = split /\t/;
		$fm++; # convert to 0-based to 1-based
		foreach my $w_acc (split ';', $inf){
		    my ($w, $acc) = split ",", $w_acc;
		    my $to = $fm + $w - 1;
		    $ACC{$fm}{$to} = -$acc;
		}
		
	    }
	    else{
		die "Unexpected line in raccess result ($_)\n";
	    }
	}
	$fh->close();
    }
    
    my $atg = $L_utr + 1;
    my $l = ($atg - 31) + 12;
    my $r = ($atg - 31) + 46;
    
    # delete file
    unlink ($tmp_seq_file);
    unlink ($tmp_out_file);
    
    if(!defined $ACC{$l}{$r}){
	die "accC value could not be calculated.\n";
    }
    
    return $ACC{$l}{$r};
}


sub readOneSeq{
    my $seq;
    
    my $fh = new FileHandle($_[0]) || die;
    
    chomp ($_ = <$fh>);
    if(!/^>/){
	die "Unexpected fasta header $_\n";
    }
    while(<$fh>){
	chomp;
	if(/^>/){
	    die "Multifasta format is prohibitted ($_)\n";
	}
	$seq .= $_;
    }
    $fh->close();
    
    return $seq;
}


sub calcEmopecScore{
    my $my_leader = $_[0];
    
    my $max     = -1e10;
    my $max_sd  = "";
    my $max_ds  = -1e10;
    #my $max_raw = -1e10;
    for(my $i = 1; $i <= 10; $i++){
	
	my $sd  = substr($my_leader, -6 - $i, 6);
	my $val = `python ./calc_emopec_score.py $sd $i`;
	chomp($val);
	
	if($max < $val){
	    $max    = $val;
	    $max_sd = $sd;
	    $max_ds = $i; # no use
    
	    #$max_raw = $raw;
	    
	}
	
    }
    #print "$sd_id\t$max\t$max_raw\t$max_sd\t$max_ds\n";
    return ($max, $max_sd);

}

