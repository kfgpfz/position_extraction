use Getopt::Std;

getopt ('icpqo'); #i input sam, c chr.. , p snp pos, q min snp base quality, o outfile.txt


open (FD_IN, "$opt_i") || die "Error opening $opt_i";
open (FD_OUT, ">>$opt_o") || die "Error opening $opt_o";


@w = ();
%bases = ();

%bases1 = ();
%bases5 = ();
%bases10 = ();
%bases_1 = ();
%bases_5 = ();
%bases_10 = ();

%BadReads = ();


while($line = <FD_IN>) 
{
    if($line =~ /^@/)
    {
        next;
    }
    else
    {
	@w = split(/\t/, $line);
        $chr = $w[2];
        if($chr ne $opt_c){next;}
        
        $start = $w[3];
        if($opt_p > $start && $opt_p < $start + 70)
        {
           $offset_zerobased = $opt_p - $start;
            
           $nt1 = substr($w[9], $offset_zerobased + 1, 1);
	    $nt5 = substr($w[9], $offset_zerobased + 5, 1);
	    $nt10 = substr($w[9], $offset_zerobased + 10, 1);
	    
	    $nt_1 = substr($w[9], $offset_zerobased - 1, 1);
	    $nt_5 = substr($w[9], $offset_zerobased - 5, 1);
	    $nt_10 = substr($w[9], $offset_zerobased - 10, 1);
	    
           $bases1{$nt1}++;
	    $bases5{$nt5}++;
	    $bases10{$nt10}++;
	    
	    $bases_1{$nt_1}++;
	    $bases_5{$nt_5}++;
	    $bases_10{$nt_10}++;
        }
	
	$qline = $w[10];
	@w2 = split(//, $qline);
	$SNP_pos_qual = ord($w2[$offset_zerobased]) - 33;
	#print "SNPq  $SNP_pos_qual\n";
        if($SNP_pos_qual < $opt_q)
	{
	    $BadReads{$w[0]} = 1;
	}
        
    }
}
close(FD_IN);


$b1 = "";
$b5 = "";
$b10 = "";
$b_1 = "";
$b_5 = "";
$b_10 = "";

$max_val = 0;
foreach $key (keys %bases1)
{    
    if($bases1{$key} > $max_val)
    {
	$max_val = $bases1{$key};
	$b1 = $key;
    }
}

$max_val = 0;
foreach $key (keys %bases5)
{    
    if($bases5{$key} > $max_val)
    {
	$max_val = $bases5{$key};
	$b5 = $key;
    }
}

$max_val = 0;
foreach $key (keys %bases10)
{
    if($bases10{$key} > $max_val)
    {
	$max_val = $bases10{$key};
	$b10 = $key;
    }
}

$max_val = 0;
foreach $key (keys %bases_1)
{
    if($bases_1{$key} > $max_val)
    {
	$max_val = $bases_1{$key};
	$b_1 = $key;
    }
}

$max_val = 0;
foreach $key (keys %bases_5)
{
    if($bases_5{$key} > $max_val)
    {
	$max_val = $bases_5{$key};
	$b_5 = $key;
    }
}

$max_val = 0;
foreach $key (keys %bases_10)
{
    if($bases_10{$key} > $max_val)
    {
	$max_val = $bases_10{$key};
	$b_10 = $key;
    }
}


$discarded = 0;

open (FD_IN, "$opt_i") || die "Error opening $opt_i";

while($line = <FD_IN>) 
{
    if($line =~ /^@/)
    {
        next;
    }
    else
    {
	@w = split(/\t/, $line);
        $chr = $w[2];
        if($chr ne $opt_c){next;}
        
        $start = $w[3];
        if($opt_p > $start && $opt_p < $start + 70)
        {
            $offset_zerobased = $opt_p - $start;
            
            #check if the quality is acceptable
            if($BadReads{$w[0]} == 1){next;}
            
            
            #print "$w[9]\n";
            $nt = substr($w[9], $offset_zerobased, 1);
	    
	    $nt1 = substr($w[9], $offset_zerobased + 1, 1);
	    $nt5 = substr($w[9], $offset_zerobased + 5, 1);
	    $nt10 = substr($w[9], $offset_zerobased + 10, 1);
	    
	    $nt_1 = substr($w[9], $offset_zerobased - 1, 1);
	    $nt_5 = substr($w[9], $offset_zerobased - 5, 1);
	    $nt_10 = substr($w[9], $offset_zerobased - 10, 1);
	    
	    $score = 0;
	    
	    if($nt1 eq $b1){$score++;}
	    if($nt5 eq $b5){$score++;}
	    if($nt10 eq $b10){$score++;}
	    if($nt_1 eq $b_1){$score++;}
	    if($nt_5 eq $b_5){$score++;}
	    if($nt_10 eq $b_10){$score++;}
	    
	    if($score >= 5)
	    {
		$bases{$nt}++;
	    }
            else
	    {
		$discarded++;
		#print "1: $nt1   $b1\n";
		#print "5: $nt5   $b5\n";
		#print "-1: $nt_1   $b_1\n";
		#print "-5: $nt_5   $b_5\n";
		
	    }
        }
	
    }
}
close(FD_IN);


print FD_OUT "$opt_i\n";
print FD_OUT "$opt_c:$opt_p\n";
foreach $key (keys %bases)
{
    print FD_OUT "$key\t$bases{$key}\n";
}
print FD_OUT "frame shifts: $discarded\n";
print FD_OUT "\n---------------------\n";

close(FD_OUT);
