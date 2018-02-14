open(IN0,"<$ARGV[0]");
open(IN1,"<$ARGV[1]");

<IN0>;
<IN1>;

$counter=0;
$meanh=0.0;
$stdevh=0.0;
$meangam=0.0;
$stdevgam=0.0;

while(<IN0>){
	$_=~/([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+(.*)/;
	$id0=$1;
	$h0=$2;
	$gam0=$3;
	$head0=$4;
	$tail0=$5;
	chomp($tail0);

	$line=<IN1>;
	$line=~/([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+(.*)/;
	$id1=$1;
	$h1=$2;
	$gam1=$3;
	$head1=$4;
	$tail1=$5;
	chomp($tail1);
		
#	print "$h0,$gam0,$h1,$gam1\n";

	$counter=$counter+1;
	$meanh=$meanh+abs($h0-$h1);
	$meangam=$meangam+abs($gam0-$gam1);
	$stdevh=$stdevh+($h0-$h1)*($h0-$h1);
	$stdevgam=$stdevgam+($gam0-$gam1)*($gam0-$gam1);
}

$v1=$meanh/$counter;
$v2=$stdevh/$counter - $v1*$v1; 
$v3=$meangam/$counter; 
$v4=$stdevgam/$counter - $v3*$v3; 

print "$v1,$v2,$v3,$v4\n";
#print ($meanh/$counter),",",($stdevh/$counter - ($meanh/$counter)*($meanh/$counter)),",",($meangam/$counter),",",($stdevgam/$counter - ($meangam/$counter)*($meangam/$counter)),"\n";

