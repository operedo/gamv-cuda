$nx=420;
$ny=600;
$nz=400;

open(IN,"<$ARGV[0]");
#open(TOP,"<$ARGV[0]");
#
#%top=();
#<TOP>;
#for($j=0;$j<$ny;$j++){
#for($i=0;$i<$nx;$i++){
#	$line=<TOP>;
#	chomp($line);
#	$top{$i.",".$j}=$line;
#}
#}



#print "Test\n";
#print "4\n";
#print "X\n";
#print "Y\n";
#print "Z\n";
#print "Var\n";
for($k=0;$k<$nz;$k++){
for($j=0;$j<$ny;$j++){
for($i=0;$i<$nx;$i++){
	$line=<IN>;
	chomp($line);
	#print 25*$i,",",25*$j,",",$top{$i.",".$j}-10*$k,",\n";
	print 1*$i," ",1*$j," ",1*$k," ",$line,"\n";
}
}
}
