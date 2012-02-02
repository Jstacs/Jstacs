open(IN,$ARGV[0]);

my @lines = <IN>;

open(OUT,">".$ARGV[0]);

for my $line (@lines){
	if($line =~ /^\s*package\s+supplementary\.cookbook\.recipes/){
		print OUT "//",$line;
	}else{
		print OUT $line;
	}
}

close(OUT);