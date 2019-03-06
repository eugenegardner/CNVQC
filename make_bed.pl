#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $path = $ARGV[0];

open (FILES, "-|", "ls $path/split*") || die "Cannot open file: $!";

print "#!/bin/bash\n";

foreach (<FILES>) {

	chomp $_;
	
	my ($name, $dir, $suffix) = fileparse($_);
	
	print "tail -n +2 $_ | perl -ane \'chomp \$_; print \"\$F\[1\]\\t\$F\[2\]\\t\$F\[2\]\\t\$F\[0\]\\t\$F\[3\]\\t\$F\[4\]\\n\";\' | bedtools sort -i stdin > /home/wpcgk/Biobank_analysis/split_files/$name.sorted.bed\n";

	print "bgzip /home/wpcgk/Biobank_analysis/split_files/$name.sorted.bed\n";

	print "tabix -p bed /home/wpcgk/Biobank_analysis/split_files/$name.sorted.bed.gz\n";

}
