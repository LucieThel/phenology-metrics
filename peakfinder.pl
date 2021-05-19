#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Text::CSV;

# Reference : Zerbe P., Clauss M., Codron D., Bingaman Lackey L., Rensch E., Streich J. W., Hatt J.-M. & Müller D. W. (2012). Reproductive seasonality in captive wild ruminants: implications for biogeographical adaptation, photoperiodic control, and life history. Biological Reviews, 87, 965-990. DOI: 10.1111/j.1469-185x.2012.00238.x

my $help=0;
my $target=80;
my $column=4;
my $numberOfRecs=75;
my $filename;

my @rows=();
my @cols;
my $sum=0;

my $trigger;
my $triggered;
my $windowSize;
my @windowSums;
my $r;
my $inner;

my$maxIndex;
my $maxval;



GetOptions (
		"t=i" => \$target,			# trigger percentage
		"c=i" => \$column, 			# number of the column containing the data
		"n=i" => \$numberOfRecs, 	        # number of records to analyse
		"help" => \$help
		) or pod2usage(1);

pod2usage(1) if $help;
if ($#ARGV < 0) {
	pod2usage(1)
}

my $csv = Text::CSV->new({ sep_char => ';' });		# define field separator

$filename = $ARGV[0];
open IN, $filename or die "Can't open file  $filename\n";
print "\nAnalyzing $numberOfRecs records in $filename\n";
while(my $line = <IN>) {
	chomp($line);
	if ($csv->parse($line)) {
		@cols = $csv->fields();
	} else {
		warn "Line could not be parsed: $line\n";
	}
  	if ($cols[$column]=~ /([\d\.]+)/) {
		$rows[$#rows+1] = $1;
		$sum+=int($1);
		last if($#rows+1==$numberOfRecs);
	} else {
		warn "wrong format in col $column : cols[$column] \n";
	}
}
close IN;

$trigger = $sum*$target/100.;
print $#rows+1," records, total= ",int($sum), " animals, trigger at ", $target,"% = ",$trigger,"\n";

$triggered=0;

if ($numberOfRecs > $#rows+1) {
	warn "Number of records (".$numberOfRecs.") to analyse is bigger  than number of records in data (". ($#rows+1).")";
}


for($windowSize=1; $windowSize<=$numberOfRecs; $windowSize++) {
	@windowSums = ();
	for($r=0;$r<$numberOfRecs-$windowSize;$r++) {
		$windowSums[$r]=0;
	 	$inner=$r;
	 	while ($inner <$r+$windowSize ) {
			if ($inner > $#rows) {
				warn "Error : no value for the row at $inner"
			}
			$windowSums[$r] += $rows[$inner];
			$inner ++;
		}
     	if($windowSums[$r] >= $trigger) {
			$triggered++;
		}
	}
	last if($triggered>0);
}

if($triggered>0) {
	print "\nWindow found at size of ", $windowSize, " data blocks\n";
	$maxIndex=-1;
	if($triggered>1) {
    	$maxval=0;
    	for($r=0;$r<$numberOfRecs-$windowSize;$r++) {
			if($windowSums[$r]>=$trigger && $maxval<$windowSums[$r]) {
				$maxval=$windowSums[$r];
				$maxIndex = $r;
			}
		}
	}
  	for($r=0;$r<$numberOfRecs-$windowSize;$r++) {
		if($windowSums[$r]>=$trigger) {
			printf "%s %d births (prop of %1.1f%%) from data block %d to %d.\n",($r==$maxIndex)?"**":"--", $windowSums[$r], $windowSums[$r]/$sum*100, $r, $r+$windowSize-1;
     	}
  	}
}
else {
  print "\nNo window found\n";
}


__END__

=head1 SYNOPSIS

perl peakfinder_2.pl [--t 80 ] [--c 5 ] [--n 75] input_file

=head1 OPTIONS

=over 8

=item B<--help>
Print a brief help message and exits.

=item B<--t>
    set  trigger percentage (integer, default: 80)

=item B<--c>
    set  number of the column containing the data (integer, default: 5)

=item B<--n>
	 set   number of records (integer, default: 75)

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut
