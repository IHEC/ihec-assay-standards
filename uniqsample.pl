#!/usr/bin/perl
use strict;
use Getopt::Long;

my $gap = 5000000;
my $outdir = "";
my $help = 0;
my $nouniq;

GetOptions("gap=i" => \$gap,
	   "outdir=s" => \$outdir,
	   "nouniq" => \$nouniq,
	   "help" => \$help);

my $file = shift(@ARGV);

if($help || $file eq "")
	{
		print "Usage: uniqsample.pl [options] infile\ninfile: Sorted bed file of reads\nOptions:\n --gap n\tCalculate unique read target counts at multiples of n [default 5000000]\n --outdir dir\tSet output directory
 --nouniq\tDo not print non-duplicates files\n --help\t\tPrint this help page\n"; 
		exit 1;
	}

chomp $file;
my $target = 1e12;   #Set initial target to more than total number of reads
my $selleft = $target;
my $totleft = $target;
my $filtcount = $target;

open(IN, "$file") or die "$!";

#ORIG file has duplicates removed
#OUT file does not have duplicates removed

my $origdupfile = $file;
if($outdir)
	{ $origdupfile =~ s/.*\//$outdir\//; }
$origdupfile =~ s/\.?[0-9]*\.bed/.dupsremoved.$target.bed/;

my $dupremfile = $origdupfile;
$dupremfile =~ s/$target.bed/bed/;

if($nouniq)
{
	open(ORIG, '>/dev/null') or die "$!";
}
else
{
	open(ORIG, ">$dupremfile") or die "$!";
}

#if($outdir)
#	{ $file =~ s/.*\//$outdir\//; }
#$file =~ s/\.?[0-9]*\.bed$/.$target.bed/;
open(OUT, '>/dev/null') or die "$!";

while($target > 0.5 * $gap)    #Count duplicates and print duplicate removed files
{
	print "Calculating for ".(($target ==  1e12) ? "all":$target)." duplicates-removed reads.\n";	

	my ($chr, $end, $start, $strand, $oldline); 

	my %counts;
	my $dupcount = 0;
	my $wc = 0;
	my $uniq = 0;

	while(<IN>)   #Count numbers of duplicates and print out file with duplicates removed
	{
		my $l = $_;

		if(int(rand()*$totleft) < $selleft) #select read with probability (filtcount - nselected)/(total[wc] - ntested)    =~ fraction of reads to be kept
		{	
			if($target != 1e12)
			{	print OUT "$l";	}

			chomp $l;
			my @s = split "\t", $l;	

			$wc++;
		
			if($s[2] == $end && $s[1] == $start && $s[5] eq $strand && $s[0] eq $chr)
			{
				$dupcount++;
			}

			else
			{
				if($dupcount) 
				{
					$counts{$dupcount}++;    # $counts{n} measures the number of reads with exactly n duplicates (>=1)
					if(!$nouniq) {print ORIG "$oldline\n";}
				}
				$dupcount = 1;
			}	

			($chr, $end, $start, $strand) = ($s[0], $s[2], $s[1], $s[5]);

			$oldline = $l;
			$selleft--;
		}
		$totleft--;
	}

	$counts{$dupcount}++;
	if(!$nouniq) {print ORIG "$oldline\n";}
	close(IN);
	close(ORIG);
	close(OUT);


	#Create R command for calculating total number of reads 
	# $uniq = current number of unique reads
	# $target = desired number of unique reads
	# Solve $target - $uniq + sum(n * $counts{n}) = 0 


	my $Rcomm = "uniroot(function(x){";

	foreach my $ct (sort {$a <=> $b} keys %counts)
	{
		$uniq += $counts{$ct};
		$Rcomm .= "$counts{$ct}*x^$ct+";
	}

	print "$uniq duplicates removed reads.\n"; #Print out number of unique reads from previous iteration
	
	if($target == 1e12) #Rename target file with true number of unique reads 
	{
		#my $origrename = $origdupfile;
		#$origrename =~ s/dupsremoved\.$target/$uniq/;
		#`mv $origdupfile $origrename`;
		$target = $gap * int($uniq/$gap);
	}
	else			#Set input file to filtered set using previous target number file including duplicated reads
	{
		my $uniqfile = $file;
		$uniqfile =~ s/$target/$uniq/;	
		my $origdupfilenew = $origdupfile;
		$origdupfilenew =~ s/$target/$uniq/;
		`mv $file $uniqfile`;
		if(!$nouniq)
			{`mv $origdupfile $origdupfilenew`;}
		$file = $uniqfile;

		if($target < 1.5 * $gap)
		{
			last; #Final filtering for approximately $gap unique reads
		}
		else
		{
			$target -= $gap;
		}
	}		

	
	$Rcomm .= "$target-$uniq}, c(0,1))\$root";

	my $Rresult = `Rscript -e '$Rcomm'`; #Run R command and parse result

	chomp $Rresult;
	$Rresult =~ s/^\[1\] //;

	$filtcount = int(0.9 + $wc * (1 - $Rresult)); # Convert fraction into number of reads to be selected

	#Use previous OUT file as new input file and set temporary filenames for new OUT and ORIG files

	open(IN, "$file") or die "$!";

	if($outdir)
		{ $file =~ s/.*\//$outdir\//; }
	$file =~ s/\.?[0-9]*\.bed$/.$target.bed/;
	open(OUT, ">$file") or die "$!";
	

	$origdupfile = $file;
	$origdupfile =~ s/\.?[0-9]*\.bed/.dupsremoved.$target.bed/;
	#open(ORIG, ">$origdupfile") or die "$!";
	if($nouniq)
	{
		open(ORIG, '>/dev/null') or die "$!";
	}
	else
	{
		open(ORIG, ">$origdupfile") or die "$!";
	}

	$selleft = $filtcount;
	$totleft = $wc;
}


