#!/usr/bin/perl
use strict;
use Getopt::Long;

my $gap = 5000000;
my $file = shift(@ARGV);
my $help = 0;

GetOptions("gap=i" => \$gap, "help" => \$help);

chomp $file;
my $uniqfile = $file;
$uniqfile =~ s/\.bam/.rmdup.bam/;

if($help || $file eq "")
{
	print "Usage: uniqsample.pl [options] infile\ninfile: Sorted file of reads\nOptions:\n --gap n\tCalculate unique read target counts at multiples of n [default 5000000]\n"; 
	exit 1;
}


#Samtools must be in system path
open(IN, "samtools view -h $file | ") or die "$!";   
open(UNIQ, " | samtools view -bSo $uniqfile -") or die "$!";

my ($chr, $start, $cigar, $strand, $oldline); 

my $header;
my $length = 0;
my %counts;
my %startstrandcounts;
my $lc = 0;
my $uniq = 0;
my $target;
my ($l, $ct, $op);
my @s;
my @numop;
my @cigararray;

while(<IN>)
{
	$l = $_;

	#Store header
	if($l =~ /^@/)
	{
		$header .= $l;
		print UNIQ "$l";
	}


	else
	{
		@s = split "\t", $l;

		$strand = $s[1] & 16;

		@cigararray = split(/(?<=[A-Z])/, $s[5]);

		$length = 0;

		if(!($s[1] & 1024))
		{
			print UNIQ "$l";
		}	

		if($cigararray[0] =~ /S$/)
		{
			my @nummasked = split(/(?<=[A-Z])/, $cigararray[0]);
			$s[3] -= $nummasked[0];
		}
		

		foreach $op (@cigararray)
		{
			if ($op =~ /[MDNXP]$/)
			{
				my @numop = split(/(?<=[A-Z])/, $op);
				$length += $numop[0];
			}
				
		}
		
		#if($s[3] == $start && $s[2] eq $chr)
		if($s[2] eq $chr)
		{
			$startstrandcounts{"$s[3]\t$length\t$strand"}++;
		}

		else
		{
			foreach $ct (values %startstrandcounts)
			{
				$counts{$ct}++;
			}
			
			%startstrandcounts = ();
			$chr = $s[2];
			#$start = $s[3];
		}
	}
}

foreach $ct (values %startstrandcounts)
	{
		$counts{$ct}++;
	}
	
%startstrandcounts = ();



close(IN);
close(UNIQ);

my $Rcomm = "uniroot(function(x){";

foreach my $dupct (sort {$a <=> $b} keys %counts)
{
	if($dupct)
	{
		$uniq += $counts{$dupct};
		$lc += $dupct * $counts{$dupct};
		$Rcomm .= "$counts{$dupct}*x^$dupct+";
	}
}

$Rcomm .= "-$uniq";
my @filehandles;
my @targets;
my @probs;

print "$target\t$gap\t$uniq\t".$gap * int($uniq/$gap)."\n";
for ($target = $gap * int($uniq/$gap); $target >= $gap; $target -= $gap)
{
	push @targets, $target;

	my $targetRcomm = "$Rcomm+$target}, c(0,1))\$root";

	print "$targetRcomm\n";

	my $Rresult = `Rscript -e '$targetRcomm'`;

	chomp $Rresult;

	$Rresult =~ s/^\[1\] //; 

	my $filtcount = int(0.9 + $lc * (1 - $Rresult)); # Convert fraction into number of reads to be selected

	push @probs, $Rresult;

	print "$target\t$Rresult\t$filtcount\n";

	my $targetfile = $file;
	$targetfile =~ s/\.bam$/.$target.rmdup.bam/;

	local *FILE;
	open(FILE, " | samtools view -bSo $targetfile -") or die "$!";
	#open (FILE, ">$targetfile") or die "$!";

	push(@filehandles, *FILE);
}



print "H\n";
open(IN2, "samtools view -h $file | ") or die "$!";
my $nsubs = @targets;
my @uniqcounts = (0) x $nsubs;
my $i;
my %startstranduniqline;
my $lenstrand;
my $rand;

$chr = "";
$start = 0;
$cigar = ""; 
$strand = "";
$oldline = "";

for($i = 0; $i < $nsubs; $i++)
{
	my $fh = $filehandles[$i];
	print $fh "$header";
}

while(<IN2>)
{
	#print "ere";
	$l = $_;

	if($l !~ /^@/)
	{
		@s = split "\t", $l;

		$strand = $s[1] & 16;

		@cigararray = split(/(?<=[A-Z])/, $s[5]);

		$length = 0;

		if($cigararray[0] =~ /S$/)
		{
			my @nummasked = split(/(?<=[A-Z])/, $cigararray[0]);
			$s[3] -= $nummasked[0];
		}
		

		foreach $op (@cigararray)
		{
			if ($op =~ /[MDNXP]$/)
			{
				my @numop = split(/(?<=[A-Z])/, $op);
				$length += $numop[0];
			}
		}
		
		if($s[2] eq $chr)
		{
			$startstrandcounts{"$s[3]\t$length\t$strand"}++;

			if(!($s[1] & 1024))
			{
				$startstranduniqline{"$s[3]\t$length\t$strand"} = $l;
			}
		}



		else
		{
			if($chr)
			{
			foreach $lenstrand (sort {$a <=> $b} keys %startstrandcounts)
			{
				$rand = rand();
				$i = 0;
				while($i < $nsubs && $rand > $probs[$i] ** $startstrandcounts{$lenstrand})
				{
					my $fh = $filehandles[$i];
					print $fh "$startstranduniqline{$lenstrand}";
					$uniqcounts[$i]++;
					$i++;
				}
				
			}
			}
			%startstrandcounts = ();
			%startstranduniqline = ();
			$chr = $s[2];
		}
	}
}


for($i = 0; $i < $nsubs; $i++)
{
	my $fh = $filehandles[$i];
	close $fh;
}

	

