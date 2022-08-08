#!/usr/bin/perl
# splitFASTA

# Shawn Polson - 20 Jan 2007; edit: 5 Feb 2007, 1 Mar 2009
# Accepts multi-sequence FASTA file and returns split multi-sequence file chunks of a specified size.
# Args: none (user prompted for input)
# Dependencies: none
# splitFASTA.pl infasta outputDir outputName numberofseqs
use strict;
use Getopt::Long;
use File::Basename;
use POSIX;

my $header ="
splitFASTA.pl

 VERSION: 1.4

 SYNTAX:  
       splitFASTA.pl -i inputfile -o outputDirectory -s #sequences
   OR
       splitFASTA.pl -i inputfile -o outputDirectory -f #files


 ARGUMENTS:
       -i input file (fasta format) REQUIRED
       -o output directory REQUIRED
       -s sequences per file (default 1000)
       -f number of files (overrides '-s')
       -h display help
";

my ($inputFile, $outPath, $help, $files);
my $count=1000;

GetOptions (    		"i|input=s"     =>              \$inputFile,
                                "o|output=s"    =>              \$outPath,
                                "s|sequences=s"   =>             \$count,
				"f|files=s"	=>		\$files,
                                "h|help"                =>       \$help);



# BASIC OPS
die $header if($help);
open (DAT, $inputFile) or die "Error! Cannot open file $inputFile\n";
if(! -d "$outPath") 
{	`mkdir -p $outPath`; 
	if (! -d "$outPath") { die "Error!  Could not create output directory $outPath!\n"; }
}
my $outFile=fileparse($inputFile, qr/\.[^.]*$/);
$inputFile =~ /\.([^.]*)$/;
my $suffix = $1;

#GLOBALS
my $j=0;
my $fileNumber=1;
  

# CHECK IF NUMBER OF FILES METHOD
if ($files)
{	my $i = `egrep -c '^>' $inputFile`;
	chomp $i;
	$count = ceil($i/$files);
}

open (DAT, $inputFile) or die "Cannot open file $inputFile\n";

open (OUT, "> $outPath/$outFile-$fileNumber.$suffix") or die "Error! Cannot create output file: $outPath/$outFile-$fileNumber.txt\n";

while(<DAT>)
{	if ($_ =~ /^>/) #if new sequence increase counter
	{	$j++;
	}
	if ($j > $count) #if time for new output file
	{	close(OUT);
		$fileNumber++;
		open (OUT, "> $outPath/$outFile-$fileNumber.$suffix") or die "Error! Cannot create output file: $outPath/$outFile-$fileNumber.txt\n";
		$j = 1;
	}
	
	#Output line to file
	print OUT $_;
}

close(OUT);
if ($j == 0)
{	my $sys = `rm $outPath/$outFile-$fileNumber.txt`;
	$j = $count;
}

print "File splitting completed.  $fileNumber files created.  Final file contains $j sequences.\n";
exit 0;
