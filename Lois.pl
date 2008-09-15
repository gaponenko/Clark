#!/usr/bin/perl

# Simple script summing all the root histograms from a mofia analysis.
#
# Usages:
# Lois.pl /Path/to/goodlinks Output.root
# Lois.pl setXX aYY
#
# Example:
# Lois.pl /twist/tw00y/systematics/data/set68/anal1/goodlinks set68a1_lois.root
#
# Lois.pl set68 a1
# 	It will automatically check the goodlinks directory on tw00y and create
# 	the file set68a1_lois.root


use strict;

my $InputDir;
my $Output;
my $Anum;

if ($ARGV[1] =~ /^a([0-9_]+)/ )
{
	$Anum	= $1;
}

if ( ($ARGV[0] =~ /^set/ || $ARGV[0] =~ /^gen/ ) && $ARGV[1] =~ /^a\d+/ )
{
	if ($ARGV[0] =~ /^set/ )
	{
		$InputDir	= "/twist/tw00y/systematics/data/$ARGV[0]/anal$Anum/goodlinks";
	}
	else
	{
		$InputDir	= "/twist/tw00y/systematics/mc/$ARGV[0]/anal$Anum/goodlinks";
	}
	$Output		= "$ARGV[0]$ARGV[1]_lois.root";
}
else
{
	$InputDir	= $ARGV[0];
	$Output		= $ARGV[1];
}
	
print "Analysing trees from $InputDir\n";

my @List		= `ls $InputDir/*/0*.root`;

my $ComList		= "";

foreach my $L (@List)
{
	chomp($L);
	$ComList	.= " ".$L;
}


system("/home/e614/analysis_2006_2007/treesum/common/ptta/histoAdd $Output $ComList");
# print("/home/e614/analysis_2006_2007/treesum/common/ptta/histoAdd $Output $ComList\n");

