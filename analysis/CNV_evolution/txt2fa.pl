#!/usr/bin/perl

use strict;
use warnings;

my $first = 1;
my @taxa;

my %data;


my @alphabet = qw(0 1 2 3 4 5 6 7 8 9 A B C D E F G H I J K L M N O P Q R S T U V);

my $cutoff = scalar(@alphabet);

my %count2char;
for(my $count = 0; $count < $cutoff; $count++)
{
	$count2char{$count} = $alphabet[$count];
}

while(<>)
{
	if($first)
	{
		@taxa = split /\s+/, $_;
		shift @taxa;
		$first = 0;
	}
	else
	{
		my @row = split /\s+/, $_;
		my $og = shift @row;

		$og =~ /OG(\d+)/;
		$og = $1;

		my %char;
		my $exclude = 0;
		for(my $i = 0; $i < scalar(@row); $i++)
		{
			if($row[$i] > $cutoff)
			{
				$exclude = 1;
				last;
			}
			$char{$taxa[$i]} = $count2char{$row[$i]};
		}
		if($exclude == 0)
		{
			$data{$og} = \%char;
		}
	}
}

foreach my $taxon (@taxa)
{
	print ">$taxon\n";
	foreach my $char (sort {$a <=> $b} keys %data)
	{
		print $data{$char}{$taxon};
	}
	print "\n";
}
