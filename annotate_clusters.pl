#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

local $\ = "\n";
local $, = "\t";

# Use file analysis/clustering/seqClust/clustering/RM-custom_output_tablesummary.csv
my $cluster_file = shift;

open(my $cluster_fh, $cluster_file);

my $header = <$cluster_fh>;
chomp $header;
$header = [split /\t/, $header];
splice @$header, 0, 3;

print "Cluster", "Lineage";

while(<$cluster_fh>){
     next unless s/^hits //;
     chomp;

     my ($cluster, $length, $sum, $tmp) = split /\t/, $_, 4;
     my $repeats = {};
     @$repeats{@$header} = map {s/\s//g; $_} split /\t/, $tmp;

     $cluster =~ s#clusters/dir_([^/]*)/#$1#;

     # Combine 'gypsy' and 'LTR/Gypsy' elements
     $repeats->{'LTR/Gypsy'} += $repeats->{gypsy};
     delete $repeats->{'gypsy'};

     # Remove 'Unknown' elements
     delete $repeats->{Unknown};

     # Sort elements descending by count 
     my $sort_header = [sort {$repeats->{$b} <=> $repeats->{$a}} keys %$repeats];

     # Check if first element is 'LTR'
     if($sort_header->[0] eq 'LTR'){
         # If second element is an LTR and third element is not LTR or is less
         # than half the second, add the first and second elements together and
         # remove the LTR element
         if ($sort_header->[1] =~ /^LTR/){
             if($sort_header->[2] !~ /^LTR/ || 
                $repeats->{$sort_header->[1]}/2 > $repeats->{$sort_header->[2]}){
                 
                 $repeats->{$sort_header->[1]} += $repeats->{LTR};
                 delete $repeats->{LTR};
             }
         # If second element is non-LTR and is larger than half the LTR,
         # classify cluster as Retroelement
         }elsif($sort_header->[1] eq 'non-LTR_retroposon' && 
                $repeats->{$sort_header->[0]}/2 < $repeats->{$sort_header->[1]}){
             print $cluster, 'Retroelement';
             next;
         }
     }else{
         # Remove LTR if not first element
         delete $repeats->{LTR};
     }
     # Resort elements
     $sort_header = [sort {$repeats->{$b} <=> $repeats->{$a}} keys %$repeats];

     # IF the top two elements are LTR, they're close in size, and are larger
     # than 10% of the cluster, classify cluster as 'LTR'
     if($sort_header->[0] =~ /^LTR/ && $sort_header->[1] =~ /^LTR/ && 
        $repeats->{$sort_header->[0]}/2 < $repeats->{$sort_header->[1]} &&
        ($repeats->{$sort_header->[0]}+$repeats->{$sort_header->[1]})/$sum > .1){
         print $cluster, 'LTR';
         next;
     }

     # If the first element is 'DNA' and the second element is a DNA element,
     # add the elements together and delete the DNA element
     if($sort_header->[0] eq 'DNA'){
         if ($sort_header->[1] =~ /^DNA/){
                 $repeats->{$sort_header->[1]} += $repeats->{DNA};
                 delete $repeats->{DNA};
         }
     }else{
         # Remove DNA if not first element
         delete $repeats->{DNA};
     }
     # Resort elements
     $sort_header = [sort {$repeats->{$b} <=> $repeats->{$a}} keys %$repeats];


     # If top element is larger than 10% total count and is 1.5x larger than
     # second element, classify cluster as the top element; else leave the
     # cluster unclassified (*).
     if($repeats->{$sort_header->[0]}/$sum > .1 && 
        $repeats->{$sort_header->[0]}/($repeats->{$sort_header->[1]}+1) > 1.5){
         print $cluster,  $sort_header->[0];
     }else{
         print $cluster, '*';
     }

}
close($cluster_fh);
