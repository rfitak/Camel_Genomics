#!/usr/bin/perl
use strict;
use warnings;

my @bases=("A","C","G","T");
foreach my $line ( <STDIN> ) {
   my $A=0;
   if ($line =~ /\tA:(\d+)/){
      $A=$1;
   }
   my $C=0;
   if ($line =~ /\tC:(\d+)/){
      $C=$1;
   }
   my $G=0;
   if ($line =~ /\tG:(\d+)/){
      $G=$1;
   }
   my $T=0;
   if ($line =~ /\tT:(\d+)/){
      $T=$1;
   }
   my @arr=($A,$C,$G,$T);
   my $index = 0;
   my $max = (sort { $b <=> $a } @arr)[0];
   for ( 0 .. $#arr){
      if ( $max eq $arr[$_] ){
         $index = $_;
      }
   }
   print "$bases[$index]\n"
}
