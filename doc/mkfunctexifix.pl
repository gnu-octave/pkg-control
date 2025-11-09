#!/usr/bin/env perl

## Copyright 2018-2024 John Donoghue
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## File    : mkfunctexifix.pl
## Purpose : Replace @sealso command and conditionals in the functions.texi
##           file created by myfuncdocs.py from all texifo texts for
##             1. being able to use full latex commands (@tex -> @latex)
##             2. having correct refs from the @seealso command
## Usage   : myfunctexifix.pl functions.texi 


use warnings;              # report warnings for questionable run-time code
use strict qw(refs subs);  # check at compile-time for bad programming style
use File::Basename;        # For splitting paths between dir and file
use File::Spec;            # For combining dirs and filenames into paths

sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

# Extract command line arguments
if ($#ARGV != 0) { die "USAGE: $0 functions.texi;" }

my $fname = basename ($ARGV[0]);

open TEXI, "<".$fname    or die "Unable to open $fname";
my $texi = do {local $/; <TEXI> };
close TEXI;

$texi =~ s/\@tex/\@latex/g;
$texi =~ s/\@end tex/\@end latex/g;
$texi =~ s/\@iftex/\@iflatex/g;
$texi =~ s/\@end iftex/\@end iflatex/g;
$texi =~ s/\@ifnottex/\@ifnotlatex/g;
$texi =~ s/\@end ifnottex/\@end ifnotlatex/g;

while ($texi =~ m/\@seealso/) {

$texi =~ s/\@seealso[\s]*\{([^\}]*)\}//;

  my $rtext = "\@strong\{See also:\}";
  my $ltext = "\@strong\{See also:\}";
  my @seealso = split (/\,/,$1);
  foreach my $ref (@seealso) {
    $ref =~ s/^\s+|\s+$//g;
    $ltext = $ltext." \@link\{$ref,$ref\},";
    $rtext = $rtext." \@ref\{$ref\},";
  }
  $rtext = substr $rtext, 0, -1;
  $ltext = substr $ltext, 0, -1;
  my $text = "";
  $text = $text."\@iflatex\n$ltext.\n\@end iflatex\n";
  $text = $text."\@ifnotlatex\n$rtext.\n\@end ifnotlatex\n";
  $texi = $`.$text.$';

} 

open TEXI, ">".$fname or die "Unable to open $fname";
print TEXI $texi;
close TEXI;

