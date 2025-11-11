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
## Purpose : Replace @sealso command and \( ... \) in the file created
##           by myfuncdocs.py from all texifo texts for
##             1. prevent errors in makeinfo due to \(...\)
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

# makeinfo does not support \(...\), only $...$, which in turn
# is not supported by pkg-octave-doc. Solution: use \(...\)
# in help texts and replace them here.
$texi =~ s/\\\)/\$/g;
$texi =~ s/\\\(/\$/g;

# Now replace all @seealso by individual @link or @ref commands
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
  $text = $text."\@iftex\n$ltext.\n\@end iftex\n";
  $text = $text."\@ifnottex\n$rtext.\n\@end ifnottex\n";
  $texi = $`.$text.$';

} 

open TEXI, ">".$fname or die "Unable to open $fname";
print TEXI $texi;
close TEXI;

