# How to Contribute

## Ways to contribute

There are many ways to help improving the **control** package for GNU Octave.

### Improving the documentation

If you find a function or feature with missing, erroneous, incomplete or incomprehensible help text, you can help us by contributing changes to address the issues. For writing help texts using the [Texinfo](https://www.gnu.org/software/texinfo/) format, please refer to [License and Documentation](#license-and-documentation).

### Reporting and/or fixing bugs

Bugs and feature requests for GNU Octave including the packages are tracked on [Github](https://github.com/gnu-octave/pkg-control/issues). Please feel free to start working on a fix or to help by testing. Patches for fixing the bug can be attached to a comment to the issue or - preferably - by opening a pull request in the package's repository (see [Contribution Workflow](#contribution-workflow)). If you would like to report a new bug, please go over existing bugs before in order to avoid duplicate reports.

### Discussions

If you have any questions or suggestions on how to extend or improve the package, please feel free to participate the community on 
[GNU Octave Discourse](https://octave.discourse.group/).

## License and Documentation

The **control** package is distributed under [GNU General Public License (GPL)](https://www.gnu.org/licenses/gpl-3.0.en.html) (except for the used SLICOT files). If you are submitting a few function, it should be licensed under GPLv3+ with the following header (use appropriate year, name, etc.) as shown below:

```bash
## Copyright (C) YEAR NAME <E-MAIL>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

```

New functions or features should be properly documented with the embedded help file in [Texinfo](https://www.gnu.org/software/texinfo/) format. This part should be placed outside (before) the function's main body and after the License block. The Texinfo manual can be found [here](https://www.gnu.org/software/texinfo/manual/texinfo/) although reading through existing function files should do the trick.

```bash
## -*- texinfo -*-
## @deftypefn {Function File} {} pzmap (@var{sys})
##
## Help file goes here.
##
## @end deftypefn
```

The texinfo is not printed on the command window screen as it appears in the source file.


## Coding style

### General rules

The general coding style for GNU Octave given in the [GNU Octave Wiki](https://wiki.octave.org/Octave_style_guide) should be used with the following additions:

- Limit the line length to 80 characters
- Use `LF` (unix) for end of lines, and NOT `CRLF` (windows)


### Tests

It is very helpful that function files contain tests for correct output. The tests are located at the end of the file with lines beginning with `%!`. As example, please finde a test for `pzmap ()` below.

```
%!test
%! s = tf('s');
%! g = (s-1)/(s-2)/(s-3);
%! [pol zer] = pzmap(g);
%! assert(sort(pol), [2 3]', 2*eps);
%! assert(zer, 1, eps);
```

### Demos

Although examples of using a function should already be provided in the documentation, it is always useful to have examples embedded as demos, which the user can invoke with the `demo` command.

```
>> demo pzmap
```

Like test, demos are also located at the end of the file. A small demo for `pzmap ()` is shown below.

```
%!demo
%! s = tf('s');
%! g = (s-1)/(s-2)/(s-3);
%! pzmap(g);
```

## Contribution Workflow

As in many other open-source projects the usual way to contribute to the control package is to

- fork the repository,
- make changes, like, e.g., fix a bug or add a new feature in your fork, and
- send pull requests to the original repository.

Please also refer to this [detailed description on collaborating with pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests).

Your contribution must be an independent work or derived from code that may be released under the terms of the GPL. **Under no circumstances may it be based on code from Matlab or other non-free code that you may have access to view**.

### Fork and build

[Fork](https://github.com/gnu-octave/pkg-control/fork) the pkg-control repository to your own account and clone the resulting repository. Refer to the [README](README.md) file for information about how to build the package archive which can be installed in GNU Octave.

### Pull request

When your changes are finished, commit and push the change to your forked repository on Github (make sure your fork is up to date) and create a pull request.

### Option for very small changes

If the changes are small and only affect one file, you can make the changes directly in the web interface of Github, select *Commit changes* and *Create a new branch for this commit and start a pull request*.
 
 
