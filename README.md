# Octave control package

This is the official repository for the control package for GNU Octave.

## About

The **control** package is a collection of functions for control systems design and analysis. 

As of 24.03.2023, the developemnt of the **control** package was moved from [SourceForge](https://sourceforge.net/p/octave/control/ci/default/tree/) and [Mercurial](https://en.wikipedia.org/wiki/Mercurial) to [GitHub](https://github.com/gnu-octave/pkg-control) and [Git](https://en.wikipedia.org/wiki/Git). Links related to the control package

- [License and copyright information](https://github.com/gnu-octave/pkg-control/blob/main/COPYING)
- [Releases](https://github.com/gnu-octave/pkg-control/releases)
- [Documentation](https://gnu-octave.github.io/pkg-control)

## Used Library SLICOT

Control uses some routines of the [SLICOT-Reference library](https://github.com/SLICOT/SLICOT-Reference) (Copyright (c)  1996-2025, The SLICOT Team). The sources of the used routines are included in the released control package archive `control-x.y.z.tar.gz` in the directory `src/slicot-src` and are compiled for the target system while installing the control package for Octave.

The SLICOT files are available under the *BSD 3-Clause License* which can be found

- in the file `src/slicot-src/LICENSE` (together with README files) in the package archive `control-x.y.z.tar.gz`,
- in the file `doc/SLICOT/LICENSE` (together with README files) in the package installation directory, or
- in the [SLICOT-Reference repository](https://github.com/SLICOT/SLICOT-Reference/blob/main/LICENSE).

Reference:

- Köhler, M., Saak, J., Sima, V., & Varga, A. (2025). SLICOT - Subroutine Library In COntrol Theory (Version 5.9.1) [Computer software]. DOI: [10.5281/zenodo.17523371](https://doi.org/10.5281/zenodo.17523371)

## Installing the control package

### Installing released package version

The easiest way to install the newest control package is to type

  `pkg install control`

For installing a certain version x.y.z of the control package, you may

- download the package archive file `control-x.y.z.tar.gz` of one of the [releases](https://github.com/gnu-octave/pkg-control/releases) and install it by typing
  `pkg install control-x.y.z.tar.gz` or
- directly issue the command `pkg install "https://github.com/gnu-octave/pkg-control/releases/download/control-x.y.z/control-x.y.z.tar.gz"`

### Creating and installing package archives from the sources

You can also clone this repository (using the option `--recurse-submodules` since SLICOT is included as git submodule) and build the package archive file by yourself. For this, you can use the following commands:

- `make dist`<br>
  Create the package archive file in the directory `target` which can be installed in Octave afterwards
- `make install`<br>
  Install the package
- `make help`<br>
  Show all targets for `make`

## Contributing to the control package

Information on how to contribute to the control package can be found in [this document](CONTRIBUTING.md).
