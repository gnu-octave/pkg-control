# Release Preparation for the Control Package

## Set version number etc.

1. Merge branch `default` into branch `release_preparation`
2. Update `NEWS`, possibly move old entries into `ONEWS`
3. Update version in `DESCRIPTION`
4. Update version in `src/configure.ac`
5. `hg commit` for storing above changes

## Testing

6. `make distclean` for removing all old distribution files
4. Tests
    1. `make dist`
    2. `pkg install target/control-X.Y.Z.tar.gz` in Octave
        - any warnings or even errors?
        - correct version number?
    3. `pkg test control` in Octave
    4. `pkg uninstall control` in Octave
    5. `make install`
        - any warnings or even errors?
        - correct version number?
    6. `pkg test control` in Octave

## Release the package

5. `make release`<br>
   *Alternatively*, this can be done by the following separate commands
    1. `make dist`
    4. `md5sum target/control-X.Y.Z.tar.gz`
    1. `make html`
    4. `md5sum target/control-html.tar.gz`
6. Create a [package release ticket on Sourceforge](https://sourceforge.net/p/octave/package-releases/):
    - Upload `control-X.Y.Z.tar.gz` and `control-html.tar.gz`
    - Provide the related commit in the repository
    - Provide the related md5sums
7. `hg tag 'control-X.Y.Z'`
8. Merge `release_preparation` into stable
9. Merge `stable` into `default`
