# Release Preparation for the Control Package

## Set version number etc.

1. Possibly merge branch `dev` into branch `main`
2. Update `NEWS`, possibly move old entries into `ONEWS`
3. Update version in `DESCRIPTION`
4. Update version in `src/configure.ac`
5. `git commit` for storing above changes

## Testing

6. `make distclean` for removing all old distribution files
4. Tests
    1. `make dist`
    2. `pkg install target/control-X.Y.Z.tar.gz` in Octave
        - any warnings or even errors?
        - correct version number?
    3. `pkg test control` in Octave
    4. `pkg uninstall control` in Octave
    5. `make release` (including docs)
        - any warnings or even errors?
        - correct version number?
    6. `pkg test control` in Octave

## Release the package

5. `make release`<br>
   *Alternatively*, this can be done by the following separate commands
    1. `make dist`
    4. `sha256sum target/control-X.Y.Z.tar.gz`
    1. `make install`
    2. `make docs`
3. `git commit` for the docs
4. `git tag control-X.Y.Z`
5. `git push`, docs should be automatically be published
5. Make release in Github
    - Select tag
    - Upload `control-X.Y.Z.tar.gz`
6. Update `control.yaml` in `packages` repository and make a pull request
9. Merge `main` into `dev`
