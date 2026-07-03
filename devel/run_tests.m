## Run embedded %!test blocks from source without installing the package.
##
## The canonical way (according to the docs I've found) to to actually test is
## to call 'pkg test control' but that  is coupled to docs today because it
## requires you to install first, but to install you have to build a release
## tarball with dist which depends on docs. This is a "hack" to get around that
## and just run the tests directly. This should probably change in the future
## b/c searching and grepping through a bunch of src files to test is very
## awkward.
##
## Called by 'make check-ci'. Exit code is nonzero if any tests fail.

## Register autoloads from PKG_ADD directives (normally done by pkg install)

src_path = fullfile (fileparts (mfilename ("fullpath")), "..", "src");

## Find all .cc files in src/ b/c those are compiled to .oct files. We literally
## just want this out of the cpp files so that we know what .oct file we should
## target
## // PKG_ADD: autoload ("__sl_sb03md__", "__control_slicot_functions__.oct");
cc_files = glob (fullfile (src_path, "*.cc"));

## Read every file and grep of PKG_ADD, call autoload() on it with the absolute
## path to the .oct file. Autoload registers a deffered load which tells octave
## "when someone calls funcname, load filepath to find it, the .oct file isn't
## actually loaded until that first call"
for i = 1:numel (cc_files)
  txt = fileread (cc_files{i});
  ## Use a regex to capture the 2 arguments (function name and .oct filename)
  toks = regexp (txt, 'PKG_ADD: autoload \("(\w+)", "(\w+\.oct)"\)', "tokens");
  for j = 1:numel (toks)
    autoload (toks{j}{1}, fullfile (src_path, toks{j}{2}));
  endfor
endfor

## Find all the class dirs to test ("inst", "inst/@tf", etc...)
dirs = [{"inst"}, cellstr(glob ("inst/@*"))'];

n_pass = 0;
n_total = 0;

## For every dir we collected, go through every mfile and call test
for i = 1:numel (dirs)
  mfiles = dir (fullfile (dirs{i}, "*.m"));
  for j = 1:numel (mfiles)
    f = fullfile (mfiles(j).folder, mfiles(j).name);
    ## collect output args from test call
    ## https://docs.octave.org/v11.1.0/Test-Functions.html
    [p, total, xfail, xbug, skip, rtskip] = test (f, "quiet", stdout);
    n_pass += p;
    n_total += total - xfail - xbug - skip - rtskip;
  endfor
endfor

n_fail = n_total - n_pass;
fprintf ("\nTotal: %d passed, %d failed\n", n_pass, n_fail);
exit (n_fail > 0);
