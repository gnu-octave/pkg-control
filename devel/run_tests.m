## Run embedded %!test blocks from the installed control package.
##
## Called by 'make check-ci'. Exit code is nonzero if any tests fail.

pkg load control;
info = pkg ("list", "control"){1};
dirs = {info.dir};
[~, nfail] = __run_test_suite__ (dirs, {});
exit (nfail > 0);
