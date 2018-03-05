## bug #41820: control package: rlocus errors on input when max_k is not specified for some input
##
## https://savannah.gnu.org/bugs/?41820

## graphics_toolkit gnuplot
figure
rlocus (tf ([1, 4, 4], [1, 8, 1, 8, 0]))

this works as of feb 2018

