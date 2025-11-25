
n = 1:100;

sys = arrayfun (@rss, n, "uniformoutput", false);

st = cellfun (@isstable, sys);

all_stable = all (st)
