% draft code
% needed for LTI array support in display routines
% sys(:, :, idx(1), idx(2), idx(3), idx(4), ...)

% nc = numel (sys.a)  % sys.a{...} cell array
% sc = size (sys.a)
% k = num2cell (reshape (1:nc, sc))

% cellfun (@display_scalar, k, {sc}, sys.a, sys.b, sys.c, sys.d, ..., "uniformoutput", false)

c = cell (3, 4, 2, 1, 5)

nc = numel (c)
sc = size (c)

csc = cumprod (sc)


div = [1, csc(1:end-1)]

% k = 115
%k = 120
k = 119 % mat (2,4,2,1,5) = mat (119)




mat = reshape (1:nc, sc)
arrayfun (@disp, mat, "uniformoutput", false);



idx = 1 + rem (fix ((k-1)./div), sc)

