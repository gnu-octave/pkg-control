## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2011
## Version: 0.1

function dat = iddata (y = [], u = [], tsam = [], varargin)

  if (nargin == 1 && isa (y, "iddata"))
    dat = y;
    return;
  elseif (nargin < 3)
    print_usage ();
  endif

  if (! issample (tsam, 1))
    error ("iddata: invalid sampling time");
  endif

  [p, m] = __iddata_dim__ (y, u);

  outname = repmat ({""}, p, 1);
  inname = repmat ({""}, m, 1);

  dat = struct ("y", y, "outname", {outname}, "outunit", {outname},
                "u", u, "inname", {inname}, "inunit", {inname},
                "tsam", tsam, "timeunit", {""},
                "name", "", "notes", {{}}, "userdata", []);

  dat = class (dat, "iddata");
  
  if (nargin > 3)
    dat = set (dat, varargin{:});
  endif

endfunction