% pack_name = "generate_html"
pack_name = "control"


% Load Packages
pkg load "generate_html"
pkg ("load", pack_name);

% Get list of functions  
list = pkg ("describe", pack_name);

%list

% Open output file
fid = fopen ("functions.texi", "w");

fprintf (fid, '@c This file is generated automatically\n');
fprintf (fid, '@c Do not edit\n\n');

for k = 1:numel (list {1}.provides)
        
        group = list {1}.provides{k};
        functions = group.functions;

        % fprintf (fid, '@section %s\n', group.category);
        fprintf (fid, '@chapter %s\n', group.category);
        
        for k=1:numel(functions)
                [TEXT, FORMAT] = get_help_text (functions(k));
                fun = functions{k};
                if (fun(1) == "@")
                  % fprintf (fid, '@subsection @%s\n', fun);
                  fprintf (fid, '@section @%s\n', fun);
                else
                  % fprintf (fid, '@subsection %s\n', fun);
                  fprintf (fid, '@section %s\n', fun);
                endif
                fprintf (fid,TEXT);
        end
        
end

fclose(fid);
