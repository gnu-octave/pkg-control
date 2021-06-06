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
        firstslash = 0;

        for k=1:numel(functions)
                fun = functions{k};
                [TEXT, FORMAT] = get_help_text (fun);
                if (fun(1) == "@")
                  fprintf (fid, '@section @%s\n', fun);
                  firstslash = (find (fun == '/'))(1);
                  if ! firstslash
                    error (["function_doc: unknown class ", fun]);
                  end
                  fprintf (fid, '@findex %s\n', fun(firstslash+1:end));
                else
                  fprintf (fid, '@section %s\n', fun);
                  fprintf (fid, '@findex %s\n', fun);
                endif
                fwrite (fid,TEXT);
        end
        
end

fclose(fid);
