homedir = pwd ();
develdir = fileparts (which ("generate_control_pdf"));
pdfdir = [develdir, "/pdfdoc"];
cd (pdfdir);

function_doc

for i = 1:5
  system ("pdftex -interaction batchmode control.tex");
  system ("texindex control.fn");
endfor

cd (homedir);
