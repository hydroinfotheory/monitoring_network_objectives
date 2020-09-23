function svg2pdf(filename)
%%
%filename=[pwd 'single_vennaa.svg']
% svg2pdf(filename)
%uses inkscape to convert svg file to pdf format
inkscape_path='C:\Program Files\Inkscape'
olddir=pwd;
cd(inkscape_path);
outfilename=[filename(1:end-4) '.pdf'];
doscommand=['inkscape.exe ' filename ' --export-pdf=' outfilename]

dos(doscommand);
cd(olddir);