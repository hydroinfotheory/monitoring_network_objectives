function outpath=addbox(boxname)
% addbox(boxname,type)
% adds toolbox paths to matlabpath. 
% if 'type' is specified, can also be used to add input data paths
% uses defpaths for definition of basic paths
outpath='';
path_common='./COMMON/'
try
    outpath=genpath(fullfile(path_common, boxname))
    addpath(outpath);
end
