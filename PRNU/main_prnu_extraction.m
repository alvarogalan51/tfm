th% script - run_prnu_extraction.m
%
% Script to run prnu_extraction.m To run the script the user must specify the following:
%
%   pmode:      'averagePRNU': the PRNU pattern is extracted from multiple noise
%                   residuals
%               'singlePRNU': a single PRNU pattern (in fact a noise residual)
%                   is saved for each image
%   inPath:      Folder name containing input images (flat-field or natural)
%   cameraID:    An identifier at the start of the PRNU filename relating to
%                   the camera make and model
%   imageType:   Set this to NI for natural images and FF for full frame
%   denoiseFilter: Select which denoising filter you wish to use.
%   iList:       Specify individual filenames or leave undefined to use all image
%                   files in folder.
%
% The format of the PRNU file is:
%       <cameraID>_<imageType>_<Image#>_denoiseFilter>_<alg>
%
% or for PRNU pattern averaged over multiple files
%
%       <cameraID>_<imageType>_AVE_denoiseFilter>_<alg>
%
% Requirements:
% - The Matlab Image Processing Toolbox is required.
% - The Matlab Wavelets Toolbox if the Mihcak denoising filter is used.
% - The BM2D package is needed f denoising using this method will be used.
%   The package can be downloaded from http://www.cs.tut.fi/~foi/GCF-BM3D/.
%
% H Muammar (hkmuammar@gmail.com)
% 25 January 2012
%
% THE SOFTWARE FURNISHED UNDER THIS AGREEMENT IS PROVIDED ON AN 'AS IS' BASIS, 
% WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS OR IMPLIED, INCLUDING, 
% BUT NOT LIMITED TO, ANY IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS
% FOR A PARTICULAR PURPOSE.

clear

featureStr = {'MATLAB'; 'Image_Toolbox'; 'Wavelet_Toolbox'};


index = cellfun(@(f) license('test',f),featureStr);
availableFeatures = featureStr(logical(index));

% Check that the Image Processing Toolbox is installed
features = strcmp('Image_Toolbox', availableFeatures);

if ~sum(features)
    fprintf('Image Processing Toolbox is missing.\n');
    return
end

% Check that Wavelet Toolbox is installed
features = strcmp('Wavelet_Toolbox', availableFeatures);
if ~sum(index)
    wavInst = true;
else
    wavInst = false;
end

% Set Discrete Wavelet Transform extension mode to
% periodization.
if wavInst
    if ~strcmp(dwtmode('status', 'nodisp'), 'per')
        dwtmode('per', 'nodisp');
    end
end

opMode = 'reference';
%opMode = 'correlate';

fprintf('run_prnu_extraction: Mode is %s\n', opMode);

switch opMode

% ++++++++++++++++++++++++++++++ Extract PRNU +++++++++++++++++++++++++++++

case 'reference'

% ------- User configurable parameters for 'reference' mode ----------

pmode = 'singlePRNU';      % 'singlePRNU' : compute an individual PRNU pattern (noise residual) for each image
%pmode = 'averagePRNU';     % 'averagePRNU': compute the PRNU by averaging over all noise residuals

inPath = 'D:\Dropbox\Máster en Ciberseguridad\TFM\Definitivo\Fotos\motog5';      % folder containing input images

outPath = 'D:\Dropbox\Máster en Ciberseguridad\TFM\Definitivo\Fotos\motog5\out';      % This is where the prnu reference data are saved (.mat file)

cameraID = 'Mobile';  % Identifier for camera make/model

imageType = 'FF';         % Specify argument as FF for flat-field and NI for natural image

denoiseFilter = 'sigma'; % Specify the denoising filter: 'mihcak', 'sigma', 'gaussian','bm3d'

saveWaveletCoeffs = 0;      % (0) Don't save wavelet coefficients (1) Save the wavelet coefficients. NOTE - only applies when using the Mihcak wavelet based denoising algorithm. Coefficient is ignored otherwise

overwrite = 0;          % (0) use existing denoised images; (1) Overwrite existing denoised images

% Specify image names. Two options are available:
%    1) If iList variable is commented out or not defined, then all the
%    images in the folder inPath are used.
%    2) If iList is defined then only the images listed in iList are used.
%
%iList = {'IMG_0001.JPG' 'IMG_0002.JPG' 'IMG_0003.JPG' 'IMG_0004.JPG' 'IMG_0005.JPG' 'IMG_0006.JPG' 'IMG_0007.JPG' 'IMG_0008.JPG' 'IMG_0009.JPG' 'IMG_0010.JPG'};  % Specify image names. If iList is commented out then all images in folder inPath are used

% --------- end of user configurable parameters ------

% Search for jpeg and png files
if (~exist('iList') == 1)
    %fileData = rdir([inPath, '\*.*'], 'regexp(name, [''.png'' ''|'' ''.jpg''], ''ignorecase'')', true);
    %iList = {fileData.name};
    
    fileData = dir (inPath);
    flg = [fileData.isdir];
    iList = {fileData(not(flg)).name};
end

% Create PRNU name

switch pmode
    case 'averagePRNU'
        imageID = 'AVE';
        prnuName = {[cameraID '_' imageType '_' imageID]};
    case 'singlePRNU'
        nL = length(iList);
        for i=1:nL
            [b1 imageID b2] = fileparts(iList{i}); clear b1 b2;
            prnuName{i} = {[cameraID '_' imageType '_' imageID]};
        end
end

denoiseFolderRoot = [inPath '\denoise'];

addargs.filterName = denoiseFilter;
switch denoiseFilter
    
    case 'gaussian'    % Use Gaussian LPF 
    
        addargs.denoiseFolder = [denoiseFolderRoot '\gaussian'];     % Folder to save denoised images
        addargs.filterSize = [3 3];     % Filter size
        addargs.filterSD = 0.5;         % Gaussian filter standard deviation
        addargs.overwrite = overwrite;          % 0: use denoised file if exists, 1: overwrite file
        
        switch pmode
            case 'averagePRNU'
                prnuName{1} = [prnuName{1} '_G'];
            case 'singlePRNU'
                for i=1:nL
                    prnuName{i} = {[prnuName{i}{:} '_G']};
                end
        end
        
    case 'mihcak'       % Use Mihcak wavelet based denoising algorithm
        
        if ~wavInst
            fprintf('Warning: Mihcak filter selected but Wavelet Toolbox is not available.\n');
        end
        addargs.denoiseFolder = [denoiseFolderRoot '\mihcak'];     % Folder to save denoised images
        addargs.sigma0 = 5;     % value of sigma0
        addargs.overwrite = overwrite;
        addargs.saveWaveletCoeffs = saveWaveletCoeffs;      % Save wavelet coefficients as a .mat file
        
        switch pmode
            case 'averagePRNU'
                prnuName{1} = [prnuName{1} '_M'];
            case 'singlePRNU'
                for i=1:nL
                    prnuName{i} = {[prnuName{i}{:} '_M']};
                end
        end

    case 'sigma'        % Use Sigma filter
        
        addargs.denoiseFolder = [denoiseFolderRoot '\sigma'];     % Folder to save denoised images
        addargs.windowSize = 7;         % window size to use in sigma filter
        %addargs.stdval = [2.35 1.6 2.0];           % value for standard deviation of flat field images for Kodak V550
        addargs.stdval = [2.0 2.0 2.0];           % value for standard deviation of flat field images for Kodak V550
        addargs.overwrite = overwrite;          % 0: use denoised file if exists, 1: overwrite file
        
        switch pmode
            case 'averagePRNU'
                prnuName{1} = [prnuName{1} '_S'];
            case 'singlePRNU'
                for i=1:nL
                    prnuName{i} = {[prnuName{i}{:} '_S']};
                end
        end
        
    case 'bm3d'     % Use the method of denoising by sparse 3D transform-domain collaborative filtering
        
        addargs.denoiseFolder = [denoiseFolderRoot '\bm3dFolder'];
        addargs.sigma = 0.4;        % Use a standard deviation of 1.0 (for intensities in range [0 255]
        addargs.overwrite = overwrite;      % use denoised file if it exists
        
        switch pmode
            case 'averagePRNU'
                prnuName{1} = [prnuName{1} '_B'];
            case 'singlePRNU'
                for i=1:nL
                    prnuName{i} = {[prnuName{i}{:} '_B']};
                end
        end
        
    otherwise
        
        fprintf('Unknown filter\n');

end

switch pmode
    case 'averagePRNU'
        prnuNameID = {[prnuName{1} '_L']};
        status = prnu_extraction(opMode, inPath, outPath, prnuNameID, iList, addargs);
    case 'singlePRNU'
        for i=1:nL
            prnuNameID = {[prnuName{i}{:} '_L']};
            iListID = iList(i);
            status = prnu_extraction(opMode, inPath, outPath, prnuNameID, iListID, addargs);
        end
end


% ++++++++++++++++++++++++++++++ Cross Correlation +++++++++++++++++++++++++++++

case 'correlate'

% ------- User configurable parameters for 'correlate' mode ----------
    
inPath = 'c:\images\prnu\reference';

outPath = 'c:\images\prnu\results';      % This is where the correlation results file is written

% If correlating a reference frame with a SINGLE noise residual:
%    1) Comment out rListFile
%    2) Set refFileNames to {<name of reference PRNU file>, <name of noise residual>}
% After running the script, the correlation results for the red, green and
% blue channels are written to Matlab Command Window. No files are saved.
%
% If correlating a reference PRNU frame with MULTIPLE noise residuals:
%   1) Create a text file containing two columns. The first column should contain the
%      name of the reference file, and the second the name of the noise residual
%      file. (See accompanying readme.txt).
%   2) Set rListFile to the path and name of the correlation list text file
%      created in step 1. The variable refFilenames is ignored, if it is
%      defined.
% After running the script, the cross-correlation results are written to a
% file in the folder specified by outPath.
% 

rListFile = 'c:\images\prnu\correlation-list.txt';   % Comment out if individual prnu reference filenames in refFileNames, below, are to be used instead.

refFileNames = {'Kodak-V550-S_FF_AVE_S_L.mat', 'Kodak-V550-B_NI_100_1697_S_L.mat'}; % This is ignored if rListFile above is defined

% --------- end of user configurable parameters ------

if exist('rListFile') == 1
    
    fid = fopen(rListFile, 'r');

    % Read pairs of prnu reference filenames from list.
    rList = textscan(fid, '%s %s');

    rList1 = rList{1};
    rList2 = rList{2};

    nList = size(rList1, 1);    % Number of entries

    % Create a file for wrinting correlation results
    dateNow = datestr(now, 30);
    
    cResultsFileName = fullfile(outPath, ['corlnResults-' datestr(now, 30) '.txt']);
    fidw = fopen(cResultsFileName, 'wt');
    
    fprintf(fidw, 'Correlation results on: %s\n', datestr(datenum(dateNow, 'yyyymmddTHHMMSS'), 0));
    fprintf(fidw, 'Reference folder: %s\n', inPath);
    
    for i = 1:nList
        prnuName = {rList1{i}, rList2{i}};
        [status, ref1, ref2, colrNames, corr] = prnu_extraction(opMode, inPath, outPath, prnuName);
        fprintf(fidw, '%s %s ', strtok(rList1{i}, '.'), strtok(rList2{i}, '.'));
        
        for j=1:length(colrNames)
            fprintf(fidw, '%s: %6.4f  ', colrNames{j}, corr(j));
        end
            
        fprintf(fidw, ' Sum: %6.4f\n', sum(corr));
    end

    fclose(fidw);
    
else
    prnuName = refFileNames;
    status = prnu_extraction(opMode, inPath, outPath, prnuName);

end

end