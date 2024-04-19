function ReconData = read_kspace(sourceDir, rawfileName, fast_read_options)
%% Reads in data from .dat
% ReconData = read_kspace(sourceDir, saveDir, rawfileName, fast_read_options)
% 
% Reads in raw kspace data from .dat file using 'fast_read_ve11.m according 
% to options provided by fast_read_options. Reads in metadata from .dat 
% file using multiple subfunctions (jGetFileInfo.m, fast_read_ve11.m, and 
% mapVBVD_uk.m). Outputs all data and metadata to struct (ReconData).
%
% Input
% -----
% sourceDir : string
%   Pathname of raw data file to read in
% rawfileName : string
%   Name of .dat file to read in
% fast_read_options : cell
%   Cell array of strings of size {1,N}. List of input options (varargin)
%   for 'fast_read_ve11.m'
%
% Output
% ------
% ReconData : struct 
%   Struct with fields containing the kspace data and sequence metadata
%   ReconData.ksdata :  complex-double array (k-space data)
%   ReconData.header :  struct  (metadata from 'jGetFileInfo.m',
%                                'fast_read_ve11.')
%   ReconData.noiseData : struct  (pre-scan NoiseAdjust data)
%   ReconData.dims : cell array  (dimension order of ksdata)
%   ReconData.fastreadOptions : cell array  (fast_read_options input var) 
%   ReconData.shd : struct (metadata from mapVDVB_uk.m)
%
% Example
% -------
% sourceDir = '/System/Volumes/Data/v/raid10/users/sjohnson/Matlab...
%                Code/mri_recon_framework/mri-recon/test-data/seg-epi/multichannel-phantom/';
% rawfileName = 'meas_MID00262_FID27298_BRIFU_MRTI_XA50_1p2x0p6x1p8_ETL5_paddleCoil.dat';
% fast_read_options = {'NoFillToFullFourier','ShiftDCToMatrixCenter','NoGUI',...
%                     'LastMeasDatOnly','CollapseSeg','ReadNoiseAdj'};
% >>ReconData = read_kspace(sourceDir, rawfileName, fast_read_options);
%% Created 2024-04-02 Sara L. Johnson, Sam Adams-Tew


%%%~~~ Inspect variables, specify defaults ~~~~%%%

% Specify default fast_read_ve11 options
if ~exist('fast_read_options', 'var')  
    fast_read_options = {'NoFillToFullFourier',...
                        'ShiftDCToMatrixCenter',...
                        'NoGUI',...
                        'LastMeasDatOnly',...
                        'CollapseSeg',...
                        'ReadNoiseAdj'};
end 

%%%~~~ Begin reading data ~~~~%%%

fullfilename = fullfile(sourceDir, rawfileName); 

% Get k-space 'header' from .dat metadata (J. Mendez)
fileInfo=jGetFileInfo(fullfilename);
header = fileInfo.Protocol(end); 

% Check if MRI acquisition isMultislice
if (length(header.sSliceArray.asSlice)) > 1
    isMultislice = 1;
else 
    isMultislice = 0;
end

% FAST_READ_ve11: Get k-space data from .dat
if isMultislice

    % for Multislice, get additional info from .dat metadata
    fast_read_options = cat(2, fast_read_options,...
                               {'ReadRelSliceNumber',...
                                'ReadSVector',...
                                'ReadQuaternion'}); 

    % read in kspace and Multislice metadata
    [measData, relSliceNo] = fast_read_ve11(fullfilename,...
                                            fast_read_options{:}); 

    % append Multislice metadata to 'header'
    header.multiSlice.relSliceNumber = relSliceNo; 
    header.multiSlice.sliceVectors = measData{end}.SVector; 
    header.multiSlice.slabQuaternion = measData{end}.Quaternion; 

else

    % read in kspace 
    measData = fast_read_ve11(fullfilename,...
                              fast_read_options{:}); 
end 

% mapVBVD_UK: Get additional metadata from .dat
twix_obj = mapVBVD_uk(fullfilename);

if size(twix_obj, 2) > 1
    shd = getshd(twix_obj{1}, twix_obj{end});
else
    shd = getshd(twix_obj, twix_obj);
end


%%%~~~ Store data output ~~~~%%%

% Store kspace data and metadata 
ReconData.ksdata = measData{end}.data;
ReconData.header = header; 
ReconData.noiseData = measData{1}.noiseAdj;
ReconData.dims = measData{end}.dim;
ReconData.fastreadOptions = fast_read_options; 
ReconData.shd = shd;


end 