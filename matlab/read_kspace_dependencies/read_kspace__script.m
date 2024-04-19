addpath(genpath('/System/Volumes/Data/v/raid10/users/sjohnson/Matlab Code/mri_recon_framework/mri-recon/matlab/read_kspace_dependencies/'))

sourceDir = '/System/Volumes/Data/v/raid10/users/sjohnson/Matlab Code/mri_recon_framework/mri-recon/test-data/seg-epi/multichannel-phantom/';
saveDir = '/System/Volumes/Data/v/raid10/users/sjohnson/Matlab Code/mri_recon_framework/matlab-sara-WIP/';

rawfileName = 'meas_MID00262_FID27298_BRIFU_MRTI_XA50_1p2x0p6x1p8_ETL5_paddleCoil.dat';

fast_read_options = {'NoFillToFullFourier',...
                    'ShiftDCToMatrixCenter',...
                    'NoGUI',...
                    'LastMeasDatOnly',...
                    'CollapseSeg',...
                    'ReadNoiseAdj'};


ReconData = read_kspace(sourceDir, rawfileName, fast_read_options); 

save([saveDir 'ks_' rawfileName(1:end-4)], 'ReconData'); 
