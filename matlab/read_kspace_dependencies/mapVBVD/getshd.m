%function [ RawData, shd ] = getshd( twix_obj_1,twix_obj_end,shd,par )
function [ shd ] = getshd( twix_obj_1,twix_obj_end,shd,par )
%UNTITLED Summary of this function goes here
%getshd: gets the short header variables from the open file
%   Detailed explanation goes here
    shd.foundnoise = 0;
    try
        shd.noise1 = double(twix_obj_1.noise{''});
        shd.dimvn1 = twix_obj_1.noise.sqzDims;
        shd.foundnoise = 1;
    catch
        disp('no noise scan in this data in twix_obj{1, 1}');
    end
    try
        shd.noise2 = double(twix_obj_end.noise{''});
        shd.dimvn2 = twix_obj_end.noise.sqzDims;
        shd.foundnoise = shd.foundnoise + 2;
    catch
        disp('no noise scan in this data in twix_obj_end');
    end
    %RawData = double(twix_obj_end.image{''});  %VB17. HOD Commented out
shd.dimv = lower(twix_obj_end.image.sqzDims);
shd.ncol = twix_obj_end.image.NCol;
shd.nch = twix_obj_end.image.NCha;
shd.nlin = twix_obj_end.image.NLin;
shd.npar = twix_obj_end.image.NPar;
shd.nsli = twix_obj_end.image.NSli;
shd.nav = twix_obj_end.image.NAve;
shd.nphs = twix_obj_end.image.NPhs;
shd.neco = twix_obj_end.image.NEco;
shd.nrep = twix_obj_end.image.NRep;
shd.nset = twix_obj_end.image.NSet;


          
shd.coils = cell(shd.nch);
shd.Txfrequency = twix_obj_end.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1, 1}.lFrequency;  
shd.Txfrequency1 = twix_obj_1.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1, 1}.lFrequency;  
shd.TxfreqDicom = twix_obj_end.hdr.Dicom.lFrequency;  
shd.TxfreqDicom1 = twix_obj_1.hdr.Dicom.lFrequency;  
shd.TxfreqPheonix = twix_obj_end.hdr.Phoenix.sTXSPEC.asNucleusInfo{1, 1}.lFrequency;  
shd.TxfreqPheonix1 = twix_obj_1.hdr.Phoenix.sTXSPEC.asNucleusInfo{1, 1}.lFrequency;  
shd.TR=twix_obj_end.hdr.Config.TR/1000;    %TR in ms
shd.contrasts = twix_obj_end.hdr.MeasYaps.lContrasts;  %should be number of echoes
shd.etl = twix_obj_end.hdr.MeasYaps.sFastImaging.lEPIFactor;           %lTurboFactor;       %Dicom.EchoTrainLength;
%shd.EchoSpacing    = twix_obj_end.hdr.Meas.lEchoSpacing;
shd.TE=twix_obj_end.hdr.MeasYaps.alTE{1,1}/1000; % TE in ms
vvv = cell2mat(twix_obj_end.hdr.MeasYaps.alTE);
shd.TEv= vvv(1,1:shd.contrasts)/1000; % TE in ms
shd.phaseFOV = twix_obj_end.hdr.Config.PhaseFoV;
shd.readoutFOV = twix_obj_end.hdr.Config.ReadFoV;
shd.sliceFOV = twix_obj_end.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness; %Attention currently it is assumed that all slices have same thickness
shd.BaseResolution = twix_obj_end.hdr.MeasYaps.sKSpace.lBaseResolution;
shd.PhaseRes = twix_obj_end.hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
shd.PhaseNlines= twix_obj_end.hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
shd.nslice = twix_obj_end.hdr.MeasYaps.sSliceArray.lSize;
          if(shd.npar > shd.nslice)         %probably should be whether shd.do3D is set or 'par' is set
              shd.SliceRes = shd.npar;
          else
              shd.SliceRes = shd.nslice;
          end
          shd.AccelFactPE = twix_obj_end.hdr.MeasYaps.sPat.lAccelFactPE;
          shd.dograppa = 0;
          if(shd.AccelFactPE > 1)
              shd.dograppa = 1;
              shd.noise2 = double(twix_obj_end.noise{''});
              shd.refscan = double(twix_obj_end.refscan{''});
              shd.dimvref = lower(twix_obj_end.refscan.sqzDims);
          end

cearray = repmat('aaa',[shd.nch,1]);
iadc = zeros(shd.nch,1);
try         %DLP 4/19/19 changed below to handle VB17 body coil
    aaa = char(twix_obj_end.hdr.MeasYaps.asCoilSelectMeas{1}.asList{jch}.sCoilElementID.tElement);
    for jch = 1:shd.nch
        aaa = char(twix_obj_end.hdr.MeasYaps.asCoilSelectMeas{1}.asList{jch}.sCoilElementID.tElement);
        [~,nchar] = size(aaa);
        nchce = min([nchar,3]);
        cearray(jch,1:nchce) = aaa(1,1:nchce);
    end
catch ME
    ME.message
    for jch = 1:shd.nch
        aaa = char(twix_obj_end.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1, 1}.asList{jch}.sCoilElementID.tElement);
        iadc(jch) = twix_obj_end.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1, 1}.asList{jch}.lADCChannelConnected;
        [~,nchar] = size(aaa);
        if length(aaa)==2 %% HOD 2020-20-30
            cearray(jch,1:2) = aaa(1,1:2);
        else
            cearray(jch,:) = aaa(1,1:3);
        end
        if(nchar>3) cearray(jch,:) = aaa(1,2:4);end
    end
end
[iadcs,jind] = sort(iadc);
shd.ce = char(cearray(jind,:));

shd.B0      = twix_obj_end.hdr.Dicom.flMagneticFieldStrength;
shd.FA      = twix_obj_end.hdr.MeasYaps.adFlipAngleDegree{1,1};
shd.txamp   = twix_obj_end.hdr.Dicom.flTransRefAmpl;
try shd.txampNL = twix_obj_end.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude; end
shd.TXSPEC_RFPULSE_all = twix_obj_end.hdr.MeasYaps.sTXSPEC.aRFPULSE;
shd.ntime   = shd.nrep * shd.nav;
shd.dim     = twix_obj_end.hdr.Dicom.tMRAcquisitionType; % HOD
try
    shd.nbaseline = twix_obj_end.hdr.MeasYaps.sWipMemBlock.alFree{1,3}  ;
    shd.wipmemblock = twix_obj_end.hdr.MeasYaps.sWipMemBlock.alFree{1,:}  ;
catch
    shd.nbaseline = 1;
    shd.wipmemblock = [0,0];
end
%shd.wipmemblock = twix_obj_end.hdr.MeasYaps.sWipMemBlock.alFree{1,3}  ;
%get the array variables for time stamps:
shd.dwellTime = twix_obj_end.hdr.MeasYaps.sRXSPEC.alDwellTime{1} ;
%shd.rbw = 1e9 ./ (size(RawData,1)*twix_obj_end.hdr.MeasYaps.sRXSPEC.alDwellTime{1}) ; % HOD Commented out
shd.rbw     = 1e9 ./ ((shd.ncol)*twix_obj_end.hdr.MeasYaps.sRXSPEC.alDwellTime{1}) ; % Use nCol instead of size(RawData,1)
shd.EchoSpacing = 1/shd.rbw; % Correct?!
shd.ulTimeStamp = twix_obj_end.hdr.ulTimeStamp; % DEP 01/08/19
shd.ulAcquisitionTimeStamp = twix_obj_end.hdr.ulAcquisitionTimeStamp; % DEP 01/08/19
shd.ulAcquisitionStrTimeStamp = twix_obj_end.hdr.ulAcquisitionStrTimeStamp; % DEP 01/08/19
shd.lTotalScanTimeSec = twix_obj_end.hdr.MeasYaps.lTotalScanTimeSec; % DEP 01/08/19
try
    %shd.Tnav1 = twix_obj_end.hdr.MeasYaps.sWipMemBlock.alFree{par.jindTEnav1};       %navigator TE1 in us   
    %shd.Tnav2 = twix_obj_end.hdr.MeasYaps.sWipMemBlock.alFree{par.jindTEnav2};       %navigator TE2 in us
catch ME
    ME.message
end

end

