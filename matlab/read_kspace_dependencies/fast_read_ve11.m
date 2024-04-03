% =============================================================================
function [KSpace, varargout] = fast_read_ve11(filenameIn,varargin)
% -----------------------------------------------------------------------------
% fast_read_ve11  Read a VE11 measurement data file with the property that all 
%                  lines for a given measurement have the same number of 
%                  samples.  This differs slightly from the VB reader 
%                  where only one measurement was stored per file.  VD can
%                  store multiple measurementhes.
%
% KSpace = fast_read_ve11() 
%    prompts user with gui and returns kSpace cell array, each cell storing
%               the kSpace from a particular measurement.
%
% KSpace = fast_read_ve11(filename, 'option1', 'option2', 'option3', ....)
%    reads filename and returns kspace with options applied where
%    options are strings
%
% fast_read_ve11 with no arguments or return values will list options
%
% [KSpace Other] = fast_read_ve11(filename, 'ReadOther') 
%    reads filename and returns kSpace as well as Other if it exists where
% -----------------------------------------------------------------------------

   global MDH_SCANSIZE MDH_CHANSIZE MAX_BYTES_PER_FREAD
   global index_keys
   global mdhColumns mdhBitFlag

   %% Nested functions (share parent's scope)
   function set_outputs_to_zero
      KSpace{1}.data = [];
      for iArgLocal=1:nargout-1
         varargout{iArgLocal} = [];
      end
      if ~(Control.User.noGui || Control.User.silent)
         if ~isempty(Control.GUI.hWaitBar)
            close(Control.GUI.hWaitBar);
         end
      end
   end

   if (nargin == 0) && (nargout == 0)
      UserInputDefault = define_user_input_default_struct_array;
      UserInputDefault = derive_user_input_default_fieldname(UserInputDefault);
      print_available_options(UserInputDefault);
      return
   end

   % Initializing
   KSpace{1}.data = [];
   for iArg=1:nargout-1
      varargout{iArg} = [];
   end
   nOptionsIn = nargin-1;
   if nOptionsIn > 0
      optionsIn = varargin;
   else
      optionsIn = [];
   end
   set_global_mdh_parameters(filenameIn);

   Control = get_controls_from_user_input(filenameIn, optionsIn, nargin, nargout);
   if isempty(Control)
      set_outputs_to_zero;
      return
   end
   if Control.User.readFileInfo
      if length(which('jGetFileInfo')) < 1
         disp('jGetFileInfo not in path');
         set_outputs_to_zero;
         return
      end
   end

   if ~(Control.User.noGui || Control.User.silent)
      Control.GUI.hWaitBar = waitbar(0, 'Extracting TextHeader', 'Name', 'Progress', 'Resize', 'on');
      set(findall(Control.GUI.hWaitBar), 'Units', 'Normalized');
      set(Control.GUI.hWaitBar, 'Units', 'Pixels', 'Position', [100 100 400 150]);
      movegui(Control.GUI.hWaitBar, 'center');
   end

   % Read noise if indicated
   if Control.User.readNoiseAdj
      mdhBitFlagPrev = mdhBitFlag;
      KNoise = fast_read_ve11(filenameIn, 'KeepBitMask', '2000000', 'LogicKeepOnly', 'NoFillToFullFourier', 'Silent');
      mdhBitFlag= mdhBitFlagPrev;
   end

   fid = fopen(Control.File.In.name,'r','ieee-le');

   % Confirm that we have a ve11 file
   [isVE11, MrParcRaidFileHeader, MrParcRaidFileCell] = is_ve11_file(Control, fid);
   if ~isVE11
      set_outputs_to_zero;
      return
   end

   % Swap SEG and ECO indices
   % This was introduced to handle HIFU meas dat where they used PHASCOR flags
   %    to store special navigators.  Unfortunately, they were stored in a way
   %    contrary to how PHASCOR are stored in general.  To read them, turn on
   %    the KeepBitFlag for phasecor and add 'SwapSegAndEco' to the call parameters.
   if Control.User.swapSegAndEco
      temp = mdhColumns.scan.indexECO;
      mdhColumns.scan.indexECO = mdhColumns.scan.indexSEG;
      mdhColumns.scan.indexSEG = temp;
      clear temp
   end

   % Loop through measurements
   RelativeSliceNumber = {};
   WiPMemBlock = {};
   SliceArray = {};
   indexOutTemp = uint64([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
   for iMeas = 1:MrParcRaidFileHeader.nMeas

      if Control.User.readNoiseAdj 
         KSpace{iMeas}.noiseAdj = KNoise{iMeas};
      end

      nNormalThisMeas = 0;
      if Control.User.lastMeasDatOnly && (iMeas < MrParcRaidFileHeader.nMeas)
         KSpace{iMeas}.data = [];
         continue
      end

      % Parse the text header for dimensions, scale factors, control, etc.
      [Control, MrParcRaidFileCell, Dim, CoilSelectMap, TextHeader, ...
       CoilSelectMeasMap, RelativeSliceNumberLocal, WiPMemBlockLocal, SliceArrayLocal, MeasYapsAsc] = ...
          parse_meas_text_header(Control,fid,MrParcRaidFileCell,iMeas);
      if isempty(Control)
         set_outputs_to_zero;
         return
      end
      RelativeSliceNumber{iMeas} = RelativeSliceNumberLocal;
      WiPMemBlock{iMeas} = WiPMemBlockLocal;
      SliceArray{iMeas} = SliceArrayLocal;

      % First pass to chunk this measurement into blocks with the same scan lengths
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(0.5, Control.GUI.hWaitBar, 'Measuring number of scans');
      end

      % With readout and PMU being different sizes, 
      %    we now have to chunk in blocks of same sized
      %    lines.  That requires scanning all line lengths
      %    up front.
      if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
          evalMaskArray = chunk_into_equal_line_length_blocks2(fid,MrParcRaidFileCell, iMeas, Dim, Control);
          display_eval_mask(evalMaskArray,iMeas);
          continue
      else
         [~, sizeAllScansInBytes, sizeChanInBytes, MdhBlock, ...
          nBlockSYNCDATA, ~, ~, ~, ...
          nScanInMeas, ~, nSamplesInScan] = chunk_into_equal_line_length_blocks2(fid,MrParcRaidFileCell, iMeas, Dim, Control);
         if isempty(sizeAllScansInBytes)
            KSpace{iMeas}.data = [];
            continue
         end
      end

      % Second pass, Part 1 to extract PMU
      if (nBlockSYNCDATA > 0) && MrParcRaidFileCell{iMeas}.Control.User.readPMU
         PMUOutCell{iMeas} = extract_pmu(fid,MdhBlock,Control);
         if isempty(PMUOutCell{iMeas})
            set_outputs_to_zero;
            return
         end
      end % if we're loading PMU

      % Second pass, Part 2 to assemble dimensions and read order for readouts
      [dimNames, nScan, n, kSpaceCentreLineNo, kSpaceCentrePartitionNo, ...
         order, dimOut, dimOutProd, channelIDUniq] = ...
         scan_mdh_for_dim_and_order(fid,Dim,MrParcRaidFileCell{iMeas}.Control,sizeChanInBytes,MdhBlock,nScanInMeas,nSamplesInScan);
      if isempty(dimNames)
         set_outputs_to_zero;
         return
      end

      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(1.0, Control.GUI.hWaitBar, 'Output dimensions computed.');
      end

      % Third pass to read the data
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(0.0, Control.GUI.hWaitBar, 'Allocating arrays');
      end

      % Allocate arrays
      %keyboard
      progressPercent = 0.0;
      mdhExtrasToRead = {};
      if nScan > 0 && ~Control.User.noKSpace
         Control.User.foundNormal = true;
         KSpace{iMeas}.data = complex(single(0) * zeros(n.COL, compute_product(dimOut), 'single'),single(0));
         if Control.User.readIceProgramParam
             mdhExtrasToRead = [mdhExtrasToRead {'iceProgramPara'}];
             KSpace{iMeas}.IceProgramPara = zeros(28, compute_product(dimOut),'uint16');
         end
         for readStrStub = {'SVector', 'Quaternion'}
            readStr = ['read' readStrStub{1}];
            if Control.User.(readStr)
               mdhExtrasToRead = [mdhExtrasToRead {[lower(readStrStub{1}(1)) readStrStub{1}(2:end)]}];
               KSpace{iMeas}.(readStrStub{1}).vec2index = {};
               KSpace{iMeas}.(readStrStub{1}).index2vec = {};
               if strcmp(readStrStub{1},'Quaternion')
                   KSpace{iMeas}.(readStrStub{1}).index2rot = {};
               end
               KSpace{iMeas}.(readStrStub{1}).index = uint8(0) * zeros(1, compute_product(dimOut),'uint8');
            end
         end
         for readStrStub = {'TimeStamp', 'PMUTimeStamp', 'TimeSinceLastRF', 'SequenceTime'}
            readStr = ['read' readStrStub{1}];
            if Control.User.(readStr)
               mdhExtrasToRead = [mdhExtrasToRead {[lower(readStrStub{1}(1)) readStrStub{1}(2:end)]}];
               KSpace{iMeas}.(readStrStub{1}).count = 0 * zeros(compute_product(dimOut),1);
               KSpace{iMeas}.(readStrStub{1}).value = uint32(0) * zeros(compute_product(dimOut),1,'uint32');
            end
         end
      else
          KSpace{iMeas}.data = [];
      end

      progressPercent = progressPercent + 0.25;
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(progressPercent, Control.GUI.hWaitBar, 'Allocating arrays');
      end
      
      iScanRead = uint64(0);
      timeRemain = uint64(10000000);
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         time0 = tic;
         time00 = tic;
      end

      for iBlock=1:length(MdhBlock)

         if Control.User.noKSpace
             continue
         end

         % Skip SYNCDATA
         if MdhBlock(iBlock).isSyncData || MdhBlock(iBlock).isAcqEnd || ~MdhBlock(iBlock).keepScan 
            continue
         end
         iReadout = uint64(1):uint64(MdhBlock(iBlock).nSamplesInScan);

         maxScanPerRead = max(idivide(MAX_BYTES_PER_FREAD,MdhBlock(iBlock).dmaLength,'fix'),uint64(1));
         nScanToRead = min(maxScanPerRead,MdhBlock(iBlock).nScan);

         fseek(fid,MdhBlock(iBlock).offsetInFile,'bof');
         nScanRemainInBlock = MdhBlock(iBlock).nScan;
      
         iScanReadBlock = uint64(0);
         while nScanRemainInBlock > 0

            if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
               dTimeElapsed = toc(time0);
               rate = dTimeElapsed / (double(iScanRead) + 1);
               if dTimeElapsed > 20
                  timeRemain = min(double(nScanInMeas - iScanRead) * rate,timeRemain);
               else
                  timeRemain = double(nScanInMeas - iScanRead) * rate;
               end
               %textMessage = sprintf('Pass 3 of 3\nReading scans\nNormal lines so far:%d\nTime remaining (s): %d', nNormalThisMeas, uint32(timeRemain));
               textMessage = sprintf('Pass 3 of 3\nReading scans\nTime remaining (s): %d', uint32(timeRemain));
               waitbar(double(iScanRead)/double(nScanInMeas), Control.GUI.hWaitBar, textMessage);
            end
            blockOfScans = fread(fid,[MdhBlock(iBlock).dmaLength nScanToRead],'uchar=>uchar');
            MdhParam = extract_mdh_parameters(blockOfScans, 'scan', [mdhExtrasToRead {'cutOffDataPre', 'cutOffDataPost'} index_keys]);
            if MrParcRaidFileCell{iMeas}.Control.User.collapseSeg
               MdhParam.indexSEG = MdhParam.indexSEG * uint16(0);
            end
            MdhBitFlags = extract_mdh_bitflags(blockOfScans, {'isAcqEnd', 'bKeepEIM', 'bDropEIM', 'isNoiseAdj', 'isPhaseCor', 'isRawDataCorrection','isReflect'});

            for indexStr = index_keys
               MdhParam.(indexStr{1}) = uint64(MdhParam.(indexStr{1}));
            end
            channelID = uint16(0) * zeros(MdhBlock(iBlock).nChannelUsed,nScanToRead,'uint16');
            for iCha = 1:MdhBlock(iBlock).nChannelUsed
               offset  = uint64(MDH_SCANSIZE) + uint64(iCha-1) * sizeChanInBytes;
               channelID(iCha,:) = reshape(extract_single_mdh_parameter(blockOfScans(offset+1:offset+sizeChanInBytes,:),'chan','channelID'),[1,nScanToRead]);
            end
            if Control.User.readSequenceTime
               iCha = 1;
               offset  = uint64(MDH_SCANSIZE) + uint64(iCha-1) * sizeChanInBytes;
               sequenceTime = bitshift(extract_single_mdh_parameter(blockOfScans(offset+1:offset+sizeChanInBytes,:),'chan','sequenceTime'),-9);
            end
            sizeSamplesInBytes = MdhBlock(iBlock).nSamplesInScan * uint64(8);
            for iCha = 1:MdhBlock(iBlock).nChannelUsed

               channelIDOut = squeeze(channelID(iCha,:));
               offset  = uint64(MDH_SCANSIZE) + uint64(iCha-1) * sizeChanInBytes + MDH_CHANSIZE + uint64(1);

               tempArr = blockOfScans(offset:offset+sizeSamplesInBytes-1,:);
               tempArr = reshape(typecast(tempArr(:),'single'), [2*MdhBlock(iBlock).nSamplesInScan nScanToRead]);
               tempArr = complex(tempArr(1:2:end,:),tempArr(2:2:end,:));

               % Handle all the reflections at once
               indexReflect = find(MdhBitFlags.isReflect);
               if ~isempty(indexReflect)
                  tempArr(:,indexReflect) = flipud(tempArr(:,indexReflect));
               end

               % Handle Pre-cutoff
               indexPre = find(MdhParam.cutOffDataPre > 0 & ~MdhBitFlags.isNoiseAdj);
               if ~isempty(indexPre)
                  cutOffDataPreUniq = unique(MdhParam.cutOffDataPre(indexPre));
                  for iUniq = 1:numel(cutOffDataPreUniq)
                     indexCut = find(MdhParam.cutOffDataPre(indexPre) == cutOffDataPreUniq(iUniq));
                     tempArr(1:cutOffDataPreUniq(iUniq),indexPre(indexCut)) = ...
                        complex(zeros(cutOffDataPreUniq(iUniq),numel(indexCut),'single'),single(0));
                  end % loop over uniq cuts
               end

               % Handle Post-cutoff
               indexPost = find(MdhParam.cutOffDataPost > 0 & ~MdhBitFlags.isNoiseAdj);
               if ~isempty(indexPost)
                  cutOffDataPostUniq = unique(MdhParam.cutOffDataPost(indexPost));
                  for iUniq = 1:numel(cutOffDataPostUniq)
                     indexCut = find(MdhParam.cutOffDataPost(indexPost) == cutOffDataPostUniq(iUniq));
                     tempArr(end-cutOffDataPostUniq(iUniq)+1:end,indexPost(indexCut)) = ...
                        complex(zeros(cutOffDataPostUniq(iUniq),numel(indexCut),'single'),single(0));
                  end % loop over uniq cuts
               end

               maskKeep = ~MdhBitFlags.isAcqEnd;
               index = find(maskKeep);
               if ~(numel(index) > 0)
                  continue;
               end
               % Pass through to sort optimally
               if Control.User.applySortAlgorithm
                  kScanTest = zeros(numel(index),1,'uint64');
                  for kScan = 1:numel(index)
                     jScan = index(kScan);
                     indexCHA = find(channelIDUniq == channelIDOut(jScan)) - 1;
                     indexOutTemp(1) = indexCHA;
                     indexOutTemp(2) = MdhParam.indexLIN(jScan);
                     indexOutTemp(3) = MdhParam.indexSLC(jScan);
                     indexOutTemp(4) = MdhParam.indexPAR(jScan);
                     indexOutTemp(5) = MdhParam.indexACQ(jScan);
                     indexOutTemp(6) = MdhParam.indexECO(jScan);
                     indexOutTemp(7) = MdhParam.indexPHS(jScan);
                     indexOutTemp(8) = MdhParam.indexREP(jScan);
                     indexOutTemp(9) = MdhParam.indexSET(jScan);
                     indexOutTemp(10) = MdhParam.indexSEG(jScan);
                     indexOutTemp(11) = MdhParam.indexIDA(jScan);
                     indexOutTemp(12) = MdhParam.indexIDB(jScan);
                     indexOutTemp(13) = MdhParam.indexIDC(jScan);
                     indexOutTemp(14) = MdhParam.indexIDD(jScan);
                     indexOutTemp(15) = MdhParam.indexIDE(jScan);
                     indexOut = mod(indexOutTemp(order),dimOut);
                     kScanTest(kScan,1) = sum(indexOut(:) .* dimOutProd(:)) + uint64(1);
                  end
                  [~,kScanSortIndex] = sort(kScanTest);
               else
                  kScanSortIndex = 1:numel(index);
               end
               for kScan = 1:numel(index)
                  if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
                     dTimeElapsed00 = toc(time00);
                     if dTimeElapsed00 > 5
                        dTimeElapsed = toc(time0);
                        rate = dTimeElapsed / (double(iScanRead+kScan) + 1);
                        if dTimeElapsed > 20
                           timeRemain = min(double(nScanInMeas - iScanRead - kScan) * rate,timeRemain);
                        else
                           timeRemain = double(nScanInMeas - iScanRead - kScan) * rate;
                        end
                        time00 = tic;
                        %textMessage = sprintf('Pass 3 of 3\nReading scans\nNormal lines so far:%d\nTime remaining (s): %d', nNormalThisMeas, uint32(timeRemain));
                        textMessage = sprintf('Pass 3 of 3\nReading scans\nTime remaining (s): %d', uint32(timeRemain));
                        waitbar(double(iScanRead+kScan)/double(nScanInMeas), Control.GUI.hWaitBar, textMessage);
                     end
                  end
                  jScan = index(kScanSortIndex(kScan));
                  keyCHACoilSelectMeas = channelIDOut(jScan) + 1;
                  keyCHACoilSelect = CoilSelectMeasMap{keyCHACoilSelectMeas}.tElement;
                  indexCHA = find(channelIDUniq == channelIDOut(jScan)) - 1;
                  indexOutTemp(1) = indexCHA;
                  indexOutTemp(2) = MdhParam.indexLIN(jScan);
                  indexOutTemp(3) = MdhParam.indexSLC(jScan);
                  indexOutTemp(4) = MdhParam.indexPAR(jScan);
                  indexOutTemp(5) = MdhParam.indexACQ(jScan);
                  indexOutTemp(6) = MdhParam.indexECO(jScan);
                  indexOutTemp(7) = MdhParam.indexPHS(jScan);
                  indexOutTemp(8) = MdhParam.indexREP(jScan);
                  indexOutTemp(9) = MdhParam.indexSET(jScan);
                  indexOutTemp(10) = MdhParam.indexSEG(jScan);
                  indexOutTemp(11) = MdhParam.indexIDA(jScan);
                  indexOutTemp(12) = MdhParam.indexIDB(jScan);
                  indexOutTemp(13) = MdhParam.indexIDC(jScan);
                  indexOutTemp(14) = MdhParam.indexIDD(jScan);
                  indexOutTemp(15) = MdhParam.indexIDE(jScan);
                  temp = tempArr(:,jScan);
                  indexOut = mod(indexOutTemp(order),dimOut);
                  iScan = sum(indexOut(:) .* dimOutProd(:)) + uint64(1);
                  if MdhBitFlags.isRawDataCorrection(jScan) && Control.User.applyRawDataCorrection
                     temp = CoilSelectMap.(keyCHACoilSelect).rawDataCorrectionFactor*temp;
                  end
                  if Control.User.useMexInsertion
                     copy_column_a2b_single_complex(CoilSelectMeasMap{keyCHACoilSelectMeas}.flFFTCorrectionFactor*temp, ...
                                                    KSpace{iMeas}.data,double(iReadout(1)),double(iScan));
                  else
                     temp = CoilSelectMeasMap{keyCHACoilSelectMeas}.flFFTCorrectionFactor*temp;
                     KSpace{iMeas}.data(iReadout,iScan) = temp;
                  end
                  nNormalThisMeas = nNormalThisMeas + length(indexOut(:));
                  if Control.User.readIceProgramParam
                      KSpace{iMeas}.IceProgramPara(:,iScan) = MdhParam.iceProgramPara(:,jScan);
                  end
                  for readStrStub = {'SVector', 'Quaternion'}
                     readStr = ['read' readStrStub{1}];
                     mdhStr = [lower(readStrStub{1}(1)) readStrStub{1}(2:end)];
                     if Control.User.(readStr)
                        if strcmp(readStrStub{1}, 'SVector')
                            vectorStr = num2str(MdhParam.(mdhStr)((jScan-1)*3+1:jScan*3)','%f_%f_%f');
                        else
                            vectorStr = num2str(MdhParam.(mdhStr)((jScan-1)*4+1:jScan*4)','%f_%f_%f_%f');
                        end
                        vectorStr = strrep(vectorStr, '-', 'm');
                        vectorStr = strrep(vectorStr, '.', 'd');
                        vectorStr = ['v' vectorStr];
                        if length(KSpace{iMeas}.(readStrStub{1}).vec2index) > 0 && isfield(KSpace{iMeas}.(readStrStub{1}).vec2index, vectorStr)
                            iData = KSpace{iMeas}.(readStrStub{1}).vec2index.(vectorStr);
                        else
                            if length(KSpace{iMeas}.(readStrStub{1}).vec2index) > 0
                                iData = length(fieldnames(KSpace{iMeas}.(readStrStub{1}).vec2index)) + 1; 
                            else
                                iData = 1; 
                            end
                            KSpace{iMeas}.(readStrStub{1}).vec2index.(vectorStr) = iData;
                            if strcmp(readStrStub{1}, 'SVector')
                                KSpace{iMeas}.(readStrStub{1}).index2vec{iData} = MdhParam.(mdhStr)((jScan-1)*3+1:jScan*3)';
                            else
                                KSpace{iMeas}.(readStrStub{1}).index2vec{iData} = MdhParam.(mdhStr)((jScan-1)*4+1:jScan*4)';
                                [q_status, rot_matrix] = quaternion_to_rotation_matrix(MdhParam.(mdhStr)((jScan-1)*4+1:jScan*4)');
                                KSpace{iMeas}.(readStrStub{1}).index2rot{iData} = rot_matrix;
                            end
                        end
                        KSpace{iMeas}.(readStrStub{1}).index(iScan) = iData;
                     end
                  end
                  for readStrStub = {'TimeStamp', 'PMUTimeStamp', 'TimeSinceLastRF', 'SequenceTime'}
                     readStr = ['read' readStrStub{1}];
                     mdhStr = [lower(readStrStub{1}(1)) readStrStub{1}(2:end)];
                     if Control.User.(readStr)
                         KSpace{iMeas}.(readStrStub{1}).count(iScan) = KSpace{iMeas}.(readStrStub{1}).count(iScan) + 1;
                         nValues = KSpace{iMeas}.(readStrStub{1}).count(iScan);
                         if strcmp(readStrStub{1}, 'SequenceTime')
                            KSpace{iMeas}.(readStrStub{1}).value(iScan,nValues) = sequenceTime(iScan);
                         else
                            KSpace{iMeas}.(readStrStub{1}).value(iScan,nValues) = MdhParam.(mdhStr)(jScan);
                         end
                     end
                  end
               end % loop over kScan

            end % endfor over channels

            iScanRead = iScanRead + nScanToRead;
            iScanReadBlock = iScanReadBlock + nScanToRead;
            nScanRemainInBlock = MdhBlock(iBlock).nScan - iScanReadBlock;
            nScanToRead = min(nScanToRead,nScanRemainInBlock);

         end % while over sub-blocks of lines
      end % for over MdhBlocks
      clear blockOfScans;
      clear index;
      clear tempArr;
      clear channelID channelIDOut MdhParam MdhBitFlags;
      clear iScanRead iScanReadBlock nScanToRead;
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(0.5, Control.GUI.hWaitBar, 'All lines read from file.  Reshaping.');
      end

      % Reshape the results
      indexDimNonZero = find(dimOut > 1);
      if isempty(indexDimNonZero)
         indexDimNonZeroTemp = find(dimOut(1:end) > 0); 
         indexDimNonZero = zeros(1,1);
         indexDimNonZero(1) = indexDimNonZeroTemp(1);
      end
      if ~isempty(KSpace{iMeas}.data)
         KSpace{iMeas}.data = reshape(KSpace{iMeas}.data, [n.COL dimOut(indexDimNonZero)]);
         if Control.User.readIceProgramParam
            KSpace{iMeas}.IceProgramPara = reshape(KSpace{iMeas}.IceProgramPara, [28 dimOut(indexDimNonZero)]);
         end
         for readStrStub = {'SVector', 'Quaternion'}
            readStr = ['read' readStrStub{1}];
            if Control.User.(readStr)
                KSpace{iMeas}.(readStrStub{1}).index = reshape(KSpace{iMeas}.(readStrStub{1}).index, dimOut(indexDimNonZero));
            end
         end
         for readStrStub = {'TimeStamp', 'PMUTimeStamp', 'TimeSinceLastRF', 'SequenceTime'}
            readStr = ['read' readStrStub{1}];
            if Control.User.(readStr)
               dSize = max(KSpace{iMeas}.(readStrStub{1}).count(:));
               if n.CHA > 1
                  if length(indexDimNonZero(2:end)) > 1 
                     newDim = dimOut(indexDimNonZero(2:end));
                  else
                     newDim = [1 dimOut(indexDimNonZero(2:end))];
                  end
                  KSpace{iMeas}.(readStrStub{1}).count = reshape(KSpace{iMeas}.(readStrStub{1}).count(1:uint64(n.CHA):end), newDim);
                  KSpace{iMeas}.(readStrStub{1}).value = reshape(KSpace{iMeas}.(readStrStub{1}).value(1:uint64(n.CHA):end,:), [newDim dSize]);
               else
                  if length(indexDimNonZero) > 1 
                     newDim = dimOut(indexDimNonZero);
                  else
                     newDim = [1 dimOut(indexDimNonZero)];
                  end
                  KSpace{iMeas}.(readStrStub{1}).count = reshape(KSpace{iMeas}.(readStrStub{1}).count, newDim);
                  KSpace{iMeas}.(readStrStub{1}).value = reshape(KSpace{iMeas}.(readStrStub{1}).value, [newDim dSize]);
               end
            end
         end
      end

      if MrParcRaidFileCell{iMeas}.Control.User.shiftDCToMatrixCenter == true
         if nScan > 0 && ~Control.User.noKSpace
            allShift = dimOut * 0;
            colShift = idivide(n.COL - nSamplesInScan, uint64(2), 'fix');
            iLin = find(order == 2);
            if dimOut(iLin(1)) > 1
               linShift = uint64(dimOut(iLin(1))/2) - uint64(kSpaceCentreLineNo);
               allShift(iLin(1)) = linShift;
            end
            iPar = find(order == 4);
            if dimOut(iPar(1)) > 1
               parShift = uint64(dimOut(iPar(1))/2) - uint64(kSpaceCentrePartitionNo);
               allShift(iPar(1)) = parShift;
            end
            allShift = [colShift allShift(indexDimNonZero)];
            KSpace{iMeas}.data = circshift(KSpace{iMeas}.data,allShift);
         end
      end

      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(1.0, Control.GUI.hWaitBar, 'KSpace complete.');
      end

      tMsg = ['No scans read this measUID: ' MrParcRaidFileCell{iMeas}.Control.File.measUIDString];
      if nScan > 0
         tMsg = sprintf('KSpace dimensions:\n   COL=%d', n.COL);
         dimNamesOut = dimNames(order);
         for iDim=1:numel(indexDimNonZero)
            tMsg = strcat(tMsg, sprintf('\n   %s=%d', ...
                                        dimNamesOut{indexDimNonZero(iDim)}, ...
                                        dimOut(indexDimNonZero(iDim))));
         end
         KSpace{iMeas}.dim = ['COL', dimNamesOut(indexDimNonZero)];
      end

      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui, ...
                         'WAIT', tMsg,'Output KSpace Dimensions');
      end

      KSpace{iMeas}.TextHeader = TextHeader;
      KSpace{iMeas}.Protocol.Dim = Dim;
      KSpace{iMeas}.Protocol.Misc.kSpaceCentreLineNo = kSpaceCentreLineNo;
      KSpace{iMeas}.Protocol.Misc.kSpaceCentrePartitionNo = kSpaceCentrePartitionNo;
      KSpace{iMeas}.Protocol.CoilSelectMap = CoilSelectMap;
      KSpace{iMeas}.Protocol.CoilSelectMeasMap = CoilSelectMeasMap;
      KSpace{iMeas}.Protocol.MeasYapsAscConv = MeasYapsAsc;
      KSpace{iMeas}.Protocol.ChannelID = channelIDUniq;

      clear Dim CoilSelectMap CoilSelectMeasMap n MeasYapsAsc;

   end % endfor over all measurements in file
   if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
      close(Control.GUI.hWaitBar);
   end
   fclose(fid);

   if Control.User.readPMU
      varargout{Control.User.iArgOutPMU} = PMUOutCell;
   end
   if and(Control.User.foundSliceArray,Control.User.readSliceArray)
      varargout{Control.User.iArgOutSliceArray} = SliceArray;
   end
   if and(Control.User.foundWiPMemBlock,Control.User.readWiPMemBlock)
      varargout{Control.User.iArgOutWiPMemBlock} = WiPMemBlock;
   end
   if and(Control.User.foundRelSliceNumber,Control.User.readRelSliceNumber)
      varargout{Control.User.iArgOutRelSliceNumber} = RelativeSliceNumber;
   end
   if Control.User.readFileInfo
      varargout{Control.User.iArgOutFileInfo} = jGetFileInfo(filenameIn);
   end   

end

% =============================================================================
function UserInputDefault = define_user_input_default_struct_array
% -----------------------------------------------------------------------------
% Return struct array of user defaults:
%    UserInputDefault(i).varArgIn - string to expect from varargin
%    UserInputDefault(i).fieldName - defines structure field name to associate 
%                                with input.  If empty string, then same as 
%                                varArgIn with first letter lower case.
%    UserInputDefault(i).defaultValue - default value
%    UserInputDefault(i).type - currently: 'logical','hex2dec'
%    UserInputDefault(i).isArgOut - boolean
%    UserInputDefault(i).description
% -----------------------------------------------------------------------------

   userInputFieldNames = ...
         {{'varArgIn'},      {'fieldName'}, {'defaultValue'},    {'type'}, {'isArgOut'}, {'description'} };
   userInputValues = ...
      { ...
         {{'ApplyRawDataCorrection'}, {''},          {false}, {'logical'},      {false}, {'Apply raw data correction to (usually central TSE) lines based on mdh flag and gain parameters in the header.'}}, ...
         {{'ApplySortAlgorithm'},     {''},          {false}, {'logical'},      {false}, {'In some cases this algorithm might speed things up'}}, ...
         {{'ScanEvalInfoMaskOnly'},   {''},          {false}, {'logical'},      {false}, {'Scans the meas.dat file to record the set of evalinfomask bitmasks encountered.  Useful when later constructing keep/drop masks'}}, ...
         {{'DropBitMask'},            {''},     {'8000002604406'}, {'hex2dec'},      {false}, {'hex string bit mask for first evalInfoMask controlling what lines are dropped.'}}, ...
         {{'DropMatchAllBits'},       {''},          {false}, {'logical'},      {false}, {'MDH bits must match DropBitMask exactly to trigger true.  Otherwise any matched bit will trigger true.'}}, ...
         {{'KeepBitMask'},            {''},     {'fd9fbbf9'}, {'hex2dec'},      {false}, {'hex string bit mask for first evalInfoMask controlling what lines are kept.'}}, ...
         {{'KeepMatchAllBits'},       {''},          {false}, {'logical'},      {false}, {'MDH bits must match KeepBitMask exactly to trigger true.  Otherwise any matched bit will trigger true.'}}, ...
         {{'LogicKeepOnly'},          {''},          {false}, {'logical'},      {false}, {'Standard logic is (keep and not drop).  However, we can ignore the drop mask.  Lines are kept when keep condition is true.'}}, ...
         {{'LogicDropOnly'},          {''},          {false}, {'logical'},      {false}, {'Standard logic is (keep and not drop).  However, we can ignore the keep mask.  Lines are kept when drop condition is false.'}}, ...
         {{'LastMeasDatOnly'},        {''},          {false}, {'logical'},      {false}, {'Read only the last meas.dat block in the multi-measure meas.dat file.'}}, ...
         {{'NoFillToFullFourier'},    {''},          {false}, {'logical'},      {false}, {'Do NOT Fill k-space matrix out to complete FFT size in case of partial-Fourier acquisitions.'}}, ...
         {{'NoGui'},                  {''},          {false}, {'logical'},      {false}, {'Do not display GUI elements.'}}, ...
         {{'NoKSpace'},               {''},          {false}, {'logical'},      {false}, {'Do not bother reading k-space itself.  Ex: only headers.'}}, ...
         {{'ReadNoiseAdj'},           {''},          {false}, {'logical'},      {false}, {'Attempt to read noise adjust scans (bitmask=0x02000000).  Data will be attached as a .noise field to each kspace output.'}}, ...
         {{'ReadPMU'},                {''},          {false}, {'logical'},      {true},  {'Reads in the embedded PMU times and waveforms from the meas.dat header.  Usually requires gated sequences.'}}, ...
         {{'ReadRelSliceNumber'},     {''},          {false}, {'logical'},      {true},  {'Read relative slice number from header.  Necessary to sort 2D acquisitions.'}}, ...
         {{'ReadSequenceTime'},       {''},          {false}, {'logical'},      {false}, {'Sequence time (10us resolution) from the channel header. Only the first channel will be read.'}}, ...
         {{'ReadTimeStamp'},          {''},          {false}, {'logical'},      {false}, {'Time stamps (2.5ms) since midnight corresponding to KSpace lines.'}}, ...
         {{'ReadPMUTimeStamp'},       {''},          {false}, {'logical'},      {false}, {'2.5ms since last trigger (one-to-one with lines)'}}, ...
         {{'ReadTimeSinceLastRF'},    {''},          {false}, {'logical'},      {false}, {'Time stamps (2.5ms) since last RF (one-to-one with lines).'}}, ...
         {{'ReadIceProgramParam'},    {''},          {false}, {'logical'},      {false}, {'Read the free ice program para associated with each readout MDH (24 uint16 per readout + 4 so called reserved uint16).'}}, ...
         {{'ReadSVector'},            {''},          {false}, {'logical'},      {false}, {'Read the slice vector from each scan line'}}, ...
         {{'ReadQuaternion'},         {''},          {false}, {'logical'},      {false}, {'Read the quaternion from each scan line'}}, ...
         {{'ReadWiPMemBlock'},        {''},          {false}, {'logical'},      {true},  {'Extract WiPMemBlock from header.  User must provide return variable in var arg out.'}}, ...
         {{'ReadSliceArray'},         {''},          {false}, {'logical'},      {true},  {'Extract SliceArray from ASCCONV portion of header.'}}, ...
         {{'ReadFileInfo'},           {''},          {false}, {'logical'},      {true},  {'If you have jGetFileInfo in your path, it will be called to read and parse portions of the text header.'}}, ...
         {{'UseMexInsertion'},        {''},          {false}, {'logical'},      {false}, {'Use copy_column_a2b_single_complex mex program to speed up insertion of lines into kspace memory.  Requires having compiled copy_column_a2b_single_complex.'}}, ...
         {{'ResetFFTScale'},          {''},          {false}, {'logical'},      {false}, {'Reset scale for each coil (read from FFT Correction Factors) to one.'}}, ...
         {{'ShiftDCToMatrixCenter'},  {''},          {false}, {'logical'},      {false}, {'Shifts output k-space such that DC occurs at n/2+1 for all dimensions.'}}, ...
         {{'Silent'},                 {''},          {false}, {'logical'},      {false}, {'Turns off all GUI and output messages.  User must supply all necessary arguments.'}}, ...
         {{'CollapseSeg'},            {''},          {false}, {'logical'},      {false}, {'Ignore the SEG dimension, effectively collapsing along that dimension.'}}, ...
         {{'SwapSegAndEco'},          {''},          {false}, {'logical'},      {false}, {'Swap incoming segment and echo indices.  Hack to permit reading of HIFU navigators encoded with SEG indices.'}}, ...
      };

   for iUserInput = 1:length(userInputValues)
      for iField = 1:length(userInputFieldNames)
         UserInputDefault(iUserInput).(userInputFieldNames{iField}{1}) = userInputValues{iUserInput}{iField}{1};
      end     
   end

end

% =============================================================================
function UserInputDefault = derive_user_input_default_fieldname(UserInputDefault)
% -----------------------------------------------------------------------------
% Derives fieldnames to used later for User Input struct.  
% Names based on variation of possible user inputs.
%    ex.  User input 'ReadSomething'  becomes field name 'readSomething'
% -----------------------------------------------------------------------------

   % Setup fieldNames
   for iUserInput = 1:length(UserInputDefault)
      if isempty(UserInputDefault(iUserInput).fieldName)
         fieldName = UserInputDefault(iUserInput).varArgIn;
         fieldName(1) = lower(fieldName(1));
         UserInputDefault(iUserInput).fieldName = fieldName;
      end
   end

end

% =============================================================================
function User = convert_user_input_default_to_struct(UserInputDefault)
% -----------------------------------------------------------------------------
% Set up user input structure with default values.
% -----------------------------------------------------------------------------

   for iUserInput = 1:length(UserInputDefault)
      if strcmp(UserInputDefault(iUserInput).type, 'logical')
         User.(UserInputDefault(iUserInput).fieldName) = UserInputDefault(iUserInput).defaultValue;
      else
         User.(UserInputDefault(iUserInput).fieldName) = uint64(hex2dec(UserInputDefault(iUserInput).defaultValue));
      end
      if strcmp(UserInputDefault(iUserInput).varArgIn(1:4),'Read')
         fieldName = ['found',UserInputDefault(iUserInput).varArgIn(5:end)];
      else
         fieldName = ['found',UserInputDefault(iUserInput).varArgIn];
      end
      User.(fieldName) = false;
   end

end

% =============================================================================
function filenameIn = check_filename_input(filenameIn, nArgIn, Control)
% -----------------------------------------------------------------------------
% Check filenameIn and query if empty.
% -----------------------------------------------------------------------------

   if nArgIn == 0
      if ~Control.User.noGui
         filenameIn = uigetfile('*.dat','Select File to Read');
      else
         if ~Control.User.silent
            textMessage = 'Must provide file name to read';
            display_message(Control.User.noGui,'ERROR',textMessage);
         end
         return
      end
   end

end

% =============================================================================
function filenameInExists = check_filename_exists(filenameIn, Control)
% -----------------------------------------------------------------------------
% Check for file existence.
% -----------------------------------------------------------------------------
   filenameInExists = (exist(filenameIn, 'file') == 2);
   if ~filenameInExists
      if ~Control.User.silent
         textMessage = ['File does not exist: ' filenameIn];
         display_message(Control.User.noGui,'ERROR',textMessage,'Error');
      end
   end
end

% =============================================================================
function Control = modify_with_inputs_outputs(optionsIn, Control, UserInputDefault)
% -----------------------------------------------------------------------------
% Adjust controls based on user inputs/outputs
% -----------------------------------------------------------------------------

   Control.nArgOut = 0;

   iOption = 1;
   while iOption <= length(optionsIn)

      flagFoundMatch = false;
      for iUserInput = 1:length(UserInputDefault)
         flagFoundMatch = flagFoundMatch || strcmpi(optionsIn{iOption},UserInputDefault(iUserInput).varArgIn);
         if flagFoundMatch
            break;
         end
      end
      if flagFoundMatch
         fieldName = UserInputDefault(iUserInput).fieldName;
         fieldNameFound = ['found',UserInputDefault(iUserInput).varArgIn];
         if strcmp(UserInputDefault(iUserInput).type, 'logical')
            Control.User.(fieldName) = true;
            if UserInputDefault(iUserInput).isArgOut
               fieldName = ['iArgOut',UserInputDefault(iUserInput).varArgIn(5:end)];
               Control.nArgOut = Control.nArgOut + 1;
               Control.User.(fieldName) = Control.nArgOut;
            end
         end
         if strcmp(UserInputDefault(iUserInput).type, 'hex2dec')
            iOption = iOption + 1;
            Control.User.(fieldName) = uint64(hex2dec(optionsIn{iOption}));
            Control.User.(fieldNameFound) = true;
         end
      end
      iOption = iOption + 1;
   end % while over options

end

% =============================================================================
function goodToGo = are_inputs_and_outputs_consistent(nOptionsIn, nArgOutMain, Control)
% -----------------------------------------------------------------------------
% Check for consistentcy between argin and argout when a given argin might
%    lead to more data in the output.
% -----------------------------------------------------------------------------

   goodToGo = true;
   if (nOptionsIn > 0) && ((nArgOutMain-1) ~= Control.nArgOut)
      goodToGo = false;
      if ~Control.User.silent
         textMessage = 'You must have an equal number of kSpace output variables for all the kSpace you expect to read.'; 
         display_message(Control.User.noGui,'ERROR',textMessage,'Insufficient Outputs');
      end
   end

   if Control.User.shiftDCToMatrixCenter && (Control.User.readTimeStamp || Control.User.readPMUTimeStamp || Control.User.readTimeSinceLastRF || Control.User.readSequenceTime)
      goodToGo = false;
      if ~Control.User.silent
         textMessage = 'You cannot turn on ShiftDCToMatrixCenter when reading time stamps';
         display_message(Control.User.noGui,'ERROR',textMessage,'Insufficient Outputs');
      end
   end

end

% =============================================================================
function Control = get_controls_from_user_input(filenameIn, optionsIn, numArgIn, numArgOut)
% -----------------------------------------------------------------------------
% Setup main control structure using user inputs.
% -----------------------------------------------------------------------------

   global mdhBitFlag

   Control.GUI.hWaitBar = []; % Initialize to blank

   nOptionsIn = numArgIn - 1;

   UserInputDefault = define_user_input_default_struct_array;
   UserInputDefault = derive_user_input_default_fieldname(UserInputDefault);

   Control.User = convert_user_input_default_to_struct(UserInputDefault);
   Control = modify_with_inputs_outputs(optionsIn, Control, UserInputDefault);

   Control.File.In.name = check_filename_input(filenameIn, nargin, Control);
   if ~check_filename_exists(Control.File.In.name, Control)
      Control = [];
      return
   end

   if ~are_inputs_and_outputs_consistent(nOptionsIn, numArgOut, Control)
      Control = [];
      return 
   end

   % Speed up insertion operations with a mex
   if Control.User.useMexInsertion
      Control.User.useMexInsertion = (exist('copy_column_a2b_single_complex') == 3);
   end

   % Adjust masks for scanevalinfomask run
   if Control.User.scanEvalInfoMaskOnly
      Control.User.foundKeepBitMask = true;
      Control.User.foundDropBitMask = false;
      Control.User.keepMatchAllBits = false;
      Control.User.keepBitMask = uint64(hex2dec('fffffffffffff'));
      Control.User.dropBitMask = uint64(0);
      Control.User.logicKeepOnly = true;
      Control.User.logicDropOnly = false;
   end

   % Address user provided keep/drop flags
   if Control.User.foundDropBitMask || Control.User.foundKeepBitMask
      if Control.User.foundDropBitMask && Control.User.foundKeepBitMask
         mdhBitFlag.bDropEIM.f = Control.User.dropBitMask;
         mdhBitFlag.bKeepEIM.f = Control.User.keepBitMask;
      elseif Control.User.foundDropBitMask
         mdhBitFlag.bDropEIM.f = Control.User.dropBitMask;
      else 
         mdhBitFlag.bKeepEIM.f = Control.User.keepBitMask;
      end
   end

   % Check logic flags
   if ((Control.User.logicKeepOnly + Control.User.logicDropOnly) > 1)
      Control = [];
      return 
   end
   if Control.User.keepMatchAllBits
      mdhBitFlag.bKeepEIM.l = 'bitandeq';
   else
      mdhBitFlag.bKeepEIM.l = 'bitandnz';
   end
   if Control.User.dropMatchAllBits
      mdhBitFlag.bDropEIM.l = 'bitandeq';
   else
      mdhBitFlag.bDropEIM.l = 'bitandnz';
   end
   if Control.User.logicKeepOnly || Control.User.logicDropOnly
      if Control.User.logicKeepOnly
         mdhBitFlag.bDropEIM.l = 'alwaysfalse';
      else % Control.User.logicDropOnly
         mdhBitFlag.bKeepEIM.l = 'alwaystrue';
      end
   end
   if ~Control.User.silent
      print_setup_status(Control.File.In.name, UserInputDefault, Control);
   end

end

% =============================================================================
function [isVE11, MrParcRaidFileHeader, MrParcRaidFileCell] = is_ve11_file(Control, fid)
% -----------------------------------------------------------------------------
% Check that this is a VE11 file.
% -----------------------------------------------------------------------------

   MrParcRaidFileHeader = [];
   MrParcRaidFileCell = {};
   isVE11 = 0;
   FileHeader.totalSize  = fread(fid,1,'uint32=>uint32');
   FileHeader.nMeas = fread(fid,1,'uint32=>uint32');
   if not(and(FileHeader.totalSize < 10000, FileHeader.nMeas <= 64))
      if ~Control.User.silent
         display_message(Control.User.noGui,'ERROR',['Appears to be VB data: ' Control.File.In.name],'Error');
      end
      return 
   else
      MrParcRaidFileHeader = FileHeader;
      isVE11 = 1;
   end

   % Extract measurement parameters
   for iMeas = 1:MrParcRaidFileHeader.nMeas

      MrParcRaidFileCell{iMeas}.measID = fread(fid,1,'uint32=>uint32');
      MrParcRaidFileCell{iMeas}.fileID = fread(fid,1,'uint32=>uint32');
      MrParcRaidFileCell{iMeas}.offsetInFile = fread(fid,1,'uint64=>uint64');
      MrParcRaidFileCell{iMeas}.totalSize = fread(fid,1,'uint64=>uint64');
      MrParcRaidFileCell{iMeas}.patName = fread(fid,64,'uchar=>char');
      MrParcRaidFileCell{iMeas}.protName = fread(fid,64,'uchar=>char');
      MrParcRaidFileCell{iMeas}.Control = Control;

   end

end

% =============================================================================
function [Control, MrParcRaidFileCell, Dim, CoilSelectMap, TextHeader, ...
          CoilSelectMeasMap, RelativeSliceNumber, WiPMemBlock, SliceArray, MeasYapsAsc] = ...
             parse_meas_text_header(Control,fid,MrParcRaidFileCell,iMeas)
% -----------------------------------------------------------------------------
% Parse the text portion of the measurement to extract required dimensions,
%    scaling factors, slice ordering, etc.
% -----------------------------------------------------------------------------

   % Skip to start of this measurement
   fseek(fid,MrParcRaidFileCell{iMeas}.offsetInFile,'bof');

   % TextHeader extraction
   Dim.dMeasHeaderSize = fread(fid,1,'int32');
   Dim.ullMeasHeaderSize = uint64(Dim.dMeasHeaderSize);
   textHeader = fread(fid,Dim.dMeasHeaderSize-4,'uchar=>char');
   textHeader = textHeader';
   TextHeader.complete = textHeader;
   clear textHeader;

   findThisStart = 'MeasYaps';
   findThisEnd = 'Phoenix';
   pStart = strfind(TextHeader.complete,findThisStart) + length(findThisStart) + 5;
   pEnd = strfind(TextHeader.complete,findThisEnd) - 3;
   % Vida has >1 occurence of Phoenix which must now be differentiated
   if length(pEnd) > 1
      indexGood = find((pEnd-pStart) > 0);
      pEnd = pEnd(indexGood(1));
   end
   TextHeader.measYapsAscConv = TextHeader.complete(pStart:pEnd);

   %             {{'address in Siemens XProt',       'which header',   'regular expression to find value',                                'conversion', default}} 
   searchCells = {{'Dim.Raw.readoutOversampleFactor' 'complete'        '.*flReadoutOSFactor[^\s]*\s*{\s*<Precision>\s*[0-9]*\s*([0-9.]+)' 'str2num'    2.0} ...
                  {'Dim.Raw.nCha'                    'measYapsAscConv' '\[0\]\.asList\[([0-9]+)\]\.lRxChannelConnected\s*=\s*'            'length'        } ...
                  {'Dim.Recon.nFourierPartitions'    'complete'        '.*iNoOfFourierPartitions[^\s]*\s*{\s*([0-9.]+)'                   'str2num'      1} ...
                  {'Dim.Recon.phaseEncodeFTLength'   'complete'        '.*iPEFTLength[^\s]*\s*{\s*([0-9.]+)'                              'str2num'       } ...
                  {'Dim.Recon.partitionFTLength'     'complete'        '.*i3DFTLength[^\s]*\s*{\s*([0-9.]+)'                              'str2num'       } ...
                  {'Control.File.scanDimension'      'measYapsAscConv' 'sKSpace.ucDimension[^=]+=\s*([0-9]+)\s*'                          'str2num'       } ...
                  {'Dim.Recon.nLin'                  'complete'        'iPEFTLen[^{}]+{\s*([0-9]+)\s*}'                                   'str2double'    } ...
                  {'Control.File.measUIDString'      'complete'        'MeasUID[^{}]+{\s*([0-9]+)\s*}'                                    'string'        } ...
                  {'Dim.Recon.nPar'                  'measYapsAscConv' 'sKSpace.lPartitions[^=]+=\s*([0-9]+)\s*'                          'str2num'       } ...
                  {'Dim.Recon.nSlc'                  'measYapsAscConv' 'sSliceArray.lSize[^=]+=\s*([0-9]+)\s*'                            'str2num'       } ...
                  {'Dim.Raw.nEco'                    'measYapsAscConv' 'lContrasts\s*=\s*([0-9]+)\s*'                                     'str2num'      1} ...
                  {'Dim.Raw.nSet'                    'measYapsAscConv' 'lSets\s*=\s*([0-9]+)\s*'                                          'str2num'      1} ...
                  {'Dim.Raw.nAcq'                    'measYapsAscConv' 'lAverages\s*=\s*([0-9]+)\s*'                                      'str2num'      1} ...
                  {'Dim.Raw.nRep'                    'measYapsAscConv' 'lRepetitions\s*=\s*([0-9]+)\s*'                                   'str2num'      1} ...
                  {'Dim.Raw.nPhs'                    'measYapsAscConv' 'sPhysioImaging.lPhases\s*=\s*([0-9]+)\s*'                         'str2num'      1} ...
                  {'Control.File.epiFactor'          'measYapsAscConv' 'sFastImaging.lEPIFactor\s*=\s*([0-9]+)\s*'                        'str2num'       } ...
                  {'Control.File.turboFactor'        'measYapsAscConv' 'sFastImaging.lTurboFactor\s*=\s*([0-9]+)\s*'                      'str2num'      1} ...
                  {'Dim.Raw.nMaxRxChannels'          'complete'        '.*iMaxNoOfRxChannels[^\s]*\s*{\s*([0-9]+)'                        'str2num'    128} ...
                  {'Dim.Raw.nPhaseCor'               'complete'        '.*lNoOfPhaseCorrScans[^\s]*\s*{\s*([0-9]+)'                       'str2num'      0} ...
                  {'Dim.Recon.nFourierColumns'       'complete'        '.*iNoOfFourierColumns[^\s]*\s*{\s*([0-9]+)'                       'str2num'       } ...
                  {'Dim.Recon.nFourierLines'         'complete'        '.*NoOfFourierLines[^\s]*\s*{\s*([0-9.]+)'                         'str2num'       } ...
                  {'Dim.Recon.readoutFTLength'       'complete'        '.*iRoFTLength[^\s]*\s*{\s*([0-9.]+)'                              'str2num'       }};
   for searchCell = searchCells
      flagHasDefault = length(searchCell{1}) > 4;
      address = strsplit(searchCell{1}{1},'.');
      searchField = searchCell{1}{2};
      searchString = searchCell{1}{3};
      conversion = searchCell{1}{4};
      if flagHasDefault
         defaultValue = searchCell{1}{5};
      end
      searchResult= regexp(TextHeader.(searchField),searchString,'tokens');
      if ~isempty(searchResult)
         if strcmp(conversion,'str2num')
            searchValue = str2num(searchResult{1}{1});
         elseif strcmp(conversion,'string')
            searchValue = searchResult{1}{1};
         elseif strcmp(conversion,'length')
            searchValue = length(searchResult);
         elseif strcmp(conversion,'str2double')
            searchValue = str2double(searchResult{1}{1});
         else
            Control = [];
            MrParcRaidFileCell = [];
            Dim = [];
            CoilSelectMap = [];
            CoilSelectMeasMap = {};
            RelativeSliceNumber = [];
            WiPMemBlock = [];
            SliceArray = [];
            MeasYapsAsc = [];
          return
         end
      elseif flagHasDefault
         searchValue = defaultValue;
      else
         continue;
      end
      if strcmp(address{1},'Dim')
         Dim.(address{2}).(address{3}) = searchValue;
      end
      if strcmp(address{1},'Control')
         MrParcRaidFileCell{iMeas}.Control.(address{2}).(address{3}) = searchValue;
      end
   end

   if MrParcRaidFileCell{iMeas}.Control.File.scanDimension == 4
      MrParcRaidFileCell{iMeas}.Control.File.is3D = true;
   else
      MrParcRaidFileCell{iMeas}.Control.File.is3D = false;
   end

%   findThis = 'AdjustSeq%/AdjCoilSensSeq';
%   p = strfind(TextHeader.measYapsAscConv,findThis) + length(findThis);
%   if Dim.Raw.nCha > 1 && ~isempty(p)
%      Dim.Raw.nCha = Dim.Raw.nCha-1;
%   end

   if Dim.Raw.nMaxRxChannels < Dim.Raw.nCha
      Dim.Raw.nCha = Dim.Raw.nMaxRxChannels;
   end

   if MrParcRaidFileCell{iMeas}.Control.File.turboFactor > 1
      Dim.Raw.nPhaseCor = 1;
   end
   if MrParcRaidFileCell{iMeas}.Control.File.epiFactor > 1
      Dim.Raw.nPhaseCor = 1;
   end

   Dim.Recon.nNonOSColumns = round(Dim.Recon.nFourierColumns/Dim.Raw.readoutOversampleFactor);

   if Dim.Recon.nFourierPartitions == 1
      Dim.Recon.partitionFTLength = 1;
   end

   %% Raw Data Correction Factors (NOTE: CURRENTLY false!!!!)
   clear CoilSelectMap;
   findThis = '{\s*{\s*{\s*"[^"]+"[^\n]+';
   % Vida started introducing extra carriage returns which requires fancier processing
   pStart = regexp(TextHeader.complete, findThis, 'start');
   %allCoilsInHeaderCell = regexp(TextHeader.complete, findThis, 'match');
   %if ~isempty(allCoilsInHeaderCell)
   if length(pStart) > 0
      pEnd = regexp(TextHeader.complete(pStart(1):end), '}[\n\s]*}[\n\s]*}', 'end');
      if length(pEnd) > 0
         %allCoilsInHeaderCell = allCoilsInHeaderCell{1};
         allCoilsInHeaderCell = TextHeader.complete(pStart(1):pStart(1)+pEnd(1)+1);
         findThis = '{\s*{\s*"(?<name>[^"]+)"\s*}\s*{\s*(?<fft>[\d\.]+)\s*}\s*{\s*(?<re>[\d\.-]+)\s*}\s*{\s*(?<im>[\d\.-]+)\s*}\s*}';
         CoilStructArray = regexp(allCoilsInHeaderCell, findThis,'names');
         if length(CoilStructArray) == Dim.Raw.nCha
            for c=1:Dim.Raw.nCha
               tName = CoilStructArray(c).name;
               CoilSelect.fftScale = sscanf(CoilStructArray(c).fft,'%f');
               CoilSelect.rawDataCorrectionFactor = complex(sscanf(CoilStructArray(c).re,'%f'),sscanf(CoilStructArray(c).im,'%f'));
               CoilSelect.txtOrder = c-1;
               CoilSelectMap.(tName) = CoilSelect;
               clear CoilSelect
            end
         end
         if (length(fieldnames(CoilSelectMap)) ~= Dim.Raw.nCha) && Control.User.applyRawDataCorrection
            if ~Control.User.silent
               display_message(Control.User.noGui,'ERROR','Non-unique channel names in CoilSelect','Error');
            end
            fclose(fid);
            Control = [];
            MrParcRaidFileCell = [];
            Dim = [];
            CoilSelectMap = [];
            CoilSelectMeasMap = {};
            RelativeSliceNumber = [];
            WiPMemBlock = [];
            SliceArray = [];
            MeasYapsAsc = [];
            return
         end
      else
         keyboard
      end
   end

   RelativeSliceNumber = [];
   if Control.User.readRelSliceNumber

      findThis = sprintf('ParamLong."relSliceNumber">[\\s\\n]*{');
      pStart = regexp(TextHeader.complete,findThis);
      if ~isempty(pStart)
         pStart = pStart(1) + strfind(TextHeader.complete(pStart(1):pStart(1) + length(findThis)+20),'{');
         pStart = pStart(1);
         pEnd = strfind(TextHeader.complete(pStart(1):end),'}');
         if ~isempty(pEnd)
            pEnd = pStart + pEnd(1) - 2;
            findThis = '([^0-9\s]*[\-0-9]+)(\s+)';
            relSliceNumCell = regexp(TextHeader.complete(pStart:pEnd),findThis,'tokens');
            if ~isempty(relSliceNumCell)
               relativeSliceNumber = zeros(1,length(relSliceNumCell));
               for i=1:length(relSliceNumCell)
                  relativeSliceNumber(i) = str2double(relSliceNumCell{i}{1});
               end
               RelativeSliceNumber{iMeas} = relativeSliceNumber;
               Control.User.foundRelSliceNumber = true;
            end
         end
      end
   end

   % Meas Yaps AscConv parsing
   %MeasYapsAsc = parse_meas_yaps_ascconv(TextHeader.measYapsAscConv);
   MeasYapsAsc = [];

   % SliceArray
   if Control.User.readSliceArray
      findThis = sprintf('sSliceArray\\.lSize\\s*=\\s*([0-9]+)\\s*');
      lSizeCell = regexp(TextHeader.measYapsAscConv, findThis, 'tokens');
      if ~isempty(lSizeCell)
         Control.User.foundSliceArray = true;
         nSliceLocal = str2num(lSizeCell{1}{1});
         SliceArray.Thickness = zeros(1,nSliceLocal);
         SliceArray.PhaseFOV = zeros(1,nSliceLocal);
         SliceArray.ReadoutFOV = zeros(1,nSliceLocal);
         SliceArray.Position = zeros(3,nSliceLocal);
         SliceArray.Normal = zeros(3,nSliceLocal);
         DirectionStr = {{'dSag'}, {'dCor'}, {'dTra'}};
         for iSlice=1:nSliceLocal
            findThis = sprintf('sSliceArray\\.asSlice\\[%d\\]\\.dThickness\\s*=\\s*([.0-9]+)\\s*', iSlice-1);
            dThicknessCell = regexp(TextHeader.measYapsAscConv, findThis, 'tokens');
            if ~isempty(dThicknessCell)
               SliceArray.Thickness(iSlice) = str2num(dThicknessCell{1}{1});
            end 
            findThis = sprintf('sSliceArray\\.asSlice\\[%d\\]\\.dPhaseFOV\\s*=\\s*([.0-9]+)\\s*', iSlice-1);
            dPhaseFOVCell = regexp(TextHeader.measYapsAscConv, findThis, 'tokens');
            if ~isempty(dPhaseFOVCell)
               SliceArray.PhaseFOV(iSlice) = str2num(dPhaseFOVCell{1}{1});
            end 
            findThis = sprintf('sSliceArray\\.asSlice\\[%d\\]\\.dReadoutFOV\\s*=\\s*([.0-9]+)\\s*', iSlice-1);
            dReadoutFOVCell = regexp(TextHeader.measYapsAscConv, findThis, 'tokens');
            if ~isempty(dReadoutFOVCell)
               SliceArray.ReadoutFOV(iSlice) = str2num(dReadoutFOVCell{1}{1});
            end 
            for iDir=1:3
               findThis = sprintf('sSliceArray\\.asSlice\\[%d\\]\\.sPosition\\.%s\\s*=\\s*([^\\s]+)\\s*',iSlice,DirectionStr{iDir}{1});
               dPosition = regexp(TextHeader.measYapsAscConv, findThis, 'tokens');
               if ~isempty(dPosition)
                  SliceArray.Position(iDir,iSlice) = str2num(dPosition{1}{1});
               end
               findThis = sprintf('sSliceArray\\.asSlice\\[%d\\]\\.sNormal\\.%s\\s*=\\s*([^\\s]+)\\s*',iSlice,DirectionStr{iDir}{1});
               dNormal = regexp(TextHeader.measYapsAscConv, findThis, 'tokens');
               if ~isempty(dNormal)
                  SliceArray.Normal(iDir,iSlice) = str2num(dNormal{1}{1});
               end
            end
         end
      end
   else
      SliceArray = [];
   end

   % WiPMemBlock
   if Control.User.readWiPMemBlock

      findThis = sprintf('ParamMap."sWiPMemBlock">[\\s\\n]*{');
      pStartWiPMemBlock = regexp(TextHeader.complete, findThis, 'ignorecase');
      WiPMemBlock.adRes = [];
      WiPMemBlock.tFree = [];
      WiPMemBlock.adFree = [];
      WiPMemBlock.lCSatBW = [];
      WiPMemBlock.alFree = [];
      if ~isempty(pStartWiPMemBlock)
         pStartWiPMemBlock = pStartWiPMemBlock(1) + strfind(TextHeader.complete(pStartWiPMemBlock(1):pStartWiPMemBlock(1) + length(findThis)+20),'{');
         pStartWiPMemBlock = pStartWiPMemBlock(1);
         findThis = sprintf('<ParamLong."alFree">[\\s\\n]*{');
         pStartAlFree = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
         if ~isempty(pStartAlFree)
            pStartAlFree = pStartWiPMemBlock + pStartAlFree(1); 
            pStartAlFree = pStartAlFree + strfind(TextHeader.complete(pStartAlFree(1):pStartAlFree(1) + length(findThis)+20),'{');
            pStartAlFree = pStartAlFree(1);
            pEnd_alFree = strfind(TextHeader.complete(pStartAlFree(1):end),'}');
            if ~isempty(pEnd_alFree)
               pEnd_alFree = pStartAlFree + pEnd_alFree(1) - 2;
               % Skip any precision/default stuff
               pLessThan = strfind(TextHeader.complete(pStartAlFree:pEnd_alFree),'<');
               if ~isempty(pLessThan)
                  pStartAlFree = pLessThan(end) + pStartAlFree;
                  pStartAlFree = pStartAlFree + strfind(TextHeader.complete(pStartAlFree:pEnd_alFree),sprintf('\n'));
                  pStartAlFree = pStartAlFree(1);
               end
               alFreeCell = regexp(TextHeader.complete(pStartAlFree:pEnd_alFree),'([\-0-9]+)\s+','tokens');
               if ~isempty(alFreeCell)
                  Control.User.foundWiPMemBlock = true;
                  WiPMemBlock.alFree = zeros(1,length(alFreeCell),'int32');
                  for i_alFree=1:length(alFreeCell)
                     WiPMemBlock.alFree(i_alFree) = int32(str2num(alFreeCell{i_alFree}{1}));
                  end
               end
            end
         end
         findThis = sprintf('<ParamLong."lCSatBW">[\\s\\n]*{');
         pStartLCSatBW = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
         if ~isempty(pStartLCSatBW)
            pStartLCSatBW = pStartWiPMemBlock + pStartLCSatBW(1); 
            pStartLCSatBW = pStartLCSatBW + strfind(TextHeader.complete(pStartLCSatBW(1):pStartLCSatBW(1) + length(findThis)+20),'{');
            pStartLCSatBW = pStartLCSatBW(1);
            pEnd_lCSatBW = strfind(TextHeader.complete(pStartLCSatBW(1):end),'}');
            if ~isempty(pEnd_lCSatBW)
               pEnd_lCSatBW = pStartLCSatBW + pEnd_lCSatBW(1) - 2;
               % Skip any precision/default stuff
               pLessThan = strfind(TextHeader.complete(pStartLCSatBW:pEnd_lCSatBW),'<');
               if ~isempty(pLessThan)
                  pStartLCSatBW = pLessThan(end) + pStartLCSatBW;
                  pStartLCSatBW = pStartLCSatBW + strfind(TextHeader.complete(pStartLCSatBW:pEnd_lCSatBW),sprintf('\n'));
                  pStartLCSatBW = pStartLCSatBW(1);
               end
               lCSatBWCell = regexp(TextHeader.complete(pStartLCSatBW:pEnd_lCSatBW),'([\-0-9]+)\s+','tokens');
               if ~isempty(lCSatBWCell)
                  WiPMemBlock.lCSatBW = zeros(1,length(lCSatBWCell),'int32');
                  Control.User.foundWiPMemBlock = true;
                  for i_lCSatBW=1:length(lCSatBWCell)
                     WiPMemBlock.lCSatBW(i_lCSatBW) = int32(str2num(lCSatBWCell{i_lCSatBW}{1}));
                  end
               end
            end
         end
         findThis = sprintf('<ParamDouble."adFree">[\\s\\n]*{');
         pStartAdFree = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
         if ~isempty(pStartAdFree)
            pStartAdFree = pStartWiPMemBlock + pStartAdFree(1); 
            pStartAdFree = pStartAdFree + strfind(TextHeader.complete(pStartAdFree(1):pStartAdFree(1) + length(findThis)+20),'{');
            pStartAdFree = pStartAdFree(1);
            pEnd_adFree = strfind(TextHeader.complete(pStartAdFree(1):end),'}');
            if ~isempty(pEnd_adFree)
               pEnd_adFree = pStartAdFree + pEnd_adFree(1) - 2;
               % Skip any precision/default stuff
               pLessThan = strfind(TextHeader.complete(pStartAdFree:pEnd_adFree),'<');
               if ~isempty(pLessThan)
                  pStartAdFree = pLessThan(end) + pStartAdFree;
                  pStartAdFree = pStartAdFree + strfind(TextHeader.complete(pStartAdFree:pEnd_adFree),sprintf('\n'));
                  pStartAdFree = pStartAdFree(1);
               end
               adFreeCell = regexp(TextHeader.complete(pStartAdFree:pEnd_adFree),'([\-0-9\.]+)\s+','tokens');
               if ~isempty(adFreeCell)
                  Control.User.foundWiPMemBlock = true;
                  WiPMemBlock.adFree = zeros(1,length(adFreeCell));
                  for i_adFree=1:length(adFreeCell)
                     WiPMemBlock.adFree(i_adFree) = str2num(adFreeCell{i_adFree}{1});
                  end
               end
            end
         end
         findThis = sprintf('<ParamDouble."adRes">[\\s\\n]*{');
         pStartAdRes = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
         if ~isempty(pStartAdRes)
            pStartAdRes = pStartWiPMemBlock + pStartAdRes(1); 
            pStartAdRes = pStartAdRes + strfind(TextHeader.complete(pStartAdRes(1):pStartAdRes(1) + length(findThis)+20),'{');
            pStartAdRes = pStartAdRes(1);
            pEnd_adRes = strfind(TextHeader.complete(pStartAdRes(1):end),'}');
            if ~isempty(pEnd_adRes)
               pEnd_adRes = pStartAdRes + pEnd_adRes(1) - 2;
               % Skip any precision/default stuff
               pLessThan = strfind(TextHeader.complete(pStartAdRes:pEnd_adRes),'<');
               if ~isempty(pLessThan)
                  pStartAdRes = pLessThan(end) + pStartAdRes;
                  pStartAdRes = pStartAdRes + strfind(TextHeader.complete(pStartAdRes:pEnd_adRes),sprintf('\n'));
                  pStartAdRes = pStartAdRes(1);
               end
               adResCell = regexp(TextHeader.complete(pStartAdRes:pEnd_adRes),'([\-0-9\.]+)\s+','tokens');
               if ~isempty(adResCell)
                  Control.User.foundWiPMemBlock = true;
                  WiPMemBlock.adRes = zeros(1,length(adResCell));
                  for i_adRes=1:length(adResCell)
                     WiPMemBlock.adRes(i_adRes) = str2num(adResCell{i_adRes}{1});
                  end
               end
            end
         end
         findThis = sprintf('<ParamString."tFree">[\\s\\n]*{');
         pStartTFree = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
         if ~isempty(pStartTFree)
            pStartTFree = pStartWiPMemBlock + pStartTFree(1); 
            pStartTFree = pStartTFree + strfind(TextHeader.complete(pStartTFree(1):pStartTFree(1) + length(findThis)+20),'{');
            pStartTFree = pStartTFree(1);
            pEnd_tFree = strfind(TextHeader.complete(pStartTFree(1):end),'}');
            if ~isempty(pEnd_tFree)
               pEnd_tFree = pStartTFree + pEnd_tFree(1) - 2;
               % Skip any precision/default stuff
               pLessThan = strfind(TextHeader.complete(pStartTFree:pEnd_tFree),'<');
               if ~isempty(pLessThan)
                  pStartTFree = pLessThan(end) + pStartTFree;
                  pStartTFree = pStartTFree + strfind(TextHeader.complete(pStartTFree:pEnd_tFree),sprintf('\n'));
                  pStartTFree = pStartTFree(1);
               end
               tFreeCell = regexp(TextHeader.complete(pStartTFree:pEnd_tFree),'"([^"])"\s+','tokens');
               if ~isempty(tFreeCell)
                  Control.User.foundWiPMemBlock = true;
                  WiPMemBlock.tFree = zeros(1,length(tFreeCell));
                  for i_tFree=1:length(tFreeCell)
                     WiPMemBlock.tFree(i_tFree) = str2num(tFreeCell{i_tFree}{1});
                  end
               end
            end
         end
      end % if wipmemblock found
   else
       WiPMemBlock = [];
   end

   %% FFT Correction Factors
   CoilSelectMeasMap = {};
   for c=1:Dim.Raw.nCha
      CoilSelectMeas.textOrder = c-1;
      findThis = 'asList\[';
      findThis = sprintf('%s%d',findThis,c-1);
      findThis = [findThis,'\]\.sCoilElementID.tCoilID\s*=\s*"([^"]+)"\s*'];
      p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
      CoilSelectMeas.tCoilID = p{1}{1};
      findThis = 'asList\[';
      findThis = sprintf('%s%d',findThis,c-1);
      findThis = [findThis,'\]\.sCoilElementID.tElement\s*=\s*"([^"]+)"\s*'];
      p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
      CoilSelectMeas.tElement = p{1}{1};
      findThis = 'asList\[';
      findThis = sprintf('%s%d',findThis,c-1);
      findThis = [findThis,'\]\.lRxChannelConnected\s*=\s*([0-9.]+)\s*'];
      p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
      CoilSelectMeas.lRxChannel = str2num(p{1}{1});
      findThis = 'asList\[';
      findThis = sprintf('%s%d',findThis,c-1);
      findThis = [findThis,'\]\.lADCChannelConnected\s*=\s*([0-9.]+)\s*'];
      p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
      lADCChannel = str2num(p{1}{1});
      if MrParcRaidFileCell{iMeas}.Control.User.resetFFTScale
         CoilSelectMeas.flFFTCorrectionFactor = single(1.0); 
      else
         findThis = 'aFFT_SCALE\[';
         findThis = sprintf('%s%d',findThis,c-1);
         findThis = [findThis,'\]\.flFactor\s*=\s*([0-9.]+)\s*'];
         p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
         CoilSelectMeas.flFFTCorrectionFactor = single(str2double(p{1}{1}));
      end
      %keySet = sprintf('k%06d',lADCChannel);
      %CoilSelectMeasMap.(keySet) = CoilSelectMeas;
      CoilSelectMeasMap{lADCChannel} = CoilSelectMeas;
      clear CoilSelectMeas;
   end

   %clear keySet valueSet TextHeader.measYapsAscConv;
   clear keySet valueSet;

   if Dim.Recon.nFourierLines == Dim.Recon.nLin && Dim.Recon.nFourierPartitions == Dim.Recon.nPar
      MrParcRaidFileCell{iMeas}.Control.User.noFillToFullFourier = true;
   end

   if ~MrParcRaidFileCell{iMeas}.Control.File.is3D
      Dim.Recon.nPar = Dim.Recon.nFourierPartitions;
   end

   %clear TextHeader; % TextHeader is now returned

   if Dim.Recon.nFourierLines > Dim.Recon.nLin
       Dim.Recon.nLin = Dim.Recon.nFourierLines;
   end

   if Dim.Recon.nFourierPartitions > Dim.Recon.nPar
       Dim.Recon.nPar = Dim.Recon.nFourierPartitions;
   end

end

% =============================================================================
function textMessage = print_setup_status(filenameIn, UserInputDefault, Control)
% -----------------------------------------------------------------------------
% Send text to the user based on GUI or print statement.
% -----------------------------------------------------------------------------
   textMessage = sprintf('Filename=%s', filenameIn);
   for iArg = 1:length(UserInputDefault)
      if strcmp(UserInputDefault(iArg).type, 'logical')
         textMessage = sprintf('%s\n%s=%d', ...
                               textMessage, ...
                               UserInputDefault(iArg).varArgIn, ...
                               Control.User.(UserInputDefault(iArg).fieldName));
      else
         textMessage = sprintf('%s\n%s=''%08x'' (hex string)', ...
                               textMessage, ...
                               UserInputDefault(iArg).varArgIn, ...
                               Control.User.(UserInputDefault(iArg).fieldName));
      end
   end
   display_message(Control.User.noGui,'WAIT',textMessage,'File Will Be Read With These Parameters');
end

% =============================================================================
function textMessage = print_available_options(UserInputDefault)
% -----------------------------------------------------------------------------
% Based on default inputs, indicate to user what is possible.
% -----------------------------------------------------------------------------
   iGroup = 1;
   stepGroup = 8;
   for iArg = 1:stepGroup:length(UserInputDefault)
      if iArg == 1
         textMessage = sprintf('[kSpace, Other] = fast_read_ve11(fileName, ''option1'', ''option2'', ...');
         textMessage = sprintf('%s\n\nfileName is the name of a meas.dat to read and',textMessage);
         textMessage = sprintf('%s\n\nOptions are as follows:', textMessage);
      else
         textMessage = sprintf('Options are as follows:');
      end
      for jArg=iArg:min((iArg+stepGroup),length(UserInputDefault))
         textMessage = sprintf('%s\n\n''%s'': (default', ...
                               textMessage, ...
                               UserInputDefault(jArg).varArgIn);
         switch UserInputDefault(jArg).type
            case 'logical'
               if UserInputDefault(jArg).defaultValue
                  textMessage = sprintf('%s true)', textMessage);
               else
                  textMessage = sprintf('%s false)', textMessage);
               end
            case 'hex2dec'
               textMessage = sprintf('%s ''%08x'' (hex string))', textMessage, UserInputDefault(jArg).defaultValue);
            otherwise
               textMessage = sprintf('%s unknown)', textMessage);
         end
         textMessage = sprintf('%s %s', textMessage, UserInputDefault(jArg).description);
      end
      if usejava('jvm') && ~feature('ShowFigureWindows')
         display(textMessage);
      else
         msgbox(textMessage, sprintf('Usage (Page %d)', iGroup));
      end
      iGroup = iGroup + 1;
   end
end

% =============================================================================
function MdhBlockMeta = append_to_mdhblockmeta(iBlock, iBlockMeta, MdhBlock, varargin)
% -----------------------------------------------------------------------------
% Append to existing MdhBlockMeta list.
% -----------------------------------------------------------------------------
   copyField = {'isSyncData', 'isAcqEnd', 'keepScan', 'dmaLength', 'nScan', 'nSamplesInScan', 'nChannelUsed'};
   %if length(varargin) > 0
   if ~isempty(varargin)
      MdhBlockMeta = varargin{1};
   end
   MdhBlockMeta(iBlockMeta).nBlock = 1;
   MdhBlockMeta(iBlockMeta).iBlock(MdhBlockMeta(iBlockMeta).nBlock) = iBlock;
   for iField=1:length(copyField)
      MdhBlockMeta(iBlockMeta).(copyField{iField}) = MdhBlock(iBlock).(copyField{iField});
   end
end

% =============================================================================
function set_global_mdh_parameters(filenameIn)
% -----------------------------------------------------------------------------
% Handle some global parameters.
% -----------------------------------------------------------------------------

   global MDH_SCANSIZE MDH_CHANSIZE MDH_PMUSIZE MAX_BYTES_PER_FREAD
   global mdhColumns mdhBitFlag
   global MDH_MAP_EIM
   global index_keys
   global MDH_HEX2DEC_FFFFFF

   MDH_HEX2DEC_FFFFFF = hex2dec('FFFFFF');

   index_keys = {'indexLIN','indexACQ','indexSLC','indexPAR','indexECO','indexPHS','indexREP','indexSET','indexSEG', 'indexIDA', 'indexIDB', 'indexIDC', 'indexIDD', 'indexIDE'};

   % EvalInfoMask
   MDH_EIM_ACQEND            = 2^uint64(0);  % last scan 
   MDH_EIM_RTFEEDBACK        = 2^uint64(1);  % Realtime feedback scan
   MDH_EIM_HPFEEDBACK        = 2^uint64(2);  % High perfomance feedback scan
   MDH_EIM_ONLINE            = 2^uint64(3);  % processing should be done online
   MDH_EIM_OFFLINE           = 2^uint64(4);  % processing should be done offline
   MDH_EIM_SYNCDATA          = 2^uint64(5);  % readout contains synchroneous data
   MDH_EIM_LASTSCANINCONCAT  = 2^uint64(8);  % Flag for last scan in concatination
   MDH_EIM_RAWDATACORRECTION = 2^uint64(10); % Correct the rawadata with the rawdata correction factor
   MDH_EIM_LASTSCANINMEAS    = 2^uint64(11); % Flag for last scan in measurement
   MDH_EIM_SCANSCALEFACTOR   = 2^uint64(12); % Flag for scan specific additional scale factor
   MDH_EIM_2NDHADAMARPSE     = 2^uint64(13); % 2nd RF exitation of HADAMAR
   MDH_EIM_REFPHASESTABSCAN  = 2^uint64(14); % reference phase stabilization scan
   MDH_EIM_PHASESTABSCAN     = 2^uint64(15); % phase stabilization scan
   MDH_EIM_D3FFT             = 2^uint64(16); % execute 3D FFT
   MDH_EIM_SIGNREV           = 2^uint64(17); % sign reversal
   MDH_EIM_PHASEFFT          = 2^uint64(18); % execute phase fft
   MDH_EIM_SWAPPED           = 2^uint64(19); % swapped phase/readout direction
   MDH_EIM_POSTSHAREDLINE    = 2^uint64(20); % shared line
   MDH_EIM_PHASCOR           = 2^uint64(21); % phase correction data
   MDH_EIM_PATREFSCAN        = 2^uint64(22); % additonal scan for PAT reference line/partition
   MDH_EIM_PATREFANDIMASCAN  = 2^uint64(23); % additonal scan for PAT reference line/partition that is also used as image scan
   MDH_EIM_REFLECT           = 2^uint64(24); % reflect line
   MDH_EIM_NOISEADJSCAN      = 2^uint64(25); % noise adjust scan 
   MDH_EIM_SHARENOW          = 2^uint64(26); % all lines are acquired from the actual and previous e.g. phases
   MDH_EIM_LASTMEASUREDLINE  = 2^uint64(27); % indicates that the current line is the last measured line of all succeeding e.g. phases
   MDH_EIM_FIRSTSCANINSLICE  = 2^uint64(28); % indicates first scan in slice (needed for time stamps)
   MDH_EIM_LASTSCANINSLICE   = 2^uint64(29); % indicates  last scan in slice (needed for time stamps)
   MDH_EIM_TREFFECTIVEBEGIN  = 2^uint64(30); % indicates the begin time stamp for TReff (triggered measurement)
   MDH_EIM_TREFFECTIVEEND    = 2^uint64(31); % indicates the   end time stamp for TReff (triggered measurement)
   MDH_EIM_MDS_REF_POSITION          = 2^uint64(32); % indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
   MDH_EIM_SLC_AVERAGED              = 2^uint64(33); % indicates avveraged slice for slice partial averaging scheme
   MDH_EIM_TAGFLAG1                  = 2^uint64(34); % adjust scan 
   MDH_EIM_CT_NORMALIZE              = 2^uint64(35); % Marks scans used to calcate correction maps for TimCT-Prescan normalize
   MDH_EIM_SCAN_FIRST                = 2^uint64(36); % Marks the first scan of a particar map
   MDH_EIM_SCAN_LAST                 = 2^uint64(37); % Marks the last scan of a particar map
   MDH_EIM_FIRST_SCAN_IN_BLADE       = 2^uint64(40); % Marks the first line of a blade
   MDH_EIM_LAST_SCAN_IN_BLADE        = 2^uint64(41); % Marks the last line of a blade
   MDH_EIM_LAST_BLADE_IN_TR          = 2^uint64(42); % Set for all lines of the last BLADE in each TR interval
   MDH_EIM_PACE                      = 2^uint64(44); % Distinguishes PACE scans from non PACE scans.
   MDH_EIM_RETRO_LASTPHASE           = 2^uint64(45); % Marks the last phase in a heartbeat
   MDH_EIM_RETRO_ENDOFMEAS           = 2^uint64(46); % Marks an ADC at the end of the measurement
   MDH_EIM_RETRO_REPEATTHISHEARTBEAT = 2^uint64(47); % Repeat the current heartbeat when this bit is found
   MDH_EIM_RETRO_REPEATPREVHEARTBEAT = 2^uint64(48); % Repeat the previous heartbeat when this bit is found
   MDH_EIM_RETRO_ABORTSCANNOW        = 2^uint64(49); % Just abort everything
   MDH_EIM_RETRO_LASTHEARTBEAT       = 2^uint64(50); % This adc is from the last heartbeat (a dummy)
   MDH_EIM_RETRO_DUMMYSCAN           = 2^uint64(51); % This adc is just a dummy scan, throw it away
   MDH_EIM_RETRO_ARRDETDISABLED      = 2^uint64(52); % Disable all arrhythmia detection when this bit is found

   % Construct the map
   keySet = {};
   valueSet = uint64([]);
   iKey = 1;
   whosStruct = whos();
   variableNames = {whosStruct.name};
   for iVar=1:length(variableNames)
      if length(variableNames{iVar}) < 8
         continue
      end
      if strcmp(variableNames{iVar}(1:8),'MDH_EIM_')
         keySet{iKey} = variableNames{iVar}(8:end);
         valueSet(iKey) = eval(variableNames{iVar});
         iKey = iKey + 1;
      end
   end
   MDH_MAP_EIM = containers.Map(keySet,valueSet);

   MDH_DEFAULT_EIM_DROP = MDH_EIM_RTFEEDBACK + ...
                         MDH_EIM_HPFEEDBACK + ...
                         MDH_EIM_PATREFSCAN + ...
                         MDH_EIM_NOISEADJSCAN + ...
                         MDH_EIM_PHASCOR + ...
                         MDH_EIM_REFPHASESTABSCAN + ...
                         MDH_EIM_RETRO_DUMMYSCAN; %+ ...
                         %MDH_EIM_RAWDATACORRECTION;
   MDH_DEFAULT_EIM_KEEP = bitcmp(MDH_DEFAULT_EIM_DROP,'uint64');

   % Some default sizes
   MDH_SCANSIZE = uint64(192);
   MDH_CHANSIZE = uint64(32);
   MDH_PMUSIZE = uint64(60);
   MAX_BYTES_PER_FREAD = 2^uint64(27);
   fileSpecs = dir(filenameIn);
   MAX_BYTES_PER_FREAD = min(MAX_BYTES_PER_FREAD,fileSpecs.bytes); % Constrain to file size

   % Columns within the binary scan lines where
   %    particular data is found along with how
   %    it should be converted.
   mdhColumns.scan.dmaLength.c = [1 4];
   mdhColumns.scan.dmaLength.t = 'uint32';
   mdhColumns.scan.timeStamp.c = [13 16];
   mdhColumns.scan.timeStamp.t = 'uint32';
   mdhColumns.scan.pMUTimeStamp.c = [17 20];
   mdhColumns.scan.pMUTimeStamp.t = 'uint32';
   mdhColumns.scan.evalMask.c = [41 48];
   mdhColumns.scan.evalMask.t = 'uint64';
   mdhColumns.scan.nSamplesInScan.c = [49 50];
   mdhColumns.scan.nSamplesInScan.t = 'uint16';
   mdhColumns.scan.nChannelUsed.c = [51 52];
   mdhColumns.scan.nChannelUsed.t = 'uint16';
   mdhColumns.scan.indexLIN.c = [53 54];
   mdhColumns.scan.indexLIN.t = 'uint16';
   mdhColumns.scan.indexACQ.c = [55 56];
   mdhColumns.scan.indexACQ.t = 'uint16';
   mdhColumns.scan.indexSLC.c = [57 58];
   mdhColumns.scan.indexSLC.t = 'uint16';
   mdhColumns.scan.indexPAR.c = [59 60];
   mdhColumns.scan.indexPAR.t = 'uint16';
   mdhColumns.scan.indexECO.c = [61 62];
   mdhColumns.scan.indexECO.t = 'uint16';
   mdhColumns.scan.indexPHS.c = [63 64];
   mdhColumns.scan.indexPHS.t = 'uint16';
   mdhColumns.scan.indexREP.c = [65 66];
   mdhColumns.scan.indexREP.t = 'uint16';
   mdhColumns.scan.indexSET.c = [67 68];
   mdhColumns.scan.indexSET.t = 'uint16';
   mdhColumns.scan.indexSEG.c = [69 70];
   mdhColumns.scan.indexSEG.t = 'uint16';
   mdhColumns.scan.indexIDA.c = [71 72];
   mdhColumns.scan.indexIDA.t = 'uint16';
   mdhColumns.scan.indexIDB.c = [73 74];
   mdhColumns.scan.indexIDB.t = 'uint16';
   mdhColumns.scan.indexIDC.c = [75 76];
   mdhColumns.scan.indexIDC.t = 'uint16';
   mdhColumns.scan.indexIDD.c = [77 78];
   mdhColumns.scan.indexIDD.t = 'uint16';
   mdhColumns.scan.indexIDE.c = [79 80];
   mdhColumns.scan.indexIDE.t = 'uint16';
   mdhColumns.scan.cutOffDataPre.c = [81 82];
   mdhColumns.scan.cutOffDataPre.t = 'uint16';
   mdhColumns.scan.cutOffDataPost.c = [83 84];
   mdhColumns.scan.cutOffDataPost.t = 'uint16';
   mdhColumns.scan.timeSinceLastRF.c = [93 96];
   mdhColumns.scan.timeSinceLastRF.t = 'uint32';
   mdhColumns.scan.kSpaceCentreLineNo.c = [97 98];
   mdhColumns.scan.kSpaceCentreLineNo.t = 'uint16';
   mdhColumns.scan.kSpaceCentrePartitionNo.c = [99 100];
   mdhColumns.scan.kSpaceCentrePartitionNo.t = 'uint16';
   mdhColumns.scan.sVector.c = [101 112];
   mdhColumns.scan.sVector.t = 'single';
   mdhColumns.scan.quaternion.c = [113 128];
   mdhColumns.scan.quaternion.t = 'single';
   mdhColumns.scan.iceProgramPara.c = [129 184];
   mdhColumns.scan.iceProgramPara.t = 'uint16';

   mdhColumns.chan.sizeChanInBytes.c = [1 4];
   mdhColumns.chan.sizeChanInBytes.t = 'uint32';
   mdhColumns.chan.sequenceTime.c = [17 20];
   mdhColumns.chan.sequenceTime.t = 'uint32';
   mdhColumns.chan.channelID.c = [25 26];
   mdhColumns.chan.channelID.t = 'uint16';

   mdhColumns.pmu.timeStamp.c = [1 8];
   mdhColumns.pmu.timeStamp.t = 'uint32';
   mdhColumns.pmu.counter.c = [9 12];
   mdhColumns.pmu.counter.t = 'uint32';
   mdhColumns.pmu.duration.c = [13 16];
   mdhColumns.pmu.duration.t = 'uint32';
   mdhColumns.pmu.data.c = [17 0];
   mdhColumns.pmu.data.t = 'uint32';

   mdhBitFlag.isNoiseAdj.s = 'evalMask';
   mdhBitFlag.isNoiseAdj.f = MDH_EIM_NOISEADJSCAN;
   mdhBitFlag.isNoiseAdj.l = 'bitandeq';
   mdhBitFlag.isPatRefScan.s = 'evalMask';
   mdhBitFlag.isPatRefScan.f = MDH_EIM_PATREFSCAN;
   mdhBitFlag.isPatRefScan.l = 'bitandeq';
   mdhBitFlag.isPhaseCor.s = 'evalMask';
   mdhBitFlag.isPhaseCor.f = MDH_EIM_PHASCOR;
   mdhBitFlag.isPhaseCor.l = 'bitandeq';
   mdhBitFlag.isRTFeedback.s = 'evalMask';
   mdhBitFlag.isRTFeedback.f = MDH_EIM_RTFEEDBACK;
   mdhBitFlag.isRTFeedback.l = 'bitandeq';
   mdhBitFlag.isAcqEnd.s = 'evalMask';
   mdhBitFlag.isAcqEnd.f = MDH_EIM_ACQEND;
   mdhBitFlag.isAcqEnd.l = 'bitandeq';
   mdhBitFlag.isReflect.s = 'evalMask';
   mdhBitFlag.isReflect.f = MDH_EIM_REFLECT;
   mdhBitFlag.isReflect.l = 'bitandeq';
   mdhBitFlag.isRawDataCorrection.s = 'evalMask';
   mdhBitFlag.isRawDataCorrection.f = MDH_EIM_RAWDATACORRECTION;
   mdhBitFlag.isRawDataCorrection.l = 'bitandeq';
   mdhBitFlag.isSyncData.s = 'evalMask';
   mdhBitFlag.isSyncData.f = MDH_EIM_SYNCDATA;
   mdhBitFlag.isSyncData.l = 'bitandeq';
   mdhBitFlag.bKeepEIM.s = 'evalMask';
   mdhBitFlag.bKeepEIM.f = MDH_DEFAULT_EIM_KEEP;
   mdhBitFlag.bKeepEIM.l = 'bitandnz';
   mdhBitFlag.bDropEIM.s = 'evalMask';
   mdhBitFlag.bDropEIM.f = MDH_DEFAULT_EIM_DROP;
   mdhBitFlag.bDropEIM.l = 'bitandnz';

end

% =============================================================================
function display_eval_mask(evalMaskArray,iMeas)
% -----------------------------------------------------------------------------
% Display eval mask array data, converting to human readable form
% -----------------------------------------------------------------------------

   global MDH_MAP_EIM

   textMessage = 'dma, nRead, nChan, eimint, eimhex, eimtext';

   bitsEIM = values(MDH_MAP_EIM);
   namesEIM = keys(MDH_MAP_EIM);
   [~,indexSort] = sort(cell2mat(bitsEIM));
   u64Zero = uint64(0);
   tableData = {};
   for iEval =1:length(evalMaskArray)
      textMessage = sprintf('%s\n%s', textMessage, evalMaskArray{iEval});
      evalMaskData = uint64(sscanf(evalMaskArray{iEval}, '%d, %d, %d, %lu'));
      for iData=1:length(evalMaskData)
         tableData{iEval,iData} = evalMaskData(iData);
      end
      tableData{iEval,length(evalMaskData)+1} = sprintf('0x%013X', evalMaskData(end));
      iMatch = 0;
      for iFlag=1:length(bitsEIM)
         iSort = indexSort(iFlag);
         if bitand(evalMaskData(end),bitsEIM{iSort},'uint64') > u64Zero
            if iMatch < 1 
               namesStr = sprintf('%s(0x%013X)',namesEIM{iSort}, bitsEIM{iSort});
            else
               namesStr = sprintf('%s+%s(0x%013X)', namesStr, namesEIM{iSort}, bitsEIM{iSort});
            end
            iMatch = iMatch + 1;
         end
      end
      if iMatch > 0
         textMessage = sprintf('%s, %s', textMessage, namesStr);
         tableData{iEval,length(evalMaskData)+2} = namesStr;
      else
         tableData{iEval,length(evalMaskData)+2} = '';
      end
   end
   if usejava('jvm') && ~feature('ShowFigureWindows')
      textMessage = sprintf('Unique Eval1Mask Combinations for Meas %d\n%s', iMeas, textMessage);
      display(textMessage);
   else
      uf = figure('Name', sprintf('Unique Eval1Mask Combinations for Meas %d', iMeas), ...
                  'NumberTitle', 'off', ...
                  'MenuBar', 'none', ...
                  'ToolBar', 'none', ...
                  'Position', [200 200 840 270], ...
                  'Visible', 'off');
      t = uitable(uf);
      t.Data = tableData;
      t.ColumnName = {'dma','nRead','nChan','eimval','eimhex','eimText'};
      t.Position = [20, 20, 800, 230];
      t.ColumnWidth = {'auto','auto','auto','auto',100,900};
      set(uf,'Visible','on');
   end

end

% =============================================================================
function MdhParam = extract_single_mdh_parameter(mdhBuffer, paramType, nameToExtract)
% -----------------------------------------------------------------------------
% Extract a parameter from the uint8 array based on input description.
% -----------------------------------------------------------------------------

   global mdhColumns 

   columns = mdhColumns.(paramType).(nameToExtract).c;
   temp = mdhBuffer(columns(1):columns(2),:);
   MdhParam = typecast(temp(:),mdhColumns.(paramType).(nameToExtract).t);

   return

end

% =============================================================================
function MdhParam = extract_mdh_parameters(mdhBuffer, paramType, namesToExtract)
% -----------------------------------------------------------------------------
% Extract an array of parameters, similar to extract_single_mdh_parameter.
% -----------------------------------------------------------------------------

   global mdhColumns 

   switch nargin
      case 3
         fieldNames = namesToExtract;
      case 2
         fieldNames = fieldnames(mdhColumns.(paramType));
      otherwise
         paramType = 'scan';
         fieldNames = fieldnames(mdhColumns.(paramType));
   end

   for iField = 1:length(fieldNames)
      columns = mdhColumns.(paramType).(fieldNames{iField}).c;
      if columns(2) == 0 % PMU data goes to the end of the buffer
         columns(2) = length(mdhBuffer);
      end
      temp = mdhBuffer(columns(1):columns(2),:);
      MdhParam.(fieldNames{iField}) = typecast(temp(:),mdhColumns.(paramType).(fieldNames{iField}).t);
   end

   return

end

% =============================================================================
function MdhBitFlags = extract_mdh_bitflags(mdhBuffer, flagsToExtract)
% -----------------------------------------------------------------------------
% Pull bitflags from the MDH evalMasks
% -----------------------------------------------------------------------------

   global mdhBitFlag
   
   if nargin == 1
      fieldNames = fieldnames(mdhBitFlag);
   else
      fieldNames = flagsToExtract;
   end

   MdhBitFlags.evalMask = extract_single_mdh_parameter(mdhBuffer, 'scan', 'evalMask');
   for iField = 1:length(fieldNames)
      evalMask = mdhBitFlag.(fieldNames{iField}).s;
      bitFlag = mdhBitFlag.(fieldNames{iField}).f;
      logicType = mdhBitFlag.(fieldNames{iField}).l;
      switch logicType
         case 'bitandeq'
            MdhBitFlags.(fieldNames{iField}) = bitand(MdhBitFlags.(evalMask), bitFlag, 'uint64') == bitFlag;
         case 'bitandnz'
            MdhBitFlags.(fieldNames{iField}) = bitand(MdhBitFlags.(evalMask), bitFlag, 'uint64')  > uint64(0);
         case 'alwaystrue'
            MdhBitFlags.(fieldNames{iField}) = ones(size(MdhBitFlags.(evalMask))) > 0;
         case 'alwaysfalse'
            MdhBitFlags.(fieldNames{iField}) = ones(size(MdhBitFlags.(evalMask))) < 1;
         otherwise % same as bitandeq
            MdhBitFlags.(fieldNames{iField}) = bitand(MdhBitFlags.(evalMask), bitFlag, 'uint64') == bitFlag;
      end % switch
   end

   return;

end

% =============================================================================
function MdhBlock = init_mdh_block(offsetInFile, dmaLength, sizeFirstChanInBytes, MdhBitFlags, MdhParam)
% -----------------------------------------------------------------------------
% Set up MdhBlock structure
% -----------------------------------------------------------------------------

   MdhBlock.offsetInFile = offsetInFile;
   MdhBlock.dmaLength = dmaLength;
   MdhBlock.nScan = uint64(0); % Initialize with zero value
   MdhBlock.isSyncData = MdhBitFlags.isSyncData;
   MdhBlock.isAcqEnd = MdhBitFlags.isAcqEnd;
   MdhBlock.nSamplesInScan = uint64(MdhParam.nSamplesInScan);
   MdhBlock.nChannelUsed = MdhParam.nChannelUsed;
   MdhBlock.sizeFirstChanInBytes = sizeFirstChanInBytes;
   MdhBlock.keepScan = MdhBitFlags.bKeepEIM(1) && ~MdhBitFlags.bDropEIM(1);

end

% =============================================================================
function dimOut = max_from_mdh_param(dimIn, nChannelUsed, dimNames, MdhParam, indices)
% -----------------------------------------------------------------------------
% Compute output dimensions
% -----------------------------------------------------------------------------
   dimOut = dimIn;
   dimOut(1) = max(dimIn(1), nChannelUsed-1);
   for iDim=2:numel(dimNames)
      indexName = strcat('index',dimNames{iDim});
      dimOut(iDim) = max(dimIn(iDim), max(MdhParam.(indexName)(indices(:))));
   end
   return
end

% =============================================================================
function rankOut = rank_from_mdh_param(rankIn, dimNames, MdhParam, nChannels, indices)
% -----------------------------------------------------------------------------
% Determine which dimensions vary most rapidly based on mdh values
% -----------------------------------------------------------------------------
   rankOut = rankIn;
   rankOut(1) = rankIn(1) + uint64(numel(indices(2:end))) * uint64(nChannels);
   for iDim=2:numel(dimNames)
      indexName = strcat('index',dimNames{iDim});
      diffMask = MdhParam.(indexName)(indices(2:end)) ~= MdhParam.(indexName)(indices(1:end-1));
      rankOut(iDim) = rankIn(iDim) + uint64(sum(diffMask,1));
   end
   return
end

% =============================================================================
function rankOut = rank_from_index_prev(rankIn, dimNames, MdhParam, nChannels, index, indexPrev)
% -----------------------------------------------------------------------------
% Determine column order based on indices in the mdh
% -----------------------------------------------------------------------------
   rankOut = rankIn;
   rankOut(1) = rankOut(1) + uint64(nChannels);
   for iDim=2:numel(dimNames)
      indexName = strcat('index',dimNames{iDim});
      diffMask = indexPrev.(indexName) ~= MdhParam.(indexName)(index);
      rankOut(iDim) = rankIn(iDim) + uint64(sum(diffMask,1));
   end
   return
end

% =============================================================================
function indexPrev = set_index_prev(dimNames, MdhParam, index)
% -----------------------------------------------------------------------------
% Bookeeping for tracking what dimensions vary most rapidly
% -----------------------------------------------------------------------------
   for iDim=2:numel(dimNames)
      indexName = strcat('index',dimNames{iDim});
      indexPrev.(indexName) = MdhParam.(indexName)(index);
   end
   return
end
 
% =============================================================================
function output = compute_cumulative_product(input)
% -----------------------------------------------------------------------------
% Compute the cumulative product for an input array.  Apparently before matlab
%    2013, cumprod only supports floating cumprod and I need one for integers!
%
% output = compute_cumulative_product(input) returns array of same type and
%    structure as input with cumulative product computed from first to last
%    element
%
% HISTORY:   05 Dec 2014. J. Roberts. roberts@ucair.med.utah.edu
%            UCAIR, Dept Radiology, School of Med, U of U, SLC, UT
% -----------------------------------------------------------------------------

   output = ones(size(input),class(input(1)));
   output(1) = input(1);
   for p=2:numel(input)
      output(p) = output(p-1) * input(p);
   end

end

% =============================================================================
function output = compute_product(input)
% -----------------------------------------------------------------------------
% Compute the product for an input array.  Apparently before matlab 2013, prod
%    only supports floating prod and I need one for integers!
%
% output = compute_product(input) returns array of same type and structure
%                                     as input with product 
%                                     computed from first to last element
%
% HISTORY:   05 Dec 2014. J. Roberts. roberts@ucair.med.utah.edu
%            UCAIR, Dept Radiology, School of Med, U of U, SLC, UT
% -----------------------------------------------------------------------------

   output = input(1);
   for p=2:numel(input)
      output = output * input(p);
   end

end

% =============================================================================
function display_message(noGui,typeMessage,textMessage,textExtra)
% -----------------------------------------------------------------------------

   if noGui
      disp(textMessage);
   else
      if strcmp(typeMessage,'WARNING')
         warndlg(textMessage);
      elseif strcmp(typeMessage,'ERROR')
         errordlg(textMessage,textExtra);
      else
         uiwait(msgbox(textMessage, textExtra, 'modal'));
      end
   end

end

% =============================================================================
function [evalMaskArray, sizeAllScansInBytes, sizeChanInBytes, MdhBlock, nBlockSYNCDATA, ...
          nBlockACQEND, nBlockScan, blockNScan, nScanInMeas, dmaLengthScan, ...
          nSamplesInScan] = ...
             chunk_into_equal_line_length_blocks2(fid,MrParcRaidFileCell, iMeas, Dim, Control)
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
   global MDH_SCANSIZE MAX_BYTES_PER_FREAD
   global MDH_HEX2DEC_FFFFFF

   if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
      mdhParamToExtract = {'dmaLength','nSamplesInScan','nChannelUsed','evalMask'};
   else
      mdhParamToExtract = {'dmaLength','nSamplesInScan','nChannelUsed'};
   end
   evalMaskArray = {};
   sizeAllScansInBytes = MrParcRaidFileCell{iMeas}.totalSize - Dim.ullMeasHeaderSize;
   fseek(fid,MrParcRaidFileCell{iMeas}.offsetInFile + Dim.ullMeasHeaderSize,'bof');
   % Read the first line to get the DMA length and line type
   mdhScan = fread(fid,MDH_SCANSIZE+4,'uchar=>uchar');
   MdhParam = extract_mdh_parameters(mdhScan, 'scan', mdhParamToExtract);
   if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
      evalStr = sprintf('%d, %d, %d, %u', MdhParam.dmaLength, MdhParam.nSamplesInScan, MdhParam.nChannelUsed, MdhParam.evalMask);
      if isempty(find(strcmp(evalMaskArray,evalStr)))
         evalMaskArray{length(evalMaskArray)+1} = evalStr;
      end
   end
   MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData', 'bKeepEIM', 'bDropEIM'});
   bHasACQEND = MdhBitFlags.isAcqEnd;
   bKeepScan = MdhBitFlags.bKeepEIM(1) && ~MdhBitFlags.bDropEIM(1);
   dmaLength = uint64(bitand(MdhParam.dmaLength,MDH_HEX2DEC_FFFFFF)); % Only first 3 bytes used for dmaLength
   sizeFirstChanInBytes = uint64(bitshift(extract_single_mdh_parameter(mdhScan(MDH_SCANSIZE+1:MDH_SCANSIZE+4),'chan','sizeChanInBytes'),-8));
   iRelativePosition = uint64(1);
   iBlock = 1;
   offsetInFile = MrParcRaidFileCell{iMeas}.offsetInFile + Dim.ullMeasHeaderSize;
   % NOTE: We make a very simple assumption that the same number of
   %    samples per scan will hold for scans of the same dmaLength
   % Similarly for the number of channels (which mean nothing in SYNCDATA scans)
   MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, sizeFirstChanInBytes, MdhBitFlags, MdhParam);
   while ((iRelativePosition < sizeAllScansInBytes) && ~bHasACQEND)
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         textMessage = 'Pass 1 of 3: Parsing blocks';
         waitbar(double(iRelativePosition)/double(sizeAllScansInBytes), Control.GUI.hWaitBar, textMessage);
      end
      fseek(fid,MdhBlock(iBlock).offsetInFile + MdhBlock(iBlock).nScan * dmaLength,'bof');
      bytesRemaining = sizeAllScansInBytes-iRelativePosition+uint64(1);
      if dmaLength < 2048
         nScanToRead=uint64(177);  % We don't usually get a lot short lines
      else
         nScanToRead = max(idivide(MAX_BYTES_PER_FREAD, dmaLength,'fix'),uint64(1));
      end   
      %if iBlock > 20
      %   indexKeep = find([MdhBlock.keepScan]);
      %   if length(indexKeep) > 0
      %      nScanToRead = min(nScanToRead, median([MdhBlock(indexKeep).nScan]));
      %   end
      %end
      bytesToRead = nScanToRead*dmaLength;
      if (bytesToRead > bytesRemaining)
         nScanToRead = idivide(bytesRemaining,dmaLength,'fix');
         if nScanToRead < 1
            [mdhScan, readCount] = fread(fid,MDH_SCANSIZE+4,'uchar=>uchar');
            MdhParam = extract_mdh_parameters(mdhScan, 'scan', mdhParamToExtract);
            if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
               evalStr = sprintf('%d, %d, %d, %u', MdhParam.dmaLength, MdhParam.nSamplesInScan, MdhParam.nChannelUsed, MdhParam.evalMask);
               if isempty(find(strcmp(evalMaskArray,evalStr)))
                  evalMaskArray{length(evalMaskArray)+1} = evalStr;
               end
            end
            MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData', 'bKeepEIM', 'bDropEIM'});
            bKeepScanNext = MdhBitFlags.bKeepEIM(1) && ~MdhBitFlags.bDropEIM(1);
            dmaLengthNext = uint64(bitand(MdhParam.dmaLength,MDH_HEX2DEC_FFFFFF)); % Only first 3 bytes used by dmaLength
            if readCount > MDH_SCANSIZE
               sizeFirstChanInBytesNext = uint64(bitshift(extract_single_mdh_parameter(mdhScan(MDH_SCANSIZE+1:MDH_SCANSIZE+4),'chan','sizeChanInBytes'),-8));
            else
               sizeFirstChanInBytesNext = uint64(0);
            end
            nScanToRead = uint64(1);
            if (dmaLengthNext == dmaLength) && (bKeepScan == bKeepScanNext)
               error_message
            else
               iBlock = iBlock + 1;
               offsetInFile = MdhBlock(iBlock-1).offsetInFile + MdhBlock(iBlock-1).nScan * dmaLength;
               dmaLength = dmaLengthNext;
               bKeepScan = bKeepScanNext;
               sizeFirstChanInBytes = sizeFirstChanInBytesNext;
               MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, sizeFirstChanInBytes, MdhBitFlags, MdhParam);
               fseek(fid,MdhBlock(iBlock).offsetInFile + MdhBlock(iBlock).nScan * dmaLength,'bof');
            end
         end
      end
      bytesToRead = nScanToRead*dmaLength;
      [blockOfScans, readCount] = fread(fid,bytesToRead,'uchar=>uchar');
      if MdhBlock(iBlock).isAcqEnd && (readCount < bytesToRead)
         dmaLength = readCount / nScanToRead;
      end
      iStart = int64(1);
      while (uint64(iStart) + MDH_SCANSIZE - 1) <= readCount
         mdhScan = blockOfScans(uint64(iStart):min(uint64(iStart)+MDH_SCANSIZE+4-1,readCount));
         MdhParam = extract_mdh_parameters(mdhScan, 'scan', mdhParamToExtract);
         if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
            evalStr = sprintf('%d, %d, %d, %u', MdhParam.dmaLength, MdhParam.nSamplesInScan, MdhParam.nChannelUsed, MdhParam.evalMask);
            if isempty(find(strcmp(evalMaskArray,evalStr)))
               evalMaskArray{length(evalMaskArray)+1} = evalStr;
            end
         end
         MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData', 'bKeepEIM', 'bDropEIM'});
         bKeepScanNext = MdhBitFlags.bKeepEIM(1) && ~MdhBitFlags.bDropEIM(1);
         dmaLengthNext = uint64(bitand(MdhParam.dmaLength,MDH_HEX2DEC_FFFFFF)); % Only first 3 bytes for dmaLength
         if length(mdhScan) > MDH_SCANSIZE
            sizeFirstChanInBytesNext = uint64(bitshift(extract_single_mdh_parameter(mdhScan(MDH_SCANSIZE+1:MDH_SCANSIZE+4),'chan','sizeChanInBytes'),-8));
         else
            sizeFirstChanInBytesNext = uint64(0);
         end
         bHasACQEND = bHasACQEND || MdhBlock(iBlock).isAcqEnd;
         if (dmaLengthNext == dmaLength) && (bKeepScanNext == bKeepScan)
            MdhBlock(iBlock).nScan = MdhBlock(iBlock).nScan + 1;
            iRelativePosition = iRelativePosition + dmaLength;
         else
            iBlock = iBlock + 1;
            %if (iBlock > 8423)
            %   keyboard
            %end
            offsetInFile =MdhBlock(iBlock-1).offsetInFile + MdhBlock(iBlock-1).nScan * dmaLength;
            flagDMALengthChange = dmaLength ~= dmaLengthNext;
            dmaLength = dmaLengthNext;
            bKeepScan = bKeepScanNext;
            sizeFirstChanInBytes = sizeFirstChanInBytesNext;
            MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, sizeFirstChanInBytes, MdhBitFlags, MdhParam);
            %if flagDMALengthChange % No need to reread a chunk if dma length has not changed
            %   break;
            %end
            iStart = iStart - int64(dmaLength); % backup and re-read this data into a the new block
         end
         if bHasACQEND
            break;
         end
         iStart = iStart + int64(dmaLength);
      end % for over block of scans with same dmalength
   end % while over full file
   if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
      textMessage = 'Pass 1 of 3: Parsing blocks';
      waitbar(double(iRelativePosition)/double(sizeAllScansInBytes), Control.GUI.hWaitBar, textMessage);
   end
   if ((iRelativePosition < sizeAllScansInBytes) && bHasACQEND)         
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'WARNING','ACQ_END before end of file.');
         keyboard
      end
   end
   if ~bHasACQEND
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'WARNING','No ACQ_END encountered.  Data possibly incomplete or corrupted.');
      end
   end
   nBlockSYNCDATA = sum([MdhBlock.isSyncData]);
   nBlockACQEND = sum([MdhBlock.isAcqEnd]);
   indexScan = find(([MdhBlock.isAcqEnd] == 0) & ([MdhBlock.isSyncData] == 0) & ([MdhBlock.keepScan] > 0));
   nBlockScanAll = length(indexScan);
   if nBlockScanAll < 1
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'Warning',sprintf('No readouts detected for meas %d.', iMeas),'Warning');
      end
      %fclose(fid);
      sizeAllScansInBytes = [];
      sizeChanInBytes = [];
      MdhBlock = [];
      nBlockSYNCDATA = [];
      nBlockACQEND = [];
      nBlockScan = [];
      blockNScan = [];
      nScanInMeas = [];
      dmaLengthScan = [];
      nSamplesInScan = [];
      return
   end

   nSamplesInScanUnique = unique([MdhBlock(indexScan).nSamplesInScan]);
   nScanInMeas = uint64(sum([MdhBlock(indexScan).nScan]));
   if (length(nSamplesInScanUnique) == 1)
      nBlockScan = nBlockScanAll;
      blockNScan = [MdhBlock(indexScan).nScan];
      if (length(unique([MdhBlock(indexScan).dmaLength])) ~= 1) && ~MrParcRaidFileCell{iMeas}.Control.User.noKSpace
         if ~MrParcRaidFileCell{iMeas}.Control.User.silent
            if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
               display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'Warning',sprintf('Multiple different dmaLength in meas %d',iMeas),'Warning');
            else
               display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'ERROR',sprintf('Multiple different dmaLength in meas %d',iMeas),'Error');
               keyboard
               fclose(fid);
            end
         end
         sizeAllScansInBytes = [];
         sizeChanInBytes = [];
         MdhBlock = [];
         nBlockSYNCDATA = [];
         nBlockACQEND = [];
         nBlockScan = [];
         blockNScan = [];
         nScanInMeas = [];
         dmaLengthScan = [];
         nSamplesInScan = [];
         return
      end
      dmaLengthScan = MdhBlock(indexScan(1)).dmaLength;
      nSamplesInScan = nSamplesInScanUnique(1);
      sizeChanInBytes = MdhBlock(indexScan(1)).sizeFirstChanInBytes;
   else
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
            display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'Warning',sprintf('Multiple different samples per scan in meas %d',iMeas),'Warning');
         else
            if ~MrParcRaidFileCell{iMeas}.Control.User.noKSpace
               display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'ERROR',sprintf('Multiple different samples per scan in meas %d',iMeas),'Error');
               keyboard
            end
         end
      end
   end
end

% =============================================================================
function [evalMaskArray, sizeAllScansInBytes, sizeChanInBytes, MdhBlock, nBlockSYNCDATA, ...
          nBlockACQEND, nBlockScan, blockNScan, nScanInMeas, dmaLengthScan, ...
          nSamplesInScan] = ...
             chunk_into_equal_line_length_blocks(fid,MrParcRaidFileCell, iMeas, Dim, Control)
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
   global MDH_SCANSIZE MAX_BYTES_PER_FREAD
   global MDH_HEX2DEC_FFFFFF

   if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
      mdhParamToExtract = {'dmaLength','nSamplesInScan','nChannelUsed','evalMask'};
   else
      mdhParamToExtract = {'dmaLength','nSamplesInScan','nChannelUsed'};
   end
   evalMaskArray = {};
   sizeAllScansInBytes = MrParcRaidFileCell{iMeas}.totalSize - Dim.ullMeasHeaderSize;
   fseek(fid,MrParcRaidFileCell{iMeas}.offsetInFile + Dim.ullMeasHeaderSize,'bof');
   % Read the first line to get the DMA length and line type
   mdhScan = fread(fid,MDH_SCANSIZE+4,'uchar=>uchar');
   MdhParam = extract_mdh_parameters(mdhScan, 'scan', mdhParamToExtract);
   if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
      evalStr = sprintf('%d, %d, %d, %u', MdhParam.dmaLength, MdhParam.nSamplesInScan, MdhParam.nChannelUsed, MdhParam.evalMask);
      if isempty(find(strcmp(evalMaskArray,evalStr)))
         evalMaskArray{length(evalMaskArray)+1} = evalStr;
      end
   end
   MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData', 'bKeepEIM', 'bDropEIM'});
   bHasACQEND = MdhBitFlags.isAcqEnd;
   bKeepScan = MdhBitFlags.bKeepEIM(1) && ~MdhBitFlags.bDropEIM(1);
   dmaLength = uint64(bitand(MdhParam.dmaLength,MDH_HEX2DEC_FFFFFF)); % Only first 3 bytes used for dmaLength
   sizeFirstChanInBytes = uint64(bitshift(extract_single_mdh_parameter(mdhScan(MDH_SCANSIZE+1:MDH_SCANSIZE+4),'chan','sizeChanInBytes'),-8));
   iRelativePosition = uint64(1);
   iBlock = 1;
   offsetInFile = MrParcRaidFileCell{iMeas}.offsetInFile + Dim.ullMeasHeaderSize;
   % NOTE: We make a very simple assumption that the same number of
   %    samples per scan will hold for scans of the same dmaLength
   % Similarly for the number of channels (which mean nothing in SYNCDATA scans)
   MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, sizeFirstChanInBytes, MdhBitFlags, MdhParam);
   while ((iRelativePosition < sizeAllScansInBytes) && ~bHasACQEND)
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         textMessage = 'Pass 1 of 3: Parsing blocks';
         waitbar(double(iRelativePosition)/double(sizeAllScansInBytes), Control.GUI.hWaitBar, textMessage);
      end
      fseek(fid,MdhBlock(iBlock).offsetInFile + MdhBlock(iBlock).nScan * dmaLength,'bof');
      bytesRemaining = sizeAllScansInBytes-iRelativePosition+uint64(1);
      if dmaLength < 2048
         nScanToRead=uint64(177);  % We don't usually get a lot short lines
      else
         nScanToRead = max(idivide(MAX_BYTES_PER_FREAD, dmaLength,'fix'),uint64(1));
      end   
      %if iBlock > 20
      %   indexKeep = find([MdhBlock.keepScan]);
      %   if length(indexKeep) > 0
      %      nScanToRead = min(nScanToRead, median([MdhBlock(indexKeep).nScan]));
      %   end
      %end
      bytesToRead = nScanToRead*dmaLength;
      if (bytesToRead > bytesRemaining)
         nScanToRead = idivide(bytesRemaining,dmaLength,'fix');
         if nScanToRead < 1
            [mdhScan, readCount] = fread(fid,MDH_SCANSIZE+4,'uchar=>uchar');
            MdhParam = extract_mdh_parameters(mdhScan, 'scan', mdhParamToExtract);
            if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
               evalStr = sprintf('%d, %d, %d, %u', MdhParam.dmaLength, MdhParam.nSamplesInScan, MdhParam.nChannelUsed, MdhParam.evalMask);
               if isempty(find(strcmp(evalMaskArray,evalStr)))
                  evalMaskArray{length(evalMaskArray)+1} = evalStr;
               end
            end
            MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData', 'bKeepEIM', 'bDropEIM'});
            bKeepScanNext = MdhBitFlags.bKeepEIM(1) && ~MdhBitFlags.bDropEIM(1);
            dmaLengthNext = uint64(bitand(MdhParam.dmaLength,MDH_HEX2DEC_FFFFFF)); % Only first 3 bytes used by dmaLength
            if readCount > MDH_SCANSIZE
               sizeFirstChanInBytesNext = uint64(bitshift(extract_single_mdh_parameter(mdhScan(MDH_SCANSIZE+1:MDH_SCANSIZE+4),'chan','sizeChanInBytes'),-8));
            else
               sizeFirstChanInBytesNext = uint64(0);
            end
            nScanToRead = uint64(1);
            if (dmaLengthNext == dmaLength) && (bKeepScan == bKeepScanNext)
               error_message
            else
               iBlock = iBlock + 1;
               offsetInFile = MdhBlock(iBlock-1).offsetInFile + MdhBlock(iBlock-1).nScan * dmaLength;
               dmaLength = dmaLengthNext;
               bKeepScan = bKeepScanNext;
               sizeFirstChanInBytes = sizeFirstChanInBytesNext;
               MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, sizeFirstChanInBytes, MdhBitFlags, MdhParam);
               fseek(fid,MdhBlock(iBlock).offsetInFile + MdhBlock(iBlock).nScan * dmaLength,'bof');
            end
         end
      end
      bytesToRead = dmaLength*nScanToRead;
      [blockOfScans, readCount] = fread(fid,bytesToRead,'uchar=>uchar');
      if MdhBlock(iBlock).isAcqEnd && (readCount < bytesToRead)
         dmaLength = uint64(readCount / nScanToRead);
      end
      iStart = int64(1);
      while (uint64(iStart) + MDH_SCANSIZE - 1) <= readCount
         mdhScan = blockOfScans(uint64(iStart):min(uint64(iStart)+MDH_SCANSIZE+4-1,readCount));
         MdhParam = extract_mdh_parameters(mdhScan, 'scan', mdhParamToExtract);
         if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
            evalStr = sprintf('%d, %d, %d, %u', MdhParam.dmaLength, MdhParam.nSamplesInScan, MdhParam.nChannelUsed, MdhParam.evalMask);
            if isempty(find(strcmp(evalMaskArray,evalStr)))
               evalMaskArray{length(evalMaskArray)+1} = evalStr;
            end
         end
         MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData', 'bKeepEIM', 'bDropEIM'});
         bKeepScanNext = MdhBitFlags.bKeepEIM(1) && ~MdhBitFlags.bDropEIM(1);
         dmaLengthNext = uint64(bitand(MdhParam.dmaLength,MDH_HEX2DEC_FFFFFF)); % Only first 3 bytes for dmaLength
         if length(mdhScan) > MDH_SCANSIZE
            sizeFirstChanInBytesNext = uint64(bitshift(extract_single_mdh_parameter(mdhScan(MDH_SCANSIZE+1:MDH_SCANSIZE+4),'chan','sizeChanInBytes'),-8));
         else
            sizeFirstChanInBytesNext = uint64(0);
         end
         bHasACQEND = bHasACQEND || MdhBlock(iBlock).isAcqEnd;
         if (dmaLengthNext == dmaLength) && (bKeepScanNext == bKeepScan)
            MdhBlock(iBlock).nScan = MdhBlock(iBlock).nScan + 1;
            iRelativePosition = iRelativePosition + dmaLength;
         else
            iBlock = iBlock + 1;
            offsetInFile =MdhBlock(iBlock-1).offsetInFile + MdhBlock(iBlock-1).nScan * dmaLength;
            flagDMALengthChange = dmaLength ~= dmaLengthNext;
            dmaLength = dmaLengthNext;
            bKeepScan = bKeepScanNext;
            sizeFirstChanInBytes = sizeFirstChanInBytesNext;
            MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, sizeFirstChanInBytes, MdhBitFlags, MdhParam);
            if flagDMALengthChange % No need to reread a chunk if dma length has not changed
               break;
            end
            iStart = iStart - int64(dmaLength); % backup and re-read this data into a the new block
         end
         if bHasACQEND
            break;
         end
         iStart = iStart + int64(dmaLength);
      end % for over block of scans with same dmalength
   end % while over full file
   if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
      textMessage = 'Pass 1 of 3: Parsing blocks';
      waitbar(double(iRelativePosition)/double(sizeAllScansInBytes), Control.GUI.hWaitBar, textMessage);
   end
   if ((iRelativePosition < sizeAllScansInBytes) && bHasACQEND)         
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'WARNING','ACQ_END before end of file.');
      end
   end
   if ~bHasACQEND
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'WARNING','No ACQ_END encountered.  Data possibly incomplete or corrupted.');
      end
   end
   nBlockSYNCDATA = sum([MdhBlock.isSyncData]);
   nBlockACQEND = sum([MdhBlock.isAcqEnd]);
   indexScan = find(([MdhBlock.isAcqEnd] == 0) & ([MdhBlock.isSyncData] == 0) & ([MdhBlock.keepScan] > 0));
   nBlockScanAll = length(indexScan);
   if nBlockScanAll < 1
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'Warning',sprintf('No readouts detected for meas %d.', iMeas),'Warning');
      end
      %fclose(fid);
      sizeAllScansInBytes = [];
      sizeChanInBytes = [];
      MdhBlock = [];
      nBlockSYNCDATA = [];
      nBlockACQEND = [];
      nBlockScan = [];
      blockNScan = [];
      nScanInMeas = [];
      dmaLengthScan = [];
      nSamplesInScan = [];
      return
   end

   nSamplesInScanUnique = unique([MdhBlock(indexScan).nSamplesInScan]);
   nScanInMeas = uint64(sum([MdhBlock(indexScan).nScan]));
   if (length(nSamplesInScanUnique) == 1)
      nBlockScan = nBlockScanAll;
      blockNScan = [MdhBlock(indexScan).nScan];
      if (length(unique([MdhBlock(indexScan).dmaLength])) ~= 1)
         if ~MrParcRaidFileCell{iMeas}.Control.User.silent
            if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
               display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'Warning',sprintf('Multiple different dmaLength in meas %d',iMeas),'Warning');
            else
               display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'ERROR',sprintf('Multiple different dmaLength in meas %d',iMeas),'Error');
               keyboard
               fclose(fid);
            end
         end
         sizeAllScansInBytes = [];
         sizeChanInBytes = [];
         MdhBlock = [];
         nBlockSYNCDATA = [];
         nBlockACQEND = [];
         nBlockScan = [];
         blockNScan = [];
         nScanInMeas = [];
         dmaLengthScan = [];
         nSamplesInScan = [];
         return
      end
      dmaLengthScan = MdhBlock(indexScan(1)).dmaLength;
      nSamplesInScan = nSamplesInScanUnique(1);
      sizeChanInBytes = MdhBlock(indexScan(1)).sizeFirstChanInBytes;
   else
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         if MrParcRaidFileCell{iMeas}.Control.User.scanEvalInfoMaskOnly
            display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'Warning',sprintf('Multiple different samples per scan in meas %d',iMeas),'Warning');
         else
            display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'ERROR',sprintf('Multiple different samples per scan in meas %d',iMeas),'Error');
            keyboard
         end
      end
   end
end

% =============================================================================
function PMUOut = extract_pmu(fid,MdhBlock,Control)
% -----------------------------------------------------------------------------
% Read the PMU data from SYNCDATA lines
% -----------------------------------------------------------------------------

   global MDH_SCANSIZE MDH_PMUSIZE

   PMU.TimeStamp.nColumn = 0;
   PMU.TimeStamp.raw = reshape([],2,0);
   PMU.TimeStamp.data = reshape([],2,0);
   Waveform.typeName = '';
   Waveform.nColumn = 0;
   Waveform.nRow = 0;
   Waveform.period = 0;
   Waveform.data = [];
   Waveform.trigger = [];
   pmuTypes = {'ECG1', 'ECG2', 'ECG3', 'ECG4', 'PULS', 'RESP', 'EXT1', 'EXT2'};
   PMU.Waveforms = repmat(Waveform,1,length(pmuTypes));
   for iType=1:length(pmuTypes)
      PMU.Waveforms(iType).typeName = pmuTypes(iType);
   end
   for iBlock = 1:length(MdhBlock)

      % Skip non-SYNCDATA
      if ~MdhBlock(iBlock).isSyncData 
         continue
      end

      fseek(fid,MdhBlock(iBlock).offsetInFile,'bof');
      nScanToRead = MdhBlock(iBlock).nScan;
      blockOfScans = fread(fid,[MdhBlock(iBlock).dmaLength nScanToRead],'uchar=>uchar');
  
      offset = MDH_SCANSIZE + MDH_PMUSIZE;

      % PMU Header
      MdhParam = extract_mdh_parameters(blockOfScans(offset+1:end,:), 'pmu', {'timeStamp', 'duration', 'data'});
      pmuTimeStamp = reshape(MdhParam.timeStamp,2,nScanToRead);
      pmuDuration = MdhParam.duration;

      PMU.TimeStamp.raw(:,PMU.TimeStamp.nColumn+1:PMU.TimeStamp.nColumn + nScanToRead) = pmuTimeStamp;
      PMU.TimeStamp.data(:,PMU.TimeStamp.nColumn+1:PMU.TimeStamp.nColumn + nScanToRead) = double(pmuTimeStamp) * 2.5d-3;
      PMU.TimeStamp.nColumn = PMU.TimeStamp.nColumn + nScanToRead;

      %pmuDataOld = blockOfScans(offset+17:end,:);
      nPmuPoints = idivide(MdhBlock(iBlock).dmaLength - offset - uint64(16),uint64(4),'fix');
      pmuData = reshape(MdhParam.data,nPmuPoints,nScanToRead);

      iPmuPoint = 1;
      while iPmuPoint <= nPmuPoints

         pmuType = bitshift(bitand(1.0*pmuData(iPmuPoint,:),hex2dec('010F0000'),'uint32'),-16) - 256;
         if length(unique(pmuType)) ~= 1
            if ~Control.User.silent
               display_message(Control.User.noGui,'ERROR','Unexpected PMU type data','PMU Problem');
            end
            fclose(fid);
            PMUOut = [];
            return;
         else
            pmuType = unique(pmuType);
         end
         % Check for end of data magic
         if (pmuType < 1) || (pmuType > 8)
            break 
         end
         pmuPeriod = pmuData(iPmuPoint+1,:);
         nPeriod = idivide(pmuDuration(:),pmuPeriod(:),'fix');
         if (length(unique(nPeriod)) ~= 1) || (length(unique(pmuPeriod(:))) ~= 1)
            if ~Control.User.silent
               display_message(Control.User.noGui,'ERROR','Unexpected PMU period or point data','PMU Problem');
            end
            fclose(fid);
            PMUOut = [];
            return;
         else
            nPeriod = unique(nPeriod);
            pmuPeriod = unique(pmuPeriod);
         end
         
         pmuPoint = uint32(iPmuPoint);
         pmuThisData = pmuData(pmuPoint+2:pmuPoint+2+nPeriod-1,:);
         % Waveforms are the lower 12 bits
         pmuWaveform = bitand(1.0*pmuThisData,hex2dec('00000fff'),'uint32');
         % Triggers in higher bits
         pmuTriggers = bitand(1.0*pmuThisData,hex2dec('fffff000'),'uint32');

         if PMU.Waveforms(pmuType).period == 0
            PMU.Waveforms(pmuType).period = (double(pmuPeriod) * 100d-6);
         else
            if PMU.Waveforms(pmuType).period ~= (double(pmuPeriod)*100d-6)
               if ~Control.User.silent
                  display_message(Control.User.noGui,'ERROR','Unexpected PMU period','PMU Problem');
               end
               fclose(fid);
               PMUOut = [];
               return;

            end
         end
         PMU.Waveforms(pmuType).data(1:nPeriod,PMU.Waveforms(pmuType).nColumn+1:PMU.Waveforms(pmuType).nColumn + nScanToRead) = pmuWaveform;
         PMU.Waveforms(pmuType).trigger(1:nPeriod,PMU.Waveforms(pmuType).nColumn+1:PMU.Waveforms(pmuType).nColumn + nScanToRead) = pmuTriggers;
         PMU.Waveforms(pmuType).nRow(PMU.Waveforms(pmuType).nColumn+1:PMU.Waveforms(pmuType).nColumn + nScanToRead) = nPeriod;
         PMU.Waveforms(pmuType).nColumn = PMU.Waveforms(pmuType).nColumn + nScanToRead;

         % Need to jump our pointer forward
         iPmuPoint = iPmuPoint + 2 + double(nPeriod);

      end

   end % for over blocks

   % Consolidate values
   WaveformOut.typeName = '';
   WaveformOut.nPoint = 0;
   WaveformOut.time = [];
   WaveformOut.data = [];
   WaveformOut.trigger = [];
   pmuTypes = {'ECG1', 'ECG2', 'ECG3', 'ECG4', 'PULS', 'RESP', 'EXT1', 'EXT2'};
   PMUOut.TimeStamp = PMU.TimeStamp;
   PMUOut.Waveforms = repmat(WaveformOut,1,length(pmuTypes));
   for iType=1:length(pmuTypes)
      PMUOut.Waveforms(iType).typeName = pmuTypes(iType);
      nRowTotal = sum(PMU.Waveforms(iType).nRow(:));
      PMUOut.Waveforms(iType).data = uint32(0) * zeros(1,nRowTotal,'uint32');
      PMUOut.Waveforms(iType).trigger = uint32(0)*zeros(1,nRowTotal,'uint32');
      PMUOut.Waveforms(iType).time = 0 * zeros(1,nRowTotal);
      iStart = 1;
      PMUOut.Waveforms(iType).nPoint = nRowTotal;
      for iCol=1:PMU.Waveforms(iType).nColumn

         nRow = PMU.Waveforms(iType).nRow(iCol);
         timeRelative = (0:(nRow-1)) * PMU.Waveforms(iType).period;
         timeStart = PMU.TimeStamp.data(2,iCol);
         timeValues = timeStart + timeRelative;
         waveValues = PMU.Waveforms(iType).data(1:nRow,iCol);
         triggers = PMU.Waveforms(iType).trigger(1:nRow,iCol);

         PMUOut.Waveforms(iType).time(iStart:iStart+nRow-1) = timeValues(:);
         PMUOut.Waveforms(iType).data(iStart:iStart+nRow-1) = waveValues(:);
         PMUOut.Waveforms(iType).trigger(iStart:iStart+nRow-1) = triggers(:);

         iStart = iStart + nRow;

      end
   end

   clear PMU;

end

% =============================================================================
function [dimNames, nScan, n, kSpaceCentreLineNo, kSpaceCentrePartitionNo, ...
          order, dimOut, dimOutProd, channelIDUniq] = ...
         scan_mdh_for_dim_and_order(fid,Dim,Control,sizeChanInBytes,MdhBlock,nScanInMeas,nSamplesInScan)
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------

   global MDH_SCANSIZE MAX_BYTES_PER_FREAD
   global index_keys

   iScanRead = uint64(0);
   timeRemain = uint64(10000000);
   dimNames  = {'CHA', 'LIN', 'SLC', 'PAR', 'ACQ', 'ECO', 'PHS', 'REP', 'SET', 'SEG', 'IDA', 'IDB', 'IDC', 'IDD', 'IDE'};
   indexPrev.dumby = 0; % Placeholder
   dimValue = uint16(0)*uint16(zeros(1,length(dimNames)));
   dimRank = uint64(zeros(1,length(dimNames)));
   nScan = 0;
   if ~(Control.User.noGui || Control.User.silent)
      time0 = tic;
   end

   kSpaceCentreLineNo = [];
   kSpaceCentrePartitionNo = [];

   % Instead of looping through blocks, we partition them into meta-blocks
   %    of the same number of samples per scan
   iBlock = 1;
   iBlockMeta = 1;
   MdhBlockMeta = append_to_mdhblockmeta(iBlock, iBlockMeta, MdhBlock);
   for iBlock=2:length(MdhBlock)

      indexMatch = find((cell2mat({MdhBlockMeta.isSyncData}) == MdhBlock(iBlock).isSyncData) & ...
                        (cell2mat({MdhBlockMeta.isAcqEnd}) == MdhBlock(iBlock).isAcqEnd) & ...
                        (cell2mat({MdhBlockMeta.keepScan}) == MdhBlock(iBlock).keepScan) & ...
                        (cell2mat({MdhBlockMeta.dmaLength}) == MdhBlock(iBlock).dmaLength) & ...
                        (cell2mat({MdhBlockMeta.nSamplesInScan}) == MdhBlock(iBlock).nSamplesInScan) & ...
                        (cell2mat({MdhBlockMeta.nChannelUsed}) == MdhBlock(iBlock).nChannelUsed));
      if numel(indexMatch) > 0
         MdhBlockMeta(indexMatch(1)).nBlock = MdhBlockMeta(indexMatch(1)).nBlock + 1;
         MdhBlockMeta(indexMatch(1)).iBlock(MdhBlockMeta(indexMatch(1)).nBlock) = iBlock;
         MdhBlockMeta(indexMatch(1)).nScan = MdhBlockMeta(indexMatch(1)).nScan + MdhBlock(iBlock).nScan;
      else 
         iBlockMeta = iBlockMeta + 1;
         MdhBlockMeta = append_to_mdhblockmeta(iBlock, iBlockMeta, MdhBlock, MdhBlockMeta);
      end
   end % for over MdhBlock

   for iBlockMeta = 1:length(MdhBlockMeta)

      % Skip SYNCDATA and blocked lines
      if MdhBlockMeta(iBlockMeta).isSyncData || MdhBlockMeta(iBlockMeta).isAcqEnd || ~MdhBlockMeta(iBlockMeta).keepScan 
         continue
      end

      for iBlockInner = 1:length(MdhBlockMeta(iBlockMeta).iBlock)

         iBlock = MdhBlockMeta(iBlockMeta).iBlock(iBlockInner);

         maxScanPerRead = max(idivide(MAX_BYTES_PER_FREAD,MdhBlock(iBlock).dmaLength,'fix'),uint64(1));
         nScanToRead = min(maxScanPerRead,MdhBlock(iBlock).nScan);

         fseek(fid,MdhBlock(iBlock).offsetInFile,'bof');
         nScanRemainInBlock = MdhBlock(iBlock).nScan;

         iScanReadBlock = uint64(0);
         while nScanRemainInBlock > 0

            if ~(Control.User.noGui || Control.User.silent)
               dTimeElapsed = toc(time0);
               rate = dTimeElapsed / (double(iScanRead) + 1);
               if dTimeElapsed > 20
                  timeRemain = min(double(nScanInMeas - iScanRead) * rate,timeRemain);
               else
                  timeRemain = double(nScanInMeas - iScanRead) * rate;
               end
               textMessage = sprintf('Pass 2 of 3\nScanning dimensions from MDH\nTime remaining (s): %d', uint32(timeRemain));
               waitbar(double(iScanRead)/double(nScanInMeas), Control.GUI.hWaitBar, textMessage);
            end

            blockOfScans = fread(fid,[MdhBlock(iBlock).dmaLength nScanToRead],'uchar=>uchar');
            MdhParam = extract_mdh_parameters(blockOfScans, 'scan', index_keys);
            if Control.User.collapseSeg
               MdhParam.indexSEG = MdhParam.indexSEG * uint16(0);
            end
            MdhBitFlags = extract_mdh_bitflags(blockOfScans, {'isAcqEnd', 'isSyncData', 'bKeepEIM', 'bDropEIM', 'isNoiseAdj', 'isPhaseCor', 'isRawDataCorrection', 'isReflect'});

            bHasACQEND = max(MdhBitFlags.isAcqEnd);
            if bHasACQEND
               break
            end
            channelID = uint16(0) * zeros(MdhBlock(iBlock).nChannelUsed,nScanToRead,'uint16');
            for iCha = 1:MdhBlock(iBlock).nChannelUsed
               offset  = uint64(MDH_SCANSIZE) + uint64(iCha-1) * sizeChanInBytes;
               channelID(iCha,:) = reshape(extract_single_mdh_parameter(blockOfScans(offset+1:offset+sizeChanInBytes,:),'chan','channelID'),[1,nScanToRead]);
            end
            channelIDUniq = unique(channelID(:));
            if prod(channelIDUniq == channelID(:,1)) == 0
               disp('Channels are not in monotonic order in data')
               keyboard
            end
            channelIDUniqNum = length(channelIDUniq);
            if (channelIDUniqNum ~= MdhBlock(iBlock).nChannelUsed) || (sum(channelIDUniq == channelID(:,1)) ~= MdhBlock(iBlock).nChannelUsed)
               if ~Control.User.silent
                  display_message(Control.User.noGui,'ERROR',['Mismatch between nChannel and coils found: ' int2str(channelIDUniqNum) ' ' int2str(MdhBlock(iBlock).nChannelUsed)],'Dimension Conflict');
               end
               fclose(fid);
               dimNames = [];
               nScan = [];
               n = [];
               kSpaceCentreLineNo = [];
               kSpaceCentrePartitionNo = [];
               order = [];
               dimOut = [];
               dimOutProd = [];
               channelIDUniq = [];
               return;
            end
            maskKeep = MdhBitFlags.bKeepEIM & ~MdhBitFlags.bDropEIM;
            index = find(maskKeep);
            if nScan == 0
               kSpaceCentreLineNo = extract_single_mdh_parameter(blockOfScans(:,index(1)),'scan','kSpaceCentreLineNo');
               kSpaceCentrePartitionNo = extract_single_mdh_parameter(blockOfScans(:,index(1)),'scan','kSpaceCentrePartitionNo');
            end
            nScan = nScan + numel(index);
            dimValue = max_from_mdh_param(dimValue, MdhBlock(iBlock).nChannelUsed, dimNames, MdhParam, index(:));
            if numel(index) > 1
               dimRank = rank_from_mdh_param(dimRank, dimNames, MdhParam, MdhBlock(iBlock).nChannelUsed, index(:));
            end
            if numel(index) == 1
               dimRank(1) = dimRank(1) + uint64(MdhBlock(iBlock).nChannelUsed);
            end
            if (iBlockInner > 1) && ~isempty(indexPrev) 
               dimRank = rank_from_index_prev(dimRank, ...
                                              dimNames, ...
                                              MdhParam, ...
                                              MdhBlock(iBlock).nChannelUsed, ...
                                              index(1), ...
                                              indexPrev);
            end
            if MdhBlock(iBlock).nChannelUsed < 2
               dimRank(1)=uint64(0);
            end
            if numel(index) > 0
               indexPrev = set_index_prev(dimNames, MdhParam, index(end));
            else
               if ~isempty(indexPrev)
                  indexPrev = [];
               end
            end

            iScanRead = iScanRead + nScanToRead;
            iScanReadBlock = iScanReadBlock + nScanToRead;
            nScanRemainInBlock = MdhBlock(iBlock).nScan - iScanReadBlock;
            nScanToRead = min(nScanToRead,nScanRemainInBlock);

         end % while over subsets of blocked lines

      end % for over blocks of non-syncdata

   end % over blocks of the same size

   if ~(Control.User.noGui || Control.User.silent)
      waitbar(0.5, Control.GUI.hWaitBar, 'MDH Scan complete');
   end
   clear blockOfScans;
   clear channelID MdhParam MdhBitFlags;
   clear index;
   clear iScanRead iScanReadInBlock nScanToRead;

   dimValue = dimValue .* uint16(dimRank > 0) + 1;
   for dimCell = {{'CHA' 1  false ''    ''     false ''      ''                   ''      ''    } ...
                  {'LIN' 2  false ''    ''     true  'Recon' 'nFourierLines'      'Recon' 'nLin'} ...
                  {'SLC' 3  false ''    ''     true  ''      ''                   'Recon' 'nSlc'} ...
                  {'PAR' 4  false ''    ''     true  'Recon' 'nFourierPartitions' 'Recon' 'nPar'} ...
                  {'ACQ' 5   true 'Raw' 'nAcq' false ''      ''                   ''      ''    } ...
                  {'ECO' 6   true 'Raw' 'nEco' false ''      ''                   ''      ''    } ...
                  {'PHS' 7   true 'Raw' 'nPhs' false ''      ''                   ''      ''    } ...
                  {'REP' 8   true 'Raw' 'nRep' false ''      ''                   ''      ''    } ...
                  {'SET' 9   true 'Raw' 'nSet' false ''      ''                   ''      ''    } ...
                  {'SEG' 10  false ''    ''    false ''      ''                   ''      ''    } ...
                  {'IDA' 11  false ''    ''    false ''      ''                   ''      ''    } ...
                  {'IDB' 12  false ''    ''    false ''      ''                   ''      ''    } ...
                  {'IDC' 13  false ''    ''    false ''      ''                   ''      ''    } ...
                  {'IDD' 14  false ''    ''    false ''      ''                   ''      ''    } ...
                  {'IDE' 15  false ''    ''    false ''      ''                   ''      ''    }}
      % Initialize with dimValue
      n.(dimCell{1}{1}) = dimValue(dimCell{1}{2});
      % Update with max from Dim structure
      if dimCell{1}{3} && (n.(dimCell{1}{1}) > 1)
         n.(dimCell{1}{1}) = max(n.(dimCell{1}{1}), Dim.(dimCell{1}{4}).(dimCell{1}{5}));
      end
      % Special handling in FillToFullFourier cases
      if dimCell{1}{6}
         if ~Control.User.noFillToFullFourier
            if ~isempty(dimCell{1}{7}) 
               n.(dimCell{1}{1}) = max(n.(dimCell{1}{1}), Dim.(dimCell{1}{7}).(dimCell{1}{8}));
            end
            if ~isempty(dimCell{1}{9}) > 0
               n.(dimCell{1}{1}) = max(n.(dimCell{1}{1}), Dim.(dimCell{1}{9}).(dimCell{1}{10}));
            end
         end
      end
   end
   dimOut = uint64([n.CHA n.LIN n.SLC n.PAR n.ACQ n.ECO n.PHS n.REP n.SET n.SEG n.IDA n.IDB n.IDC n.IDD n.IDE]);
   [~, order] = sort(dimRank, 2, 'descend');
   %keyboard
   dimOut = dimOut(order);
   dimOutProd = [uint64(1) compute_cumulative_product(dimOut(1:end-1))];

   if Control.User.noFillToFullFourier == true 
      n.COL = nSamplesInScan;
   else
      n.COL = max(nSamplesInScan, uint64(Dim.Recon.readoutFTLength));
   end

   clear dimValue dimRank;

end

% =============================================================================
function [status, rot_matrix] = quaternion_to_rotation_matrix(quaternion)
% -----------------------------------------------------------------------------
% PURPOSE:   Based on MdhProxy class method, converts quaternions into 
%               rotation matrix relating PRS to SagCorTra
% NOTES:     From MdhProxy.h
%        cout << "mdh::getRotMatrix: Gp = " << adRotMatrix[0][0] << " * Gsag + " <<
%                                              adRotMatrix[0][1] << " * Gcor + " <<
%                                              adRotMatrix[0][2] << " * Gtra\n";
%        cout << "mdh::getRotMatrix: Gr = " << adRotMatrix[1][0] << " * Gsag + " <<
%                                              adRotMatrix[1][1] << " * Gcor + " <<
%                                              adRotMatrix[1][2] << " * Gtra\n";
%        cout << "mdh::getRotMatrix: Gs = " << adRotMatrix[2][0] << " * Gsag + " <<
%                                              adRotMatrix[2][1] << " * Gcor + " <<
%                                              adRotMatrix[2][2] << " * Gtra\n";
%        cout << "mdh::getRotMatrix: ...finished" << endl; 
% -----------------------------------------------------------------------------

   MDH_QW = 1;
   MDH_QX = 2;
   MDH_QY = 3;
   MDH_QZ = 4;

   if sum(quaternion(:)) == 0
      status = 1;
      rot_matrix = [];
      return
   end

   ds = 2.0 / (quaternion(MDH_QW) * quaternion(MDH_QW) + ...
               quaternion(MDH_QX) * quaternion(MDH_QX) + ...
               quaternion(MDH_QY) * quaternion(MDH_QY) + ...
               quaternion(MDH_QZ) * quaternion(MDH_QZ));

   dxs = quaternion(MDH_QX) *  ds;
   dys = quaternion(MDH_QY) *  ds;
   dzs = quaternion(MDH_QZ) *  ds;

   dwx = quaternion(MDH_QW) * dxs;
   dwy = quaternion(MDH_QW) * dys;
   dwz = quaternion(MDH_QW) * dzs;

   dxx = quaternion(MDH_QX) * dxs;
   dxy = quaternion(MDH_QX) * dys;
   dxz = quaternion(MDH_QX) * dzs;

   dyy = quaternion(MDH_QY) * dys;
   dyz = quaternion(MDH_QY) * dzs;
   dzz = quaternion(MDH_QZ) * dzs;

   rot_matrix = zeros(3);
   rot_matrix(1,1) = 1.0 - (dyy + dzz);
   rot_matrix(1,2) =        dxy + dwz;
   rot_matrix(1,3) =        dxz - dwy;
   rot_matrix(2,1) =        dxy - dwz;
   rot_matrix(2,2) = 1.0 - (dxx + dzz);
   rot_matrix(2,3) =        dyz + dwx;
   rot_matrix(3,1) =        dxz + dwy;
   rot_matrix(3,2) =        dyz - dwx;
   rot_matrix(3,3) = 1.0 - (dxx + dyy);

   status = 0;

   return

end


% =============================================================================
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------

