% --------------------------------------------------------------------------------
% FileInfo=jGetFileInfo(FileName)
% Jason Mendes - UCAIR - University of Utah - 2015
% --------------------------------------------------------------------------------
% This function opens the .dat file and retrieves information about the sequence
% data stored.
% --------------------------------------------------------------------------------
function FileInfo=jGetFileInfo(FileName,PrintHeader)

if nargin < 2
    PrintHeader = 1;
end

% Prompt for filename if one is not specified
if ((exist('FileName')~=1)||isempty(FileName))
	[file_name,path_name]=uigetfile('*.dat','Select Raw Data File');
	FileName=sprintf('%s%s',path_name,file_name);
end

% Open the .dat file
FileInfo=[];
FileInfo.FileId=fopen(FileName,'r');
if (FileInfo.FileId<0)
   error('jGetFileInfo(): Could not open file specified');
end
fseek(FileInfo.FileId,0,'eof');
FileSize=ftell(FileInfo.FileId);
FileInfo.FileName=FileName;

% Find BaselineString
fseek(FileInfo.FileId,0,'bof');
while (ftell(FileInfo.FileId)<FileSize)
   NextLine=fgets(FileInfo.FileId);
   if (~isempty(findstr(NextLine,'<ParamString."tMeasuredBaselineString">')))
      Data=jReadSimpleString(FileInfo,NextLine);
      if (~isempty(findstr(lower(Data.String),'n4_v')))
         FileInfo.Version=Data.String;
         break;
      end
   elseif (~isempty(findstr(NextLine,'<ParamString."tBaselineString">')))
      Data=jReadSimpleString(FileInfo,NextLine);
      if (~isempty(findstr(lower(Data.String),'n4_v')))
         FileInfo.Version=Data.String;
         break;
      end
   end
end

% If no version is found try and guess one
if (~isfield(FileInfo,'Version'))
   fseek(FileInfo.FileId,0,'bof');
   Data=fread(FileInfo.FileId,2,'uint32');
   if ((Data(1)<10000)&&(Data(2)<=64))
      FileInfo.Version='BestMatchIsVD';
   else
      FileInfo.Version='BestMatchIsVB';
   end
end
if (~isempty(findstr(lower(FileInfo.Version),'vb')))
   FileInfo.IsVB=1;
   FileInfo.IsVD=0;
elseif (~isempty(findstr(lower(FileInfo.Version),'vd')))
   FileInfo.IsVB=0;
   FileInfo.IsVD=1;
elseif (~isempty(findstr(lower(FileInfo.Version),'ve')))
   FileInfo.IsVB=0;
   FileInfo.IsVD=1;
else
   error('jGetFileInfo(): Unknown version');
end

% Find how many measurements are in the file
if (FileInfo.IsVB)
   fseek(FileInfo.FileId,0,'bof');
   FileInfo.Header.StartPosition=0;
   FileInfo.Header.EndPosition=FileSize;
   FileInfo.Header.DataPosition=fread(FileInfo.FileId,1,'uint32');
elseif (FileInfo.IsVD)
   fseek(FileInfo.FileId,4,'bof');
   NumberOfMeasurements=fread(FileInfo.FileId,1,'uint32');
   for MeasIndex=1:NumberOfMeasurements
      FileInfo.Header(MeasIndex).MeasId=fread(FileInfo.FileId,1,'uint32');
      FileInfo.Header(MeasIndex).FileId=fread(FileInfo.FileId,1,'uint32');
      FileInfo.Header(MeasIndex).StartPosition=fread(FileInfo.FileId,1,'uint64');
      FileInfo.Header(MeasIndex).EndPosition=FileInfo.Header(MeasIndex).StartPosition+fread(FileInfo.FileId,1,'uint64');
      FileInfo.Header(MeasIndex).PatientName=sprintf('%s',fread(FileInfo.FileId,64,'uchar=>char'));
      FileInfo.Header(MeasIndex).ProtocolName=sprintf('%s',fread(FileInfo.FileId,64,'uchar=>char'));
   end
   for MeasIndex=1:NumberOfMeasurements
      if (fseek(FileInfo.FileId,FileInfo.Header(MeasIndex).StartPosition,'bof'))
         error('jGetFileInfo(): Data set not complete');
      end
      FileInfo.Header(MeasIndex).DataPosition=FileInfo.Header(MeasIndex).StartPosition+fread(FileInfo.FileId,1,'uint32');
   end
else
   error('jGetFileInfo(): Unknown file type');
end

% Load all protocols
FileInfo.Protocol=jGetProtocol(FileInfo,PrintHeader);


% --------------------------------------------------------------------------------
% FileInfo=jLoadBasicProtocol(FileName)
% Jason Mendes - UCAIR - University of Utah - 2015
% --------------------------------------------------------------------------------
% This function opens the .dat file and retrieves information about the sequence
% protocol from the ASCCONV section and some additional parameters.
% --------------------------------------------------------------------------------
function Protocol=jGetProtocol(FileInfo,PrintHeader)

% Prepare to read protocol
Protocol=[];

% There are some useful parameters not in the ASCCONV section
ParameterTag(1).String='<ParamLong."lNoOfPhaseCorrScans">';
ParameterTag(2).String='<ParamLong."relSliceNumber">';
ParameterTag(3).String='<ParamLong."NSetMeas">';
ParameterTag(4).String='<ParamLong."RampupTime">';
ParameterTag(5).String='<ParamLong."FlattopTime">';
ParameterTag(6).String='<ParamLong."DelaySamplesTime">';
ParameterTag(7).String='<ParamLong."RegridMode">';
ParameterTag(8).String='<ParamDouble."ADCDuration">';
ParameterTag(9).String='<ParamLong."alRegridRampupTime">';
ParameterTag(10).String='<ParamLong."alRegridFlattopTime">';
ParameterTag(11).String='<ParamLong."alRegridDelaySamplesTime">';
ParameterTag(12).String='<ParamLong."alRegridMode">';
ParameterTag(13).String='<ParamDouble."aflRegridADCDuration">';
ParameterTag(14).String='<ParamLong."NoOfFourierPartitions">';
ParameterTag(15).String='<ParamString."PatientPosition">';

ParameterTag(16).String='<ParamString."PatientBirthDay">';

% Find protocol for each requested sequence
for SearchIndex=1:length(FileInfo.Header)
   % Read in the header
   HeaderStart=FileInfo.Header(SearchIndex).StartPosition;
   HeaderEnd=FileInfo.Header(SearchIndex).DataPosition;
   if (fseek(FileInfo.FileId,HeaderStart,'bof'))
      error('jGetProtocol(): Data header not complete');
   end
   Header=fread(FileInfo.FileId,HeaderEnd-HeaderStart,'char=>char')';
   
%    %Svedin: Start
%    if PrintHeader
%        fhid = fopen([FileInfo.FileName(1:end-4) '_hdr.txt'], 'w');
%        fprintf(fhid, '%c', Header);
%        fclose(fhid);
%    end
%    %Svedin: End
   
   % Find ASCCONV section
   StartIndex=max(findstr(Header,'ASCCONV BEGIN'));
   EndIndex=max(findstr(Header,'ASCCONV END'));
   if ((~isempty(StartIndex))&&(~isempty(EndIndex)))
      String=Header(StartIndex:EndIndex);
      String(find(String==';'))=':'; % Remove original ';'
      String=regexprep(String,'\n',';\n');
      String=String((find(String==10,1)+1):find(String==10,1,'last'));
      String(find(String<=32))=[]; % Remove spaces
      String(find(String=='_'))=[]; % '_' cause problems
      String=regexprep(String,'[','('); % Matlab uses round brackets for indices
      String=regexprep(String,']',')'); % Matlab uses round brackets for indices
      String=regexprep(String,'"',''''); % Matlab uses round brackets for indices
      String=regexprep(String,'(','(1+'); % Matlab is not zero based
      String=regexprep(String,'0x(\w*)','hex2dec(''$1'')'); % Convert hex values
      LineEnd=[0 find(String==';')];
      for Index=2:length(LineEnd)
         EvalString=String((LineEnd(Index-1)+1):LineEnd(Index));
         if isempty(findstr(lower(EvalString),'attribute'))
            try
               eval(sprintf('Protocol(%d).%s',SearchIndex,EvalString));
            catch
               fprintf('Error adding %s\n',EvalString');
            end
         end
      end
   end   
   % Read in extra parameters
   for Index=1:length(ParameterTag);
      StartIndex=min(findstr(Header,ParameterTag(Index).String));
      EndIndex=find(Header(StartIndex:end)=='}',1);
      if ((~isempty(StartIndex))&&(~isempty(EndIndex)))
         CurrentString=Header((1:EndIndex)+StartIndex);
         NameIndex=find(CurrentString=='"',2);
         ValueIndex=(find(CurrentString=='{',1,'last')+1):(find(CurrentString=='}',1,'last')-1);
         if (length(NameIndex==2)&&length(ValueIndex>0))
            NameString=CurrentString((NameIndex(1)+1):(NameIndex(2)-1));
            ValueString=CurrentString(ValueIndex);
            if (~isempty(findstr(lower(CurrentString),'double')))
               ValueString=ValueString((min(findstr(lower(ValueString),'<precision>'))+12):end);
               if (~isempty(find(ValueString==10,1)))
                   ValueString=ValueString((find(ValueString==10,1)+1):end);
               end
               ValueString=sprintf('%s',ValueString(find(ValueString<=32,1):end));
            elseif (~isempty(findstr(lower(CurrentString),'long')))
               ValueString=sprintf('%s',ValueString);
            elseif (~isempty(findstr(lower(CurrentString),'string')))
               if (length(find(ValueString=='"'))<2)
                  ValueString='''''';
               else
                  ValueString(find(ValueString=='"'))='''';
               end
            else
               % Not supported at this time
               ValueString='[]';
            end
            EvalString=sprintf('Protocol(%d).%s=[%s];',SearchIndex,NameString,ValueString);
            EvalString(find(EvalString==10))=32;
            try
               eval(EvalString);
            catch
               fprintf('Error %s\n',EvalString');
               pause;
            end
         end
      end
   end
   % Compatibility issues
    if (isfield(Protocol(SearchIndex),'ADCDuration'))
        if (~isempty(Protocol(SearchIndex).ADCDuration))
            Protocol(SearchIndex).aflRegridADCDuration=Protocol(SearchIndex).ADCDuration;
        end
    end
    if (isfield(Protocol(SearchIndex),'RampupTime'))
        if (~isempty(Protocol(SearchIndex).RampupTime))
            Protocol(SearchIndex).alRegridRampupTime=Protocol(SearchIndex).RampupTime;
        end
    end
    if (isfield(Protocol(SearchIndex),'FlattopTime'))
        if (~isempty(Protocol(SearchIndex).FlattopTime))
            Protocol(SearchIndex).alRegridFlattopTime=Protocol(SearchIndex).FlattopTime;
        end
    end
    if (isfield(Protocol(SearchIndex),'DelaySamplesTime'))
        if (~isempty(Protocol(SearchIndex).DelaySamplesTime))
            Protocol(SearchIndex).alRegridDelaySamplesTime=Protocol(SearchIndex).DelaySamplesTime;
        end
    end
    if (isfield(Protocol(SearchIndex),'RegridMode'))
        if (~isempty(Protocol(SearchIndex).RegridMode))
            Protocol(SearchIndex).alRegridMode=Protocol(SearchIndex).RegridMode;
        end
    end
    if (isfield(Protocol(SearchIndex),'sWiPMemBlock'))
        if (~isempty(Protocol(SearchIndex).sWiPMemBlock))
            Protocol(SearchIndex).sWipMemBlock=Protocol(SearchIndex).sWiPMemBlock;
        end
    end
   % Fix problem with partition count for 2D
   if (isfield(Protocol(SearchIndex),'sKSpace'))
      if (Protocol(SearchIndex).sKSpace.ucDimension<3)
         Protocol(SearchIndex).sKSpace.lPartitions=1;
      end
   end
   % Single repetition is not recorded
   if (~isfield(Protocol(SearchIndex),'lRepetitions'))
      Protocol(SearchIndex).lRepetitions=0;
   end
   if (isempty(Protocol(SearchIndex).lRepetitions))
       Protocol(SearchIndex).lRepetitions=0;
   end      
    % Number of coils
    if (isfield(Protocol(SearchIndex),'sCoilSelectMeas'))
        Protocol(SearchIndex).asCoilSelectMeas=Protocol(SearchIndex).sCoilSelectMeas;
    end
    if (isfield(Protocol(SearchIndex),'asCoilSelectMeas'))
        if (FileInfo.IsVB)
            Protocol(SearchIndex).lNumberOfChannels=length(Protocol(SearchIndex).asCoilSelectMeas(1).asList);
        elseif (FileInfo.IsVD)
            Protocol(SearchIndex).lNumberOfChannels=length(Protocol(SearchIndex).sCoilSelectMeas(1).aRxCoilSelectData(1).asList);
        end
    end
   % Find location of PMU and image lines
   if (FileInfo.IsVB)
      fseek(FileInfo.FileId,FileInfo.Header(SearchIndex).DataPosition,'bof');
      Data=fread(FileInfo.FileId,3,'uchar=>uchar');
      LengthDMA=double(typecast([Data(1:3);0],'uint32'));
      Protocol(SearchIndex).OffsetPMU=[];
      Protocol(SearchIndex).OffsetImage=(FileInfo.Header(SearchIndex).DataPosition): LengthDMA:(FileInfo.Header(SearchIndex).EndPosition-LengthDMA);
   elseif (FileInfo.IsVD)
      for LoopIndex=1:2
         CountPMU=0;
         CountImage=0;
         CurrentPosition=FileInfo.Header(SearchIndex).DataPosition;
         while (CurrentPosition<=(FileInfo.Header(SearchIndex).EndPosition-192))
            fseek(FileInfo.FileId,CurrentPosition,'bof');
            Data=fread(FileInfo.FileId,3,'uchar=>uchar');
            LengthDMA=double(typecast([Data(1:3);0],'uint32'));
            fseek(FileInfo.FileId,CurrentPosition+40,'bof');
            Data=fread(FileInfo.FileId,4,'uchar=>uchar');
            Mask=double(typecast(Data,'uint32'));
            if (min(bitand(Mask,2^5),1))
               CountPMU=CountPMU+1;
               if (LoopIndex==2)
                  Protocol(SearchIndex).OffsetPMU(CountPMU)=CurrentPosition;
               end
            else
               if (LengthDMA>(Protocol(SearchIndex).lNumberOfChannels*32+192))
                  CountImage=CountImage+1;
                  if (LoopIndex==2)
                     Protocol(SearchIndex).OffsetImage(CountImage)=CurrentPosition;
                  end
               end
            end
            CurrentPosition=CurrentPosition+LengthDMA;
         end
         if (LoopIndex==1)
            Protocol(SearchIndex).OffsetPMU=zeros(1,CountPMU);
            Protocol(SearchIndex).OffsetImage=zeros(1,CountImage);
         end
      end
   end
   % Get readout length
   
end


% --------------------------------------------------------------------------------
% Data=jReadSimpleString(FileInfo,InitialString)
% Jason Mendes - UCAIR - University of Utah - 2014
% Update: Returns both data string name and value (October 2015)
% --------------------------------------------------------------------------------
% This function reads in a string parameter value from a .dat file.
% Returns Data.Name and Data.String containing the string name and value.
% --------------------------------------------------------------------------------
function Data=jReadSimpleString(FileInfo,InitialString)

% Read in data until first closing bracket is found
FinalString=InitialString;
Data=[];
while (isempty(find(FinalString=='}',1)))
   NextString=fgets(FileInfo.FileId,1000);
   FinalString=sprintf('%s %s',FinalString,NextString);
end

% Get field name
ValueIndex=find(FinalString=='"',2,'first');
if (length(ValueIndex)<2)
   DataName=sprintf('UnknownName%d',round(rand(1)*1000));
else
   DataName=FinalString(ValueIndex(1)+1:ValueIndex(2)-1);
end
DataName(find(DataName<33))=[]; % Ignore line feeds, etc...

% Get the values
ValueIndex=find(FinalString=='{',1,'last')+1:find(FinalString=='}',1,'last')-1;
StringText=FinalString(ValueIndex);
ValueIndex=find(StringText=='"',2,'first');
if (length(ValueIndex)<2)
   DataString='Error';
else
   DataString=StringText(ValueIndex(1)+1:ValueIndex(2)-1);
end

% Return values
Data.Name=DataName;
Data.String=DataString;