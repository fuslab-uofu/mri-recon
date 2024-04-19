function [Time, StrTime] = AcqTime(hdr, Str)
%added by Dylan E. Palomino
%the Siemens time stamps are converted (to hour:minute:second format) and exported as strings to put in "shd.ulAcquisitionStrTimeStamp"
SecTot = double(hdr.ulTimeStamp)*0.0025; % In Seconds
MinTot = SecTot/60;
HrTot = MinTot/60;
Hr = HrTot-mod(HrTot,1);
MinRem = mod(HrTot,1)*60;
Min = MinRem-mod(MinRem,1);
SecRem = mod(MinRem,1)*60;
Sec = SecRem-mod(SecRem,1);
Sec1 = round(mod(SecRem,1),2);
Time(:,1) = Hr;
Time(:,2) = Min;
Time(:,3) = SecRem;
if Str == 1
    StrHr = num2str(Hr);
    for j = 1:length(StrHr)
        if StrHr(j,1) == ' '
            StrHr(j,1) = '0';
        end
        if size(StrHr,2) == 1 % DEP 02/01/19
            StrHr(:,2) = StrHr(:,1); 
            StrHr(:,1) = '0';
        end 
    end
    StrMin = num2str(Min);
    for j = 1:length(StrMin)
        if StrMin(j,1) == ' '
            StrMin(j,1) = '0';
        end
    end
     StrSec = num2str(SecRem,'%0.4f');
    for j = 1:length(StrSec)
        if StrSec(j,1) == ' '
            StrSec(j,1) = '0';
        end
    end
    for j = 1:length(Time)
        StrTime(j,:) = [StrHr(j,:), ':', StrMin(j,:), ':', StrSec(j,:)];
    end
end
% TimeStr = num2str(Time);
% for j = 1:length(SecTot)
%     StrHr(j) = num2str(Hr(j));
%     StrMin = num2str(Min);
%     StrSec = num2str(Sec);
%     if length(StrHr) == 1
%         StrHr = ['0',StrHr]
%     end
%     if length(StrMin) == 1
%         StrMin = ['0',StrMin];
%     end
%     if length(StrSec) == 1
%         StrSec = ['0',StrSec];
%     end
%     Time(j,:) = [StrHr,':',StrMin,':',StrSec];
% end
end