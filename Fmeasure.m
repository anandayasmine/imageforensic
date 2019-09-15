function [ output ] = Fmeasure(maskS,maskP)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

maskP = maskP(:);
maskS = maskS(:);

maskP1 = maskP;
maskP0 = (maskP==0); %reverse maskP

measure = struct;
measure.TP = sum(maskS & maskP1);
measure.FP = sum(maskS & maskP0);
measure.FN = sum(not(maskS) & maskP1);
measure.TN = sum(not(maskS) & maskP0);

bil_pembagi     = (measure.TP+measure.FP+measure.TP+measure.FN);
bil_pembagi_ACC = (measure.TP+measure.FP+measure.TN+measure.FN);
if(bil_pembagi == 0)
    bil_pembagi = bil_pembagi+0.1;
    measure.FM  = (2* measure.TP./bil_pembagi)*100;
    measure.ACC = ((measure.TP+measure.TN)./bil_pembagi_ACC)*100;
else
    measure.FM  = (2* measure.TP./bil_pembagi)*100;
    measure.ACC = ((measure.TP+measure.TN)./bil_pembagi_ACC)*100;
end
output = measure;

end

