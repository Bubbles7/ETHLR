function [AC,MIhat,recall, precision, Fmeasure] = compute_metrics(V,gnd)
nClass = length(unique(gnd));
label = kmeans(V,nClass);
res = bestMap(gnd,label);
%=============  evaluate AC: accuracy ==============
AC = length(find(gnd == res))/length(gnd);
MIhat = MutualInfo(gnd,label);%ÅÐ¶Ï¾ÛÀàÐ§¹û
% disp(['Clustering in the HNMF AC and MI: ',num2str(AC), ',' , num2str(MIhat)]);
[recall, precision, Fmeasure]=Value_compute(gnd,res);
end