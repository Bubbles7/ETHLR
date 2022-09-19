function [recall, precision, Fmeasure]=Value_compute(gnd,label)
[m,n] = size(gnd);
nclass = size(unique(gnd),1);
count=zeros(nclass,nclass);
gnd(m+1)=0;
for i=1:m
    count(gnd(i),label(i))= count(gnd(i),label(i))+1;
end
l=sum(count,1);
h=sum(count,2);
for i=1:nclass
    recall(i)=count(i,i)/l(i);
    precision(i)=count(i,i)/h(i);
end
    recall=mean(recall);
    precision=mean(precision);
    Fmeasure = 2/(1/precision + 1/recall);
end



        
