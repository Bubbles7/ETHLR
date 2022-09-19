%% Load Setting
clear
addpath([pwd, '/Dataset']);
addpath([pwd, '/funs']);


%% Load Dataset
load('PAAD_CHOL_ESCA.mat');                              % raw sample
%% PA-CH-ES
X = PAAD_CHOL_ESCA;
labels = PAAD_CHOL_ESCA_Samplecategory;


[n1,n2,n3] = size(X);  
n=min(n1,n2);

%% Hyper-Parameters
% Weighted vector of weighted tensor Schatten p-norm
w = [];
w = [w; 9.6*ones(1,1)];
w = [w; 20*ones(1,1)];
w = [w; 69*ones(n-2,1)];
% The power p of weighted tensor Schatten p-norm£¬p in (0,1]
p = 0.97;
beta = 10^(-2.36);
lambda = 10^(-2.47);


Xn = X;
rhos = 0.7;
ind = find(rand(n1*n2*n3,1)<rhos);
Xn(ind) = rand(length(ind),1);
nClass = length(unique(labels));


%% hyper-graph
for i=1:n2
    Xnn = reshape(Xn(:,i,:),n1,n3);
    G(:,i,:)=ConstructHCFW(Xnn',n3,6);
end
                
%% Optimizataion
[D, E, obj, err, iter]  = ethlr_tnn_lp(Xn, lambda, w, p, beta, G);
DD = permute(D,[3 1 2]);
XE = [DD(:,:,1),DD(:,:,2),DD(:,:,3)]; 
rand('twister',5489);
[AC,MIhat,recall, precision, Fmeasure] = compute_metrics(XE,labels);

