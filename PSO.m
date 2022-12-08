clc;
clear all;
addpath([pwd, '/Dataset']);
addpath([pwd, '/funs']);
tic

%% Load Dataset
load('PAAD_CHOL_ESCA.mat');                              % raw sample
%% PA-CH-ES
X = PAAD_CHOL_ESCA;
labels = PAAD_CHOL_ESCA_Samplecategory;

[n1,n2,n3] = size(X);                   % sample size
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
                 

%% Initialize the population
N = 10;                         % Number of initial populations
d = 6;                          
ger =  10;                         
limit = [0.1,   100; 
    0.1    100;
    0.1     100;
    0.1,         1; 
    -7, 7;
    -7,     7];
vlimit = [-0.01,   0.01; 
    -0.01,   0.01; 
    -0.01,   0.01; 
    -0.01,   0.01; 
    -1,     1;
    -1,     1];                 % speed limit
w = 0.8;                        % inertia weight
c1 = 0.5;                       
c2 = 0.5;                        
for i = 1:d
    x(:,i) = limit(i,1) + (limit(i,2) - limit(i,1)) .* rand(N, 1); % The location of the initial population
end
v = rand(N, d);                  % The velocity of the initial population
xm = x;                          % The historical best position for each individual
ym = zeros(1, d);                % The historical best position of the population
fxm = zeros(N, 1);               % The historical best fitness of each individual
fym = -inf;                      % The historical best fitness of the population

%% Population update
iter = 1;
record = zeros(ger, 1);          % Recorder
fx = zeros(N,1);
while iter <= ger
    parfor i = 1:N
        fx(i) = tnn(x(i,:),Xn,labels, G) ; % Individual current fitness   
    end
     for i = 1:N      
        if fxm(i) < fx(i)
            fxm(i) = fx(i);     
            xm(i,:) = x(i,:);   
        end 
     end
    if fym < max(fxm)
            [fym, nmax] = max(fxm);   
            ym = xm(nmax, :);      
     end
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x); % Velocity update
    % Boundary velocity
    for ii = 1:d
        v(v(:,ii) > vlimit(ii,2),ii) = vlimit(ii,2);
        v(v(:,ii) < vlimit(ii,1),ii) = vlimit(ii,1);
    end
    x = x + v; % Location update
    % Boundary position
    for ii = 1:d
        x(x(:,ii) > limit(ii,2),ii) = limit(ii,2);
        x(x(:,ii) < limit(ii,1),ii) = limit(ii,1);
    end
    record(iter) = fym; % Maximum value record
    iter = iter+1;
end
figure(3);plot(record);title('Convergence process')

disp(['maximum value£º',num2str(fym)]);
disp(['variable values£º',num2str(ym)]);
t=toc;