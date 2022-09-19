function [L,S,obj,err,iter] = ethlr_tnn_lp(X, lambda, weight, p, beta, G)

% Solve the Tensor Robust Principal Component Analysis based on Weighted Tensor Schatten p-Norm problem by ADMM
%
%
% ---------------------------------------------
% Input:
%       X       -    d1*d2*d3 tensor
%       G       -    d3*d2*d3 tensor
%       lambda  -    >0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
% Output:
%       L       -    d1*d2*d3 tensor
%       S       -    d1*d2*d3 tensor
%       obj     -    objective function value
%       err     -    residual
%       iter    -    number of iterations
%
%
%
% Written by Na Yu

tol = 1e-8;
max_iter = 100;
rho = 1.1;
mu1 = 1e-4;
mu2 = 1e-4;
max_mu = 1e10;

if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end

dim = size(X);
L = zeros(dim);
S = L;
Y1 = L;
Y2 = L;
W = L;
[n1,n2,n3]=size(X);
for iter = 1 : max_iter
    Lk = L;
    Sk = S;
    % Update L
    H1=X-S+Y1/mu1;
    H2=W+Y2/mu2;
    L1=mu1*H1+mu2*H2;
    L2=L1/(mu1+mu2);
    [L,tnnL] = prox_tnn(L2,weight/(mu1+mu2),p);
    % Update S
    S = prox_l1(-L+X+Y1/mu1,lambda/mu1);
    % update W
    for i=1:n2
        LL=reshape((L(:,i,:)),n1,n3);
        YY=reshape((Y2(:,i,:)),n1,n3);
        GG=reshape((G(:,i,:)),n3,n3);
        I = speye((n3));
        W(:,i,:)=(mu2*LL-YY)/(2*beta*GG+mu2*I);
    end
    dY = X-L-S;
    dY2 = W-L;
    % Coverge condition
    chgL = max(abs(Lk(:)-L(:)));
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([ chgL chgS max(abs(dY(:))) ]);
    
    if chg < tol
        break;
    end
    Y1 = Y1 + mu1*dY;
    Y2 = Y2 + mu2*dY2;
    mu1 = min(rho*mu1,max_mu);
    mu2 = min(rho*mu2,max_mu); 
    obj = tnnL+lambda*norm(S(:),1);
    err = norm(dY(:));
    disp(['The ' num2str(iter)  '-th Iteration' ', mu1=' num2str(mu1) ', mu2=' num2str(mu2) ', obj=' num2str(obj) ', err=' num2str(err)]);
end