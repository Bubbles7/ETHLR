function [Lap]=ConstructHCFW(X,n,q)
% X,输入样本（n*m）维,n样本个数，返回值L为超图拉普拉斯矩阵（n*n）
%% 1.计算欧氏距离矩阵
%  X=XQ;
%  [mm,nn]=size(X);
%  n=mm;

Dis = dist2(X,X);%计算两两样本之间的欧氏距离

%% 2.寻找K个近邻
H = zeros(n,n);
E = zeros(n,q+1);
%% 2.1 选出第i类
for j = 1:n
 Dis_j = Dis(j,:)'; % 取第j个特征
 Dis_j_s = sort(Dis_j); % 对Dis_j进行排序
 [~,a,b] = intersect(Dis_j,Dis_j_s);% 计算Dis_j和Dis_j_s的交集,得到a的下标
 E(j,:) = a(1:q+1,:);
 H(a(1:q+1),j) = 1;
end

%% 3 计算权值矩阵
W = zeros(n,1);
for i = 1:n
    x = X(i,:);
    e = E(i,:);
    w = 0;
    w_m = 0;
    for j = 1:(q+1)
       x_e = X(e(j),:); 
       w_m =w_m+dist2(x,x_e);
    end
    w_m = w_m/(q+1);
    for j = 1:(q+1)
       x_e = X(e(j),:); 
        w = w+exp(-dist2(x,x_e)/(w_m^2));
%         w = w+exp(-dist2(x,x_e)/ss);
    end
    W(i,1) = w;
end
%% 4 计算顶点度矩阵
De=diag(sum(H));
W = diag(W);
Dv=diag(sum((H*W)'));
% NewW=H*W*De*H';
% NewW=max(NewW,NewW');
% Dv_I=1./diag(Dv);
% Dv_I=diag(Dv_I);
% Lap=sqrt(Dv_I)*NewW*sqrt(Dv_I);
Lap = Dv-H*W*(De^(-1))*H';
end