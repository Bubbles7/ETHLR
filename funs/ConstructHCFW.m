function [Lap]=ConstructHCFW(X,n,q)
% X,����������n*m��ά,n��������������ֵLΪ��ͼ������˹����n*n��
%% 1.����ŷ�Ͼ������
%  X=XQ;
%  [mm,nn]=size(X);
%  n=mm;

Dis = dist2(X,X);%������������֮���ŷ�Ͼ���

%% 2.Ѱ��K������
H = zeros(n,n);
E = zeros(n,q+1);
%% 2.1 ѡ����i��
for j = 1:n
 Dis_j = Dis(j,:)'; % ȡ��j������
 Dis_j_s = sort(Dis_j); % ��Dis_j��������
 [~,a,b] = intersect(Dis_j,Dis_j_s);% ����Dis_j��Dis_j_s�Ľ���,�õ�a���±�
 E(j,:) = a(1:q+1,:);
 H(a(1:q+1),j) = 1;
end

%% 3 ����Ȩֵ����
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
%% 4 ���㶥��Ⱦ���
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