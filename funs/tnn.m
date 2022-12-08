function AC = tnn(x,Xn,labels, G)
[n1,n2,n3] = size(Xn);                   % sample size
n=min(n1,n2);
w1 = x(1:3);
w = [];
w = [w; w1(1)*ones(1,1)];
w = [w; w1(2)*ones(1,1)];
w = [w; w1(3)*ones(n-2,1)];
lambda = 10^x(6);
p = x(4);
beta = 10^x(5);
[D, E, obj, err, iter]  = ethlr_tnn_lp(Xn, lambda, w, p, beta,G);
DD = permute(D,[3 1 2]);
XE = [DD(:,:,1),DD(:,:,2),DD(:,:,3)];   
[AC,MIhat,recall, precision, Fmeasure] = compute_metrics(XE,labels);
end
