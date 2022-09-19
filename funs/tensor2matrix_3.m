function [ XX ] = tensor2matrix_3( X)

[n1,n2,n3]=size(X);
XX=zeros(n1*n3,n2);
for i=1:n2
    im=X(:,i,:);
    im=reshape(im,[n1,n3]);
    XX(:,i)=im(:);
end

end

