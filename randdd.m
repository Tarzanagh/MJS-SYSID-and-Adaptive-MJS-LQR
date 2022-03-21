function x=randdd(A_e) 
x=zeros(1,20);
for i=1:20 
x(:,i)=max(abs(eig(A_e(:,:,i))));
end
end