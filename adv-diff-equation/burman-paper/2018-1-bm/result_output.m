function result_output(L2,H1,L_max)

m=length(L2);
k1=zeros(1,m-1);k2=k1;k3=k1;
for i=1:m-1
    k1(i)=log2(L2(i)/L2(i+1));
    k2(i)=log2(H1(i)/H1(i+1));
    k3(i)=log2(L_max(i)/L_max(i+1));
end

fprintf '------ L2 error and order ------\n'
disp(L2)
disp(k1)

fprintf '------ H1 error and order ------\n'
disp(H1)
disp(k2)

fprintf '------ L_max error and order ------\n'
disp(L_max)
disp(k3)





