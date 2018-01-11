function result_output(L2,Le)

m=length(L2);n=length(Le);
k1=zeros(1,m-1);k2=zeros(1,n-1);
for i=1:m-1
    k1(i)=log2(L2(i)/L2(i+1));
end
for j=1:n-1
    k2(j)=log2(Le(j)/Le(j+1));
end
fprintf '------ L2 error and order ------\n'
disp(L2)
disp(k1)

fprintf '------ L_curl error and order ------\n'
disp(Le)
disp(k2)




