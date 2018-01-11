z=zeros(33,33);
for i=1:33
    for j=1:33
        if 1<=i && i<=17 && 1<=j && j<=9
            z(i,j)=1;
        else
            if 1<=i && i<=9 && 9<=j && j<=17
                z(i,j)=1;
            else
                z(i,j)=0;
            end
        end
    end
end
z
