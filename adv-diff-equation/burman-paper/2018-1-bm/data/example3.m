function data=example3
global b
a_1=1;a_2=0;e=10^(-5); 

% QIU 
%例子3中的数据
data.g=@g;data.a=@a;data.f=@f;data.ux=@u_x;data.uy=@u_y;data.u=@u;
    
    function r=u(x,y)
    r=(x.^2/(2*a_1)+e*x/(a_1^2)+(1/(2*a_1)+e/(a_1^2))*(exp(-a_1/e)-exp(-a_1/e*(1-x)))/(1-exp(-a_1/e))).*(y.^2/(2*a_2)+e*y/(a_2^2)+(1/(2*a_2)+e/(a_2^2))*(exp(-a_2/e)-exp(-a_2/e*(1-y)))/(1-exp(-a_2/e)));
    end

    function r=a(~,~)
        r=1;
    end

    function r=g(~,~)
        r=0;
    end
function r =u_x(x,y)
r=(e/a_1^2 + x/a_1 + (a_1*exp((a_1*(x - 1))/e)*(e/a_1^2 + 1/(2*a_1)))/(e*(exp(-a_1/e) - 1))).*(y.^2/(2*a_2) + ((e/a_2^2 + 1/(2*a_2))*(exp((a_2*(y - 1))/e) - exp(-a_2/e)))/(exp(-a_2/e) - 1) + (e*y)/a_2^2);
end
function r = u_y(x,y)
r=(e/a_2^2 + y/a_2 + (a_2*exp((a_2*(y - 1))/e)*(e/a_2^2 + 1/(2*a_2)))/(e*(exp(-a_2/e) - 1))).*(x.^2/(2*a_1) + ((e/a_1^2 + 1/(2*a_1))*(exp((a_1*(x - 1))/e) - exp(-a_1/e)))/(exp(-a_1/e) - 1) + (e*x)/a_1^2);
end
function r=f(x,y)

r=a_2*(e/a_2^2 + y/a_2 + (a_2*exp((a_2*(y - 1))/e)*(e/a_2^2 + 1/(2*a_2)))/(e*(exp(-a_2/e) - 1))).*(x.^2/(2*a_1) + ((e/a_1^2 + 1/(2*a_1))*(exp((a_1*(x - 1))/e) - exp(-a_1/e)))/(exp(-a_1/e) - 1) + (e*x)/a_1^2) - e*((1/a_2 + (a_2^2*exp((a_2*(y - 1))/e)*(e/a_2^2 + 1/(2*a_2)))/(e^2*(exp(-a_2/e) - 1))).*(x.^2/(2*a_1) + ((e/a_1^2 + 1/(2*a_1))*(exp((a_1*(x - 1))/e) - exp(-a_1/e)))/(exp(-a_1/e) - 1) + (e*x)/a_1^2) + (1/a_1 + (a_1^2*exp((a_1*(x - 1))/e)*(e/a_1^2 + 1/(2*a_1)))/(e^2*(exp(-a_1/e) - 1))).*(y.^2/(2*a_2) + ((e/a_2^2 + 1/(2*a_2))*(exp((a_2*(y - 1))/e) - exp(-a_2/e)))/(exp(-a_2/e) - 1) + (e*y)/a_2^2)) + a_1*(e/a_1^2 + x/a_1 + (a_1*exp((a_1*(x - 1))/e)*(e/a_1^2 + 1/(2*a_1)))/(e*(exp(-a_1/e) - 1)))*(y^2/(2*a_2) + ((e/a_2^2 + 1/(2*a_2))*(exp((a_2*(y - 1))/e) - exp(-a_2/e)))/(exp(-a_2/e) - 1) + (e*y)/a_2^2) + b*(x.^2/(2*a_1) + ((e/a_1^2 + 1/(2*a_1))*(exp((a_1*(x - 1))/e) - exp(-a_1/e)))/(exp(-a_1/e) - 1) + (e*x)/a_1^2).*(y.^2/(2*a_2) + ((e/a_2^2 + 1/(2*a_2))*(exp((a_2*(y - 1))/e) - exp(-a_2/e)))/(exp(-a_2/e) - 1) + (e*y)/a_2^2);

end
end