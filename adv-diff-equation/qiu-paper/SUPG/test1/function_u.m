function r = function_u(x,y,a_1,a_2,e)
%ËãÀý1
% r=(x.^2/(2*a_1)+e*x/(a_1^2)+(1/(2*a_1)+e/(a_1^2))*(exp(-a_1/e)-exp(-a_1/e*(1-x)))/(1-exp(-a_1/e))).*(y.^2/(2*a_2)+e*y/(a_2^2)+(1/(2*a_2)+e/(a_2^2))*(exp(-a_2/e)-exp(-a_2/e*(1-y)))/(1-exp(-a_2/e)));
%Àý×Ó2
r=sin(pi*x).*sin(pi*y);

