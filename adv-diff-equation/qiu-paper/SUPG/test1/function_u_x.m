function r = function_u_x(x,y,a_1,a_2,e)
% test-1
% r=(e/a_1^2 + x/a_1 + (a_1*exp((a_1*(x - 1))/e)*(e/a_1^2 + 1/(2*a_1)))/(e*(exp(-a_1/e) - 1))).*(y.^2/(2*a_2) + ((e/a_2^2 + 1/(2*a_2))*(exp((a_2*(y - 1))/e) - exp(-a_2/e)))/(exp(-a_2/e) - 1) + (e*y)/a_2^2);
%test-2
r=pi*cos(pi*x).*sin(pi*y);