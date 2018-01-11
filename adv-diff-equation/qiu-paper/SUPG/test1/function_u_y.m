function r = function_u_y(x,y,a_1,a_2,e)
%test-1
% r=(e/a_2^2 + y/a_2 + (a_2*exp((a_2*(y - 1))/e)*(e/a_2^2 + 1/(2*a_2)))/(e*(exp(-a_2/e) - 1))).*(x.^2/(2*a_1) + ((e/a_1^2 + 1/(2*a_1))*(exp((a_1*(x - 1))/e) - exp(-a_1/e)))/(exp(-a_1/e) - 1) + (e*x)/a_1^2);
%test-2
r=pi*cos(pi*y).*sin(pi*x);