function r = function_u(x,y,a_1,a_2,e)

% r=exp(-(x-0.5).^2/0.2-3*(y-0.5).^2/0.2);
 
r=1/2*(1-tanh((x-0.5)/0.05));


% f=@(x,y)exp(-(x-0.5).^2/0.2-3*(y-0.5).^2/0.2);
% ezmesh(f,[0,1],[0,1])
% 
% f=@(x,y)1/2*(1-tanh((x-0.5)/0.05));
% ezmesh(f,[0,1],[0,1])