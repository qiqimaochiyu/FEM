close all
x=linspace(0,1,33);
y=x;
[x,y]=meshgrid(x,y);
%Z=exp(-(x-0.5).^2/0.2-3*(y-0.5).^2/0.2);
Z=1/2*(1-tanh((x-0.5)/0.05));
colormap('jet');
mesh(x,y,Z)
colorbar

