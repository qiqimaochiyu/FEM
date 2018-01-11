function exact_u_plot_smooth(uh1,uh2)

%square domain

global left right bottom top h

function [r1,r2]=function_u(x,y)

    r1=sin(pi*y).*cos(pi*x);
    r2=-sin(pi*x).*cos(pi*y);
end

N1=(right-left)/h(1);
N2=(top-bottom)/h(2);

z1=reshape(uh1,N1+1,N2+1);

z2=reshape(uh2,N1+1,N2+1);

close all
colormap('jet')
x=linspace(-1,1,N1+1);
y=x;
[x,y]=meshgrid(x,y);
[z,Z]=function_u(x,y);

subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');
subplot(2,2,2),mesh(x,y,z1);xlabel('FE u1');
subplot(2,2,3),contour(x,y,z);xlabel('u1 contour');
subplot(2,2,4),contour(x,y,z1);xlabel('FE u1 contour');

figure 
colormap('jet')
subplot(2,2,1),mesh(x,y,Z);xlabel('exact u2');
subplot(2,2,2),mesh(x,y,z2);xlabel('FE u2');
subplot(2,2,3),contour(x,y,Z);xlabel('u2 contour');
subplot(2,2,4),contour(x,y,z2);xlabel('FE u2 contour');
hold off
end

