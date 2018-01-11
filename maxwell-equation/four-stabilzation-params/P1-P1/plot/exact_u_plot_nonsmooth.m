function exact_u_plot_nonsmooth(index,L_type)
% Yu wei 2017-05
%mesh and contour the exact u
%index:index of the example 
%L_type:the type of the domian; 1:square; 0:L_shape;

h=16;
if L_type==1 %方形区域
    x=linspace(0,1,2*h+1);
    y=x; 
    [x,y]=meshgrid(x,y);
    [z,Z]=function_u(x,y,index);
    close all
    colormap('jet')
    subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');
    subplot(2,2,4),contour(x,y,Z,20);xlabel('exact u2');
    subplot(2,2,3),mesh(x,y,Z);xlabel('exact u2');

elseif L_type==3 %裂缝区域
    x1=linspace(-1,1,2*h+1);
    y1=linspace(0,1,h+1);
    [x,y]=meshgrid(x1,y1);
    [z,Z]=function_u(x,y,index);
    close all
    colormap('jet')
    subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');hold on
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');hold on
    subplot(2,2,4),contour(x,y,Z,20);xlabel('contour u2');hold on
    subplot(2,2,3),mesh(x,y,Z);xlabel('exact u2');hold on
    
    x2=linspace(-1,1,2*h+1);
    y2=linspace(-1,0,h+1);
    [x,y]=meshgrid(x2,y2);
    [z,Z]=function_u(x,y,index);
    subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');
    subplot(2,2,4),contour(x,y,Z,20);xlabel('contour u2');
    subplot(2,2,3),mesh(x,y,Z);xlabel('exact u2');
    

elseif L_type==2 %L型区域
    x1=linspace(-1,1,2*h+1);
    y1=linspace(0,1,h+1);
    [x,y]=meshgrid(x1,y1);
    [z,Z]=function_u(x,y,index);
    close all
    colormap('jet')
    subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');hold on
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');hold on
    subplot(2,2,4),contour(x,y,Z,20);xlabel('contour u2');hold on
    subplot(2,2,3),mesh(x,y,Z);xlabel('exact u2');hold on
    
    x2=linspace(-1,0,h+1);
    y2=x2;
    [x,y]=meshgrid(x2,y2);
    [z,Z]=function_u(x,y,index);
    subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');
    subplot(2,2,4),contour(x,y,Z,20);xlabel('contour u2');
    subplot(2,2,3),mesh(x,y,Z);xlabel('exact u2');
    
end

%% function 
function [r1,r2]=function_u(x,y,index)
    if index==1
        r1=sin(pi*y).*cos(pi*x);
        r2=-sin(pi*x).*cos(pi*y);
    elseif index==2
        r=sqrt(x.^2+y.^2);
        t=atan2(y,x);
        t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
        r1=-2/3*r.^(-1/3).*sin(t/3);
        r2=2/3*r.^(-1/3).*cos(t/3);
     elseif index==2.2
        r=sqrt(x.^2+y.^2);
        t=atan2(y,x);
        t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
        r1=-2/3*r.^(-1/3).*sin(t/3)+sin(y);
        r2=2/3*r.^(-1/3).*cos(t/3)+sin(x);
    elseif index==2.3
        r=sqrt(x.^2+y.^2);
        t=atan2(y,x);
        t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
        r1=2/3*r.^(-1/3).*cos(t/3);
        r2=2/3*r.^(-1/3).*sin(t/3);       
    elseif index==3
        r=sqrt(x.^2+y.^2);
        t=atan2(y,x);
        t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
        r1=-1/2*r.^(-1/2).*sin(t/2);
        r2=1/2*r.^(-1/2).*cos(t/2);
    elseif index==3.2
        r=sqrt(x.^2+y.^2);
        t=atan2(y,x);
        t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
        r1=-1/2*r.^(-1/2).*sin(t/2)+sin(y);
        r2=1/2*r.^(-1/2).*cos(t/2)+sin(x);
    elseif index==3.3
        r=sqrt(x.^2+y.^2);
        t=atan2(y,x);
        t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
        r1=1/2*r.^(-1/2).*cos(t/2);
        r2=1/2*r.^(-1/2).*sin(t/2);
    elseif index==4
        r1=-sin(2*pi*y).*(cos(2*pi*x)-1);
        r2=sin(2*pi*x).*(cos(2*pi*y)-1);
    end
end
end
