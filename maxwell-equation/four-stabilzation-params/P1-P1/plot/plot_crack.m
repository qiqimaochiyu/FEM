function plot_crack 
h=16;
    x1=linspace(-1,1,2*h+1);
    y1=linspace(0,1,h+1);
    [x,y]=meshgrid(x1,y1);
    [z,Z,~,~]=function_u(x,y,1);
    close all
    colormap('jet')
    subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');hold on
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');hold on
    subplot(2,2,4),contour(x,y,Z,20);xlabel('contour u2');hold on
    subplot(2,2,3),mesh(x,y,Z);xlabel('exact u2');hold on
    
    x2=linspace(-1,1,2*h+1);
    y2=linspace(-1,0,h+1);
    [x,y]=meshgrid(x2,y2);
    [~,~,z,Z]=function_u(x,y,2);
    subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');
    subplot(2,2,4),contour(x,y,Z,20);xlabel('contour u2');
    subplot(2,2,3),mesh(x,y,Z);xlabel('exact u2');
    
function [r1,r2,r3,r4]=function_u(x,y,index) 
     r=sqrt(x.^2+y.^2);
     t=atan2(y,x);
     if index==1
        t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
        r1=-1/2*r.^(-1/2).*sin(t/2).*(x+1).*(y+1);
        r2=1/2*r.^(-1/2).*cos(t/2).*(x+1).*(y+1);
        r3=0;r4=0;
     elseif index==2
        t=t+2*pi;
        r1=0;r2=0;
        r3=-1/2*r.^(-1/2).*sin(t/2).*(x+1).*(y+1);
        r4=1/2*r.^(-1/2).*cos(t/2).*(x+1).*(y+1);
     end
     
end
end
