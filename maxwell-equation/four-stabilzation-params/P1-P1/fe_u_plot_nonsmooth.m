function fe_u_plot_nonsmooth(domain_type,uh1,uh2,h)

hf=1/(2*h(1));
if domain_type==1 %方形区域
    x=linspace(0,1,2*hf+1);
    y=x; 
    [x,y]=meshgrid(x,y);
    z=reshape(uh1,2*hf+1,2*hf+1);Z=reshape(uh2,2*hf+1,2*hf+1);
    close all
    colormap('jet')
    subplot(2,2,1),mesh(x,y,z);xlabel('exact u1');
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');
    subplot(2,2,4),contour(x,y,Z,20);xlabel('exact u2');
    subplot(2,2,3),mesh(x,y,Z);xlabel('exact u2');

elseif domain_type==2 %L型区域
    x1=linspace(-1,0,hf+1);
    y1=linspace(-1,1,2*hf+1);
    [x,y]=meshgrid(x1,y1);
    tnp1=(2*hf+1)*(hf+1);
    z=reshape(uh1(1:tnp1),2*hf+1,hf+1);
    Z=reshape(uh2(1:tnp1),2*hf+1,hf+1);
    close all
    colormap('jet')
    subplot(2,2,1),mesh(x,y,z);xlabel('fe u1');hold on
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');hold on
    subplot(2,2,4),contour(x,y,Z,20);xlabel('contour u2');hold on
    subplot(2,2,3),mesh(x,y,Z);xlabel('fe u2');hold on
    
    x2=linspace(0,1,hf+1);
    y2=x2;
    [x,y]=meshgrid(x2,y2);
    tnp2=2*hf*(hf+1)+1;tnp3=(3*hf+1)*(hf+1);
    z=reshape(uh1(tnp2:tnp3),hf+1,hf+1);
    Z=reshape(uh2(tnp2:tnp3),hf+1,hf+1);
    subplot(2,2,1),mesh(x,y,z);xlabel('fe u1');
    subplot(2,2,2),contour(x,y,z,20);xlabel('contour u1');
    subplot(2,2,4),contour(x,y,Z,20);xlabel('contour u2');
    subplot(2,2,3),mesh(x,y,Z);xlabel('fe u2');
    
end