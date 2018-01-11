close all
global mesh_type
for i=4
left=0;
right=1;
bottom=0;
top=1;
mesh_type=3;

    a_1=-cos(pi*55/180);a_2=-sin(pi*55/180);
    %a_1=0.15;a_2=0.1;
    %a_1=1/2;a_2=sqrt(3)/2;
    
e=10^(-5);b=0;
h_1=[1/2,1/2]*1/2^i;
uh=SUPG_solver_triangle(left,right,bottom,top,h_1,e,a_1,a_2);
x=linspace(0,1,1/h_1(1)+1);
y=x;
[X,Y]=meshgrid(x,y);
if mesh_type<3
    Z=reshape(uh,1/h_1(1)+1,1/h_1(1)+1);
else
    uh=SUPG_solver_triangle_3(left,right,bottom,top,h_1,e,a_1,a_2);
    Z=reshape(uh(1:((1/h_1(1)+1)*(1/h_1(1)+1))),1/h_1(1)+1,1/h_1(1)+1);
end
%close all
%MESH;
figure
colormap(jet);
mesh(X,Y,Z);
view([45,15])
figure
colormap(jet);
contour(X,Y,Z,15);
colorbar;
figure;colormap(jet);B=Z(:,1/(2*h_1(1))+1)';plot(y,B);
figure;colormap(jet);C=Z(1/(2*h_1(1))+1,:);plot(x,C);
end