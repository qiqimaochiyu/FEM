format short e

left=0;
right=1;
bottom=0;
top=1;
a_1=1/2;a_2=sqrt(3/4);
b=0;e=10^(-6);

h_1=[1/16,1/16]*(1/2);
uh=poisson_solver_triangle(left,right,bottom,top,h_1,e,b,a_1,a_2);
x=linspace(0,1,1/h_1(1)+1);
y=x;
[X,Y]=meshgrid(x,y);
Z=reshape(uh,1/h_1(1)+1,1/h_1(1)+1);
mesh(X,Y,Z);
figure
contour(X,Y,Z)