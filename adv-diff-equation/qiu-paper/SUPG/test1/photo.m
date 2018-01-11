format short e
for i=1:3
left=0;
right=1;
bottom=0;
top=1;
a_1=1/2;a_2=sqrt(3/4);
e=10^(-5);

h_1=[1/4,1/4]*(1/2)^i;
uh=poisson_solver_triangle(left,right,bottom,top,h_1,e,a_1,a_2);
x=linspace(0,1,1/h_1(1)+1);
y=x;
[X,Y]=meshgrid(x,y);
Z=reshape(uh,1/h_1(1)+1,1/h_1(1)+1);
subplot(1,3,i);
mesh(X,Y,Z);
end