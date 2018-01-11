
h_1=[1/32,1/32];
left=0;
right=1;
bottom=0;
top=1;
e=0.000001;
b=32;
a_1=sqrt(1/2);
a_2=sqrt(1/2);

[uh,A]=poisson_solver_triangle(left,right,bottom,top,h_1,e,b,a_1,a_2);
  uh=32*uh;  
 uh=poisson_solver_triangle1(uh,A,left,right,bottom,top,h_1,e,b,a_1,a_2);
   

Z=reshape(uh,33,33);
x=linspace(0,1,33);
y=linspace(0,1,33);
[X,Y]=meshgrid(x,y);
mesh(X,Y,Z)

