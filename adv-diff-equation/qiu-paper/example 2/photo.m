
left=0;
right=1;
bottom=0;
top=1;
e=0.000001;
b=0;
a_1=-cos(55/180*pi);
a_2=-sin(55/180*pi);

h_1=[1/32,1/32];
uh=poisson_solver_triangle(left,right,bottom,top,h_1,e,b,a_1,a_2);
A=reshape(uh,1/h_1(1)+1,1/h_1(1)+1);
x=linspace(0,1,1/h_1(1)+1);
y=linspace(0,1,1/h_1(1)+1);
[X,Y]=meshgrid(x,y);
 mesh(X,Y,A);
colormap('jet')
%figure
contour(A,30);

% figure;
% B=A(:,17)';
% plot(y,B);
% 
% figure
% C=A(17,:);
% plot(x,C);
