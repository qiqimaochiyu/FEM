
left=0;
right=1;
bottom=0;
top=1;
e=0.0001;
b=0;
a_1=1;
a_2=0;

h_1=[1/2,1/2];
uh=poisson_solver_triangle(left,right,bottom,top,h_1,e,b,a_1,a_2);


% A=reshape(uh,33,33);
% x=linspace(0,1,33);
% y=linspace(0,1,33);
% [X,Y]=meshgrid(x,y);
% mesh(X,Y,A);
% figure;
% contour(X,Y,A);
% figure;
% B=A(:,17)';
% plot(y,B);
% figure
% C=A(17,:);
% plot(x,C);