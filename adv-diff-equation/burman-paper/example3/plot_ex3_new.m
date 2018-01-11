
global T b
close all
h_1=[1/32 1/32];h=h_1(1);t=[0.25 1.5];bt=[h h/4 h/8];
left=0;right=1;bottom=0;top=1;
e=10^(-6);h_num=3;
data=example_time;
for i=1:h_num
    T=bt(i);b=1/T;a_1=sqrt(1/2);a_2=sqrt(1/2);
    % P=[1.2 0.6 0.3];p=P(i);
    p=0.6;p_1=p;
    q=0.0125;Gauss_point_number=9;mesh_type=1;
    [A,uh,F1,F2,F3,F4]=B_P_solver_triangle(data,mesh_type,q,left,right,bottom,...
        top,h_1,e,b,p,a_1,a_2);
    uh=uh*b;
%load data3
nt=t(2)/T;
a=zeros(1089,nt);
a(:,1)=uh;
for j=2:nt
a(:,j)=B_P_solver_triangle_f(A,a(:,j-1),F1,F2,F3,F4,data,mesh_type,q,left,right,bottom,top,h_1,e,b,p,a_1,a_2)*b;
end
figure
for n=1
set(gcf, 'doublebuffer', 'on')
for k = 1:nt
Z=reshape(a(:,k),33,33)/b;
x=linspace(0,1,33);
y=linspace(0,1,33);
[X,Y]=meshgrid(x,y);
colormap(jet)
mesh(Y,X,Z)
view([135,45])
% view([-45 15]);
% contour(X,Y,Z,20)
    drawnow
end
end
m=max(max(Z));
disp(m)
end

    
 