
format short e
global mesh_type b h_1
mesh_type=3;
data=example3;
num=1;
L2_error=zeros(num,4);
L2_abs=zeros(num,4);
H1_error=zeros(num,4);
H1_abs=zeros(num,4);
H1_u=zeros(num,4);
L2_u=zeros(num,4);
r=zeros(num,3);
k=zeros(num,3);
for j=1:num
    b=1;
    for i=2
        h_1=[1/8,1/8]*(1/2)^i;
        fprintf 'start loop:' 
        disp(i)
        p=0;e=10^(-5);q=0.0125;
        left=0;right=1;bottom=0;top=1;
        a_1=sqrt(1/2);a_2=sqrt(3/4);
        Gauss_point_number=9;
uh=B_P_solver_triangle(q,left,right,bottom,top,h_1,e,b,p,a_1,a_2);
L2_error(j,i)=FE_solution_error_triangle(uh,data.u,left,right,bottom,top,h_1,0,0,Gauss_point_number);
L2_u(j,i)=FE_u_real_triangle(data.u,left,right,bottom,top,h_1,Gauss_point_number);
L2_abs(j,i)=L2_error(j,i)/L2_u(j,i);
H1_error_x=FE_solution_error_triangle(uh,data.ux,left,right,bottom,top,h_1,1,0,Gauss_point_number);
H1_error_y=FE_solution_error_triangle(uh,data.uy,left,right,bottom,top,h_1,0,1,Gauss_point_number);
H1_error(j,i)=sqrt(H1_error_x^2+H1_error_y^2+L2_error(j,i)^2);
H1_error_u_x=FE_u_real_triangle(data.ux,left,right,bottom,top,h_1,Gauss_point_number);
H1_error_u_y=FE_u_real_triangle(data.uy,left,right,bottom,top,h_1,Gauss_point_number);
H1_u(j,i)=sqrt(H1_error_u_x^2+H1_error_u_y^2+L2_u(j,i)^2);
H1_abs(j,i)=H1_error(j,i)/H1_u(j,i);
r(j,:)=log2([L2_abs(j,1)/L2_abs(j,2) L2_abs(j,2)/L2_abs(j,3) L2_abs(j,3)/L2_abs(j,4)]);
%L2_abs(3)/L2_abs(4)])
k(j,:)=log2([H1_abs(j,1)/H1_abs(j,2) H1_abs(j,2)/H1_abs(j,3) H1_abs(j,3)/H1_abs(j,4)]);
%H1_abs(3)/H1_abs(4)])
    end
end
L2_abs;H1_abs;
r;k;
x=linspace(0,1,1/h_1(1)+1);
y=x;[x,y]=meshgrid(x,y);
z1=sin(x*pi).*sin(pi*y);
H=(1/h_1(1)+1)^2;
z2=reshape(uh(1:H),1/h_1(1)+1,1/h_1(1)+1);
subplot(2,2,1);mesh(x,y,z1)
subplot(2,2,2);contour(x,y,z1)
subplot(2,2,3);mesh(x,y,z2)
subplot(2,2,4);contour(x,y,z2)

