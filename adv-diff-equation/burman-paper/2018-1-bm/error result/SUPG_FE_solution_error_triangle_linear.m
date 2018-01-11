
format short e
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
    %b=1;
    for i=1:4
        
 e=10^(-2);  
 left=0;
 right=1;
 bottom=0;
 top=1;
 a_1=sqrt(1/2);a_2=sqrt(3/4);

Gauss_point_number=9;

h_1=[1/8,1/8]*(1/2)^i;

uh=SUPG_solver_triangle(left,right,bottom,top,h_1,e,b,a_1,a_2);
L2_error(j,i)=FE_solution_error_triangle(uh,'function_u',left,right,bottom,top,h_1,0,0,Gauss_point_number,a_1,a_2,e);
L2_u(j,i)=FE_u_real_triangle('function_u',left,right,bottom,top,h_1,Gauss_point_number,a_1,a_2,e);
L2_abs(j,i)=L2_error(j,i)/L2_u(j,i);
H1_error_x=FE_solution_error_triangle(uh,'function_u_x',left,right,bottom,top,h_1,1,0,Gauss_point_number,a_1,a_2,e);
H1_error_y=FE_solution_error_triangle(uh,'function_u_y',left,right,bottom,top,h_1,0,1,Gauss_point_number,a_1,a_2,e);
H1_error(j,i)=sqrt(H1_error_x^2+H1_error_y^2+L2_error(j,i)^2);
H1_error_u_x=FE_u_real_triangle('function_u_x',left,right,bottom,top,h_1,Gauss_point_number,a_1,a_2,e);
H1_error_u_y=FE_u_real_triangle('function_u_y',left,right,bottom,top,h_1,Gauss_point_number,a_1,a_2,e);
H1_u(j,i)=sqrt(H1_error_u_x^2+H1_error_u_y^2+L2_u(j,i)^2);
H1_abs(j,i)=H1_error(j,i)/H1_u(j,i);
r(j,:)=log2([L2_abs(j,1)/L2_abs(j,2) L2_abs(j,2)/L2_abs(j,3) L2_abs(j,3)/L2_abs(j,4)]);
%L2_abs(3)/L2_abs(4)])
k(j,:)=log2([H1_abs(j,1)/H1_abs(j,2) H1_abs(j,2)/H1_abs(j,3) H1_abs(j,3)/H1_abs(j,4)]);
%H1_abs(3)/H1_abs(4)])
    end
end
L2_abs;
H1_abs;
r;
k;
