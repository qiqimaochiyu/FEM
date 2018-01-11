format short e
for i=1:4
left=0;
right=1;
bottom=0;
top=1;
a_1=1/2;a_2=sqrt(3/4);
e=10^(-8);
Gauss_point_number=9;


h_1=[1/4,1/4]*(1/2)^i;
uh=poisson_solver_triangle(left,right,bottom,top,h_1,e,a_1,a_2);
L2_error(i)=FE_solution_error_triangle(uh,'function_u',left,right,bottom,top,h_1,0,0,Gauss_point_number);
L2_u(i)=FE_u_real_triangle('function_u',left,right,bottom,top,h_1,Gauss_point_number);
L2_abs(i)=L2_error(i)/L2_u(i);
H1_error_x=FE_solution_error_triangle(uh,'function_u_x',left,right,bottom,top,h_1,1,0,Gauss_point_number);
H1_error_y=FE_solution_error_triangle(uh,'function_u_y',left,right,bottom,top,h_1,0,1,Gauss_point_number);
H1_error(i)=sqrt(H1_error_x^2+H1_error_y^2);
H1_error_u_x=FE_u_real_triangle('function_u_x',left,right,bottom,top,h_1,Gauss_point_number);
H1_error_u_y=FE_u_real_triangle('function_u_y',left,right,bottom,top,h_1,Gauss_point_number);
H1_u(i)=sqrt(H1_error_u_x^2+H1_error_u_y^2);
H1_abs(i)=H1_error(i)/H1_u(i);

end
L2_abs
H1_abs
k=[log2(L2_error(1)/L2_error(2)) log2(L2_error(2)/L2_error(3)) log2(L2_error(3)/L2_error(4))]
r=[log2(H1_error(1)/H1_error(2)) log2(H1_error(2)/H1_error(3)) log2(H1_error(3)/H1_error(4))]

