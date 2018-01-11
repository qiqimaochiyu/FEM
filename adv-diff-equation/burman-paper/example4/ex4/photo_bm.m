
load data3
global b p v t a_1 a_2

b=0;data=example_h;
a_1=1;a_2=0;e=10^(-4);
b_p=0.5;q=0.0125;
Gauss_point_number=9;
mesh_type=1;   
uh=B_P_solver_triangle(data,q,e,b,b_p,a_1,a_2); 
result_mesh;  