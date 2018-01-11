function r=Gauss_quadrature_volume_integral_FE_solution_error_triangle_2(uh_local1,accurate_function1,vertices,derivative_degree1_x,derivative_degree1_y,uh_local2,accurate_function2,derivative_degree2_x,derivative_degree2_y)

global Gauss_coefficient_reference_triangle 

Gpn=length(Gauss_coefficient_reference_triangle);

r=0;
[Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(vertices);

for i=1:Gpn
    r=r+Gauss_coefficient_local_triangle(i)*(feval(accurate_function1,Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2))-FE_solution_triangle(Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),uh_local1,vertices,derivative_degree1_x,derivative_degree1_y))*(feval(accurate_function2,Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2))-FE_solution_triangle(Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),uh_local2,vertices,derivative_degree2_x,derivative_degree2_y));
end