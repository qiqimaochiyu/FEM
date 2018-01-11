function r=FE_solution_error_triangle(uh,accurate_function,left,right,bottom,top,h_1,derivative_degree_x,derivative_degree_y,Gauss_point_number)



N1_partition=(right-left)/h_1(1);
N2_partition=(top-bottom)/h_1(2);
number_of_elements=2*N1_partition*N2_partition;

[M_partition,T_partition]=mesh_devise(left,right,bottom,top,h_1);


 T_basis=T_partition;


[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

r=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    vertices=M_partition(:,T_partition(:,n));
    uh_local=uh(T_basis(:,n));
    r=r+Gauss_quadrature_for_volume_integral_FE_solution_error_triangle(uh_local,accurate_function,vertices,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,derivative_degree_x,derivative_degree_y);
end
r=sqrt(r);