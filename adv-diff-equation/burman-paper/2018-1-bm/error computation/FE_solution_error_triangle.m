function r=FE_solution_error_triangle(mesh_type,uh,accurate_function,left,right,bottom,top,h_1,derivative_degree_x,derivative_degree_y,Gauss_point_number)


N1=(right-left)/h_1(1);
N2=(top-bottom)/h_1(2);
if mesh_type==1
    [M,T]=mesh_divise_1(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
elseif mesh_type==2
    [M,T]=mesh_divise_2(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
elseif mesh_type==3
    [M,T]=mesh_divise_3(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;

elseif mesh_type==4
    [M,T]=mesh_divise_4(left,right,bottom,top,h_1);
    number_of_elements=4*N1*N2;
end

 T_basis=T;


[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

r=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    vertices=M(:,T(1:3,n));
    uh_local=uh(T_basis(1:3,n));
    r=r+Gauss_quadrature_for_volume_integral_FE_solution_error_triangle(uh_local,accurate_function,vertices,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,derivative_degree_x,derivative_degree_y);
end
r=sqrt(r);