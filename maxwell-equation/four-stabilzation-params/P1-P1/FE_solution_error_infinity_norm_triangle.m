function r=FE_solution_error_infinity_norm_triangle(uh,accurate_function,left,right,bottom,top,h_partition,derivative_degree_x,derivative_degree_y,Gauss_point_number)


N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);
number_of_elements=4*N1_partition*N2_partition;

[M_partition,T_partition]=mesh_divise1(left,right,bottom,top,h_partition);


    T_basis=T_partition;


[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

r=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    vertices=M_partition(:,T_partition(1:3,n));
    uh_local=uh(T_basis(1:3,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    temp=max(abs(feval(accurate_function,Gauss_point_local_triangle(:,1),Gauss_point_local_triangle(:,2))-FE_solution_triangle(Gauss_point_local_triangle(:,1),Gauss_point_local_triangle(:,2),uh_local,vertices,derivative_degree_x,derivative_degree_y)));    
    if temp>r
        r=temp;
    end
end