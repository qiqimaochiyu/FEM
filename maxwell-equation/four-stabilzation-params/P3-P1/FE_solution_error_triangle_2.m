function r=FE_solution_error_triangle_2(uh1,accurate_function1,M,T,derivative_degree1_x,derivative_degree1_y,uh2,accurate_function2,derivative_degree2_x,derivative_degree2_y)

number_of_elements=length(T);
r=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    vertices=M(:,T(:,n));
    uh_local1=uh1(T(:,n));
    uh_local2=uh2(T(:,n));
    r=r+Gauss_quadrature_volume_integral_FE_solution_error_triangle_2(uh_local1,accurate_function1,vertices,derivative_degree1_x,derivative_degree1_y,uh_local2,accurate_function2,derivative_degree2_x,derivative_degree2_y);
end