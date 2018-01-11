function r=FE_solution_error_triangle(uh,accurate_function,M,T,derivative_degree_x,derivative_degree_y)
% Yu wei 2017-05
% L2 ·¶ÊýÎó²î
number_of_elements=length(T);

r=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    vertices=M(:,T(:,n));
    uh_local=uh(T(:,n));
    r=r+Gauss_quadrature_for_volume_integral_FE_solution_error_triangle(uh_local,accurate_function,vertices,derivative_degree_x,derivative_degree_y);
end
