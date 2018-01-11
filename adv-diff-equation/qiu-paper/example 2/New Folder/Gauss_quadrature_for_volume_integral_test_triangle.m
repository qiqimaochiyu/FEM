function r=Gauss_quadrature_for_volume_integral_test_triangle(coefficient_function_name,Gauss_coefficient_local,Gauss_point_local,vertices,test_basis_index,test_derivative_degree_x,test_derivative_degree_y,a_1,a_2,e,b)


Gpn=length(Gauss_coefficient_local);
r=0;
for i=1:Gpn
      r=r+Gauss_coefficient_local(i)*feval(coefficient_function_name,Gauss_point_local(i,1),Gauss_point_local(i,2),a_1,a_2,e,b)*tri_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
end