function r=Gauss_quadrature_for_volume_integral_test_triangle1(uh_local,Gauss_coefficient_local,Gauss_point_local,vertices,test_basis_index,test_derivative_degree_x,test_derivative_degree_y)


Gpn=length(Gauss_coefficient_local);
r=0;
for i=1:Gpn
      r=r+Gauss_coefficient_local(i)*FE_solution_triangle(Gauss_point_local(i,1),Gauss_point_local(i,2),uh_local,vertices,0,0)*tri_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
end
