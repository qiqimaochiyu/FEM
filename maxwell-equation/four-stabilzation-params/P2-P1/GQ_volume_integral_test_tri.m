function r=GQ_volume_integral_test_tri(accurate_function,Gauss_coefficient_local,Gauss_point_local,vertices,test_basis_index,test_degree_x,test_degree_y,degree_P)
%JYu wei 2017-05
%Gpn: the number of the Gauss points of the Gauss quadrature we are using.

Gpn=length(Gauss_coefficient_local);
r=0;
for i=1:Gpn
      fh=feval(accurate_function,Gauss_point_local(i,1),Gauss_point_local(i,2));
      r=r+Gauss_coefficient_local(i)*fh*tri_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,test_basis_index,test_degree_x,test_degree_y,degree_P);
end
