function r=GQ_volume_integral_trial_test_tri(Gauss_coefficient_local,Gauss_point_local,vertices,trial_basis_index,degree_P1,trial_degree_x,trial_degree_y,test_basis_index,degree_P2,test_degree_x,test_degree_y)
%Yu wei 2017-05
%Gpn: the number of the Gauss points of the Gauss quadrature we are using.

Gpn=length(Gauss_coefficient_local);
r=0;
for i=1:Gpn
    r=r+Gauss_coefficient_local(i)*tri_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,trial_basis_index,trial_degree_x,trial_degree_y,degree_P1)*tri_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,test_basis_index,test_degree_x,test_degree_y,degree_P2);
end