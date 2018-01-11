function r=tri_local_basis(x,y,vertices,basis_index,derivative_degree_x,derivative_degree_y)

J_11=vertices(1,2)-vertices(1,1);
J_12=vertices(1,3)-vertices(1,1);
J_21=vertices(2,2)-vertices(2,1);
J_22=vertices(2,3)-vertices(2,1);
J_det=J_11*J_22-J_12*J_21;

x_hat=(J_22*(x-vertices(1,1))-J_12*(y-vertices(2,1)))/J_det;
y_hat=(-J_21*(x-vertices(1,1))+J_11*(y-vertices(2,1)))/J_det;

if derivative_degree_x==0&&derivative_degree_y==0
    r=tri_reference_basis(x_hat,y_hat,basis_index,0,0);
elseif derivative_degree_x==1&&derivative_degree_y==0
    r=(tri_reference_basis(x_hat,y_hat,basis_index,1,0)*J_22+tri_reference_basis(x_hat,y_hat,basis_index,0,1)*(-J_21))/J_det;
elseif derivative_degree_x==0&&derivative_degree_y==1
    r=(tri_reference_basis(x_hat,y_hat,basis_index,1,0)*(-J_12)+tri_reference_basis(x_hat,y_hat,basis_index,0,1)*J_11)/J_det;
elseif derivative_degree_x==2&&derivative_degree_y==0
    r=(tri_reference_basis(x_hat,y_hat,basis_index,2,0)*J_22^2+tri_reference_basis(x_hat,y_hat,basis_index,0,2)*J_21^2+tri_reference_basis(x_hat,y_hat,basis_index,1,1)*(-2*J_21*J_22))/J_det^2;
elseif derivative_degree_x==0&&derivative_degree_y==2
    r=(tri_reference_basis(x_hat,y_hat,basis_index,2,0)*J_12^2+tri_reference_basis(x_hat,y_hat,basis_index,0,2)*J_11^2+tri_reference_basis(x_hat,y_hat,basis_index,1,1)*(-2*J_11*J_12))/J_det^2;
elseif derivative_degree_x==1&&derivative_degree_y==1
    r=(tri_reference_basis(x_hat,y_hat,basis_index,2,0)*(-J_22*J_12)+tri_reference_basis(x_hat,y_hat,basis_index,0,2)*(-J_21*J_11)+tri_reference_basis(x_hat,y_hat,basis_index,1,1)*(J_21*J_12+J_11*J_22))/J_det^2;
end