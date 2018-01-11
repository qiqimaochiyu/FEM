function r=assemble_out_line_matrix(coefficient_function_name,M_basis,T_basis,boundary_edges,T_basis_trial,T_basis_test,number_of_trial_local_basis,number_of_test_local_basis,N1_partition,N2_partition,matrix_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,trial_derivative_degree_x,trial_derivative_degree_y,test_derivative_degree_x,test_derivative_degree_y)

r=sparse(matrix_size(1),matrix_size(2));

for n=1:2*(N1_partition+N2_partition)
    
    vertices=M_basis(:,T_basis(1:3,boundary_edges(2,n)));
    end_point_1=M_basis(:,boundary_edges(3,n));
    end_point_2=M_basis(:,boundary_edges(4,n));
    
    for alpha=1:number_of_trial_local_basis
       for beta=1:number_of_test_local_basis 
            temp=Gauss_quadrature_for_line_integral_trial_test_triangle(coefficient_function_name,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2,vertices,alpha,trial_derivative_degree_x,trial_derivative_degree_y,beta,test_derivative_degree_x,test_derivative_degree_y);
            r(T_basis_test(beta,boundary_edges(2,n)),T_basis_trial(alpha,boundary_edges(2,n)))=r(T_basis_test(beta,boundary_edges(2,n)),T_basis_trial(alpha,boundary_edges(2,n)))+temp;
       end
    end
end