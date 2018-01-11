function r=assemble_matrix(coefficient_function_name,M_partition,T_partition,T_basis_trial,T_basis_test,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,trial_derivative_degree_x,trial_derivative_degree_y,test_derivative_degree_x,test_derivative_degree_y)


r=sparse(matrix_size(1),matrix_size(2));

for n=1:number_of_elements   
    vertices=M_partition(:,T_partition(1:3,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
   
    for alpha=1:number_of_trial_local_basis
       for beta=1:number_of_test_local_basis      
          temp=Gauss_quadrature_for_volume_integral_trial_test_triangle(coefficient_function_name,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,alpha,trial_derivative_degree_x,trial_derivative_degree_y,beta,test_derivative_degree_x,test_derivative_degree_y);
          r(T_basis_test(beta,n),T_basis_trial(alpha,n))=r(T_basis_test(beta,n),T_basis_trial(alpha,n))+temp;
       end
    end
end


