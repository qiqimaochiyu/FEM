function r=assemble_in_line_matrix_plus(q,coefficient_function_name,Boundary,T_index,M_basis,T_basis,T_plus,number_of_trial_local_basis,number_of_test_local_basis,matrix_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,a_1,a_2,e)
 
global b h_1
r=sparse(matrix_size(1),matrix_size(2));


for n=1:length(Boundary)
    vertices_1=M_basis(:,T_basis(1:3,Boundary(1,n)));
    vertices_2=M_basis(:,T_basis(1:3,Boundary(2,n)));
    end_point_1=M_basis(:,T_plus(1,n));
    end_point_2=M_basis(:,T_plus(2,n));
    end_point_3=M_basis(:,T_plus(3,n));
    end_point_4=M_basis(:,T_plus(4,n));
    v=[-(end_point_1(2)-end_point_2(2)) end_point_1(1)-end_point_2(1)];
    k=[end_point_4(1)-end_point_3(1) end_point_4(2)-end_point_3(2)];
    h=norm(v);
    % c=q*(b*h_1(1)^2+h^3/(h*sqrt(a_1^2+a_2^2)+e));
    c=q*(h^3/(h*sqrt(a_1^2+a_2^2)+e));
    if v*k'>=0
        nv=v/h;
    else
        nv=-v/h;
        
    end
    for alpha=1:(number_of_trial_local_basis+1)
       for beta=1:(number_of_test_local_basis+1)
            temp=(nv(1)^2)*Gauss_quadrature_for_line_integral_trial_test_triangle_plus(T_index,n,vertices_1,vertices_2,coefficient_function_name,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2,alpha,1,0,beta,1,0)...
            +(nv(2)^2)*Gauss_quadrature_for_line_integral_trial_test_triangle_plus(T_index,n,vertices_1,vertices_2,coefficient_function_name,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2,alpha,0,1,beta,0,1)...
            +(nv(1)*nv(2))*Gauss_quadrature_for_line_integral_trial_test_triangle_plus(T_index,n,vertices_1,vertices_2,coefficient_function_name,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2,alpha,1,0,beta,0,1)...
            +(nv(1)*nv(2))*Gauss_quadrature_for_line_integral_trial_test_triangle_plus(T_index,n,vertices_1,vertices_2,coefficient_function_name,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2,alpha,0,1,beta,1,0);

            r(T_plus(beta,n),T_plus(alpha,n))=r(T_plus(beta,n),T_plus(alpha,n))+temp*c;
       end
    end
   
end