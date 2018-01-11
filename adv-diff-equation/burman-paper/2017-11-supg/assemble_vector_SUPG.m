function r=assemble_vector_SUPG(e,h_1,coefficient_function_name,M_partition,T_partition,T_basis_test,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,test_derivative_degree_x,test_derivative_degree_y)
global mesh_type
r=zeros(vector_size,1);

if mesh_type==3

    for n=1:number_of_elements

        vertices=M_partition(:,T_partition(:,n));
        [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
        if mod(n,4)==0 || mod(n,4)==2
            h=sqrt(1/2)*h_1(1)/cos(15/180*pi);
        else
            h=1/2/sin(55*pi/180)*h_1(1);
        end
            c=h*(coth(h/(2*e))-2*e/h)/2;
            for beta=1:number_of_test_local_basis     
            temp=c*Gauss_quadrature_for_volume_integral_test_triangle(coefficient_function_name,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,beta,test_derivative_degree_x,test_derivative_degree_y);
            r(T_basis_test(beta,n),1)=r(T_basis_test(beta,n),1)+temp;
            end
    end
elseif mesh_type==2
    

    for n=1:number_of_elements

        vertices=M_partition(:,T_partition(1:3,n));
        [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
        if T_partition(4,n)==-1
            h=sqrt(1/2)*h_1(1)/cos(15/180*pi);
        else
            h=h_1(1)/sin(55*pi/180);
        end
            c=h*(coth(h/(2*e))-2*e/h)/2;
            for beta=1:number_of_test_local_basis     
            temp=c*Gauss_quadrature_for_volume_integral_test_triangle(coefficient_function_name,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,beta,test_derivative_degree_x,test_derivative_degree_y);
            r(T_basis_test(beta,n),1)=r(T_basis_test(beta,n),1)+temp;
            end
    end
end