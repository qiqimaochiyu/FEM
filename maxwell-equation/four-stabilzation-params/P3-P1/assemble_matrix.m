function r=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,degree_P1,trial_degree_x,trial_degree_y,degree_P2,test_degree_x,test_degree_y)


if degree_P1==1&&degree_P2==1
    r=sparse(tnp,tnp);
    T_basis_trial=T1;T_basis_test=T1;
    dof1=3;dof2=3;
elseif degree_P1==1&&degree_P2==3
    r=sparse(tnp_2,tnp);
    T_basis_trial=T1;T_basis_test=T;
    dof1=3;dof2=10;
elseif degree_P1==3&&degree_P2==1
    r=sparse(tnp,tnp_2);
    T_basis_trial=T;T_basis_test=T1;
    dof1=10;dof2=3;
elseif degree_P1==3&&degree_P2==3
    r=sparse(tnp_2,tnp_2);
    T_basis_trial=T;T_basis_test=T;
    dof1=10;dof2=10;    
end

    for n=1:number_of_elements   
        vertices=M(:,T(1:3,n));
        [Gauss_coefficient_local,Gauss_point_local]=generate_Gauss_local_triangle(vertices);
   
        for alpha=1:dof1
            for beta=1:dof2      
                temp=GQ_volume_integral_trial_test_tri(Gauss_coefficient_local,Gauss_point_local,vertices,alpha,degree_P1,trial_degree_x,trial_degree_y,beta,degree_P2,test_degree_x,test_degree_y);
                r(T_basis_test(beta,n),T_basis_trial(alpha,n))=r(T_basis_test(beta,n),T_basis_trial(alpha,n))+temp;
            end
        end
    end




