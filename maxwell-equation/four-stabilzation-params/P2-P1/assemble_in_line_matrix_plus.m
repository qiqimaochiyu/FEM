function r=assemble_in_line_matrix_plus(coefficient_function_name,M,T,T_index,T1,tnp,tnp_2,number_of_edges,boundary_edges,degree_P1,trial_degree_x,trial_degree_y,degree_P2,test_degree_x,test_degree_y)
 

if degree_P1==1&&degree_P2==1
    r=sparse(tnp,tnp);
    T_basis_trial=T1;T_basis_test=T1;
    dof1=3;dof2=3;
elseif degree_P1==1&&degree_P2==2
    r=sparse(tnp_2,tnp);
    T_basis_trial=T1;T_basis_test=T;
    dof1=3;dof2=6;
elseif degree_P1==2&&degree_P2==1
    r=sparse(tnp,tnp_2);
    T_basis_trial=T;T_basis_test=T1;
    dof1=6;dof2=3;
elseif degree_P1==2&&degree_P2==2
    r=sparse(tnp_2,tnp_2);
    T_basis_trial=T;T_basis_test=T;
    dof1=6;dof2=6;    
end


for n=1:number_of_edges
    vertices_1=M(:,T(1:3,boundary_edges(1,n)));
    vertices_2=M(:,T(1:3,boundary_edges(2,n)));
    end_point_1=M(:,T1(1,n));
    end_point_2=M(:,T1(2,n));
    v=[-(end_point_1(2)-end_point_2(2)) end_point_1(1)-end_point_2(1)];
    h=norm(v);
    for alpha=1:(dof1+1)
       for beta=1:(dof2+1)
            temp=h*GQ_line_integral_trial_test_tri_plus(T_index,n,vertices_1,vertices_2,coefficient_function_name,end_point_1,end_point_2,alpha,trial_degree_x,trial_degree_y,beta,test_degree_x,test_degree_y);
            r(T_basis_test(beta,n),T_basis_trial(alpha,n))=r(T_basis_test(beta,n),T_basis_trial(alpha,n))+temp;
       end
    end
   
end