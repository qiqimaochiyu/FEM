function r=assemble_out_line_matrix(index,coefficient_function_name,M,T,T1,tnp,tnp_2,number_of_edges,boundary_edges,degree_P1,trial_degree_x,trial_degree_y,degree_P2,test_degree_x,test_degree_y)

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
    
    vertices=M(:,T(1:3,boundary_edges(2,n)));
    end_point_1=M(:,boundary_edges(3,n));
    end_point_2=M(:,boundary_edges(4,n));
    v=[end_point_2(1)-end_point_1(1) end_point_2(2)-end_point_1(2)];
    V=v/norm(v);
    N=V;
    if index==1
        c=N(1);
    elseif index==2
        c=N(2);
    elseif index==3
        c=N(1)^2;
    elseif index==4
        c=N(2)^2;
    elseif index==5
        c=N(1)*N(2);
    end
    
    for alpha=1:dof1
       for beta=1:dof2
            temp=c*GQ_line_integral_trial_test_tri(coefficient_function_name,end_point_1,end_point_2,vertices,degree_P1,alpha,trial_degree_x,trial_degree_y,vertices,degree_P2,beta,test_degree_x,test_degree_y);
            r(T_basis_test(beta,boundary_edges(2,n)),T_basis_trial(alpha,boundary_edges(2,n)))=r(T_basis_test(beta,boundary_edges(2,n)),T_basis_trial(alpha,boundary_edges(2,n)))+temp;
       end
    end
end