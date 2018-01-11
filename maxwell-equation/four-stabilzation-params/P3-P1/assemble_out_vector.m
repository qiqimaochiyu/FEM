function r=assemble_out_vector(coefficient_function_name,M,T,T1,boundary_edges,tnp,tnp_2,number_of_edges,degree_P,test_degree_x,test_degree_y)


if degree_P==1
    r=zeros(tnp,1);
    dof=3;
    T_basis_test=T1;
elseif degree_P==3
    r=zeros(tnp_2,1);
    dof=10;
    T_basis_test=T;
end

for n=1:number_of_edges
    
    vertices=M(:,T(1:3,boundary_edges(2,n)));
    end_point_1=M(:,boundary_edges(3,n));
    end_point_2=M(:,boundary_edges(4,n));

    for beta=1:dof 
            temp=GQ_line_integral_test_tri(coefficient_function_name,end_point_1,end_point_2,vertices,beta,degree_P,test_degree_x,test_degree_y);
            r(T_basis_test(beta,boundary_edges(2,n)),1)=r(T_basis_test(beta,boundary_edges(2,n)),1)+temp;
      
    end
end




