function r=assemble_vector(accurate_function,M,T,T1,tnp,tnp_2,number_of_elements,degree_P,test_degree_x,test_degree_y)

if degree_P==1
    r=zeros(tnp,1);
    dof=3;
    T_basis_test=T1;
elseif degree_P==3
    r=zeros(tnp_2,1);
    dof=10;
    T_basis_test=T;
end

for n=1:number_of_elements
    vertices=M(:,T(1:3,n));
    [Gauss_coefficient_local,Gauss_point_local]=generate_Gauss_local_triangle(vertices);
    for beta=1:dof    
       temp=GQ_volume_integral_test_tri(accurate_function,Gauss_coefficient_local,Gauss_point_local,vertices,beta,test_degree_x,test_degree_y,degree_P);
       r(T_basis_test(beta,n),1)=r(T_basis_test(beta,n),1)+temp;
    end

end