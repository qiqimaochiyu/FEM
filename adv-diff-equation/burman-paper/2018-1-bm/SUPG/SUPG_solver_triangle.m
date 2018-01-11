function r=SUPG_solver_triangle(data,mesh_type,left,right,bottom,top,h_1,e,a_1,a_2)
%Yu wei 2018-01



N1=(right-left)/h_1(1);
N2=(top-bottom)/h_1(2); 

if mesh_type==1
    [M,T]=mesh_divise_1(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
    matrix_size=[(N1+1)*(N2+1) (N1+1)*(N2+1)];
    vector_size=(N1+1)*(N2+1);
    [boundary_nodes,~]=generate_boundary_nodes_edges(N1,N2,N1,N2);
    h=h_1(1);
    c=h*(coth(h/(2*e))-2*e/h)/2;
elseif mesh_type==2
    [M,T]=mesh_divise_2(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
    matrix_size=[(N1+1)*(N2+1) (N1+1)*(N2+1)];
    vector_size=(N1+1)*(N2+1);
    [boundary_nodes,~]=generate_boundary_nodes_edges_mesh2(T,N1,N2,N1,N2);
    h=h_1(1);
    c=h*(coth(h/(2*e))-2*e/h)/2;
elseif mesh_type==3
    [M,T]=mesh_divise_3(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
    matrix_size=[(N1+1)*(N2+1) (N1+1)*(N2+1)];
    vector_size=(N1+1)*(N2+1);
    [boundary_nodes,~]=generate_boundary_nodes_edges_mesh2(T,N1,N2,N1,N2);

elseif mesh_type==4
    [M,T]=mesh_divise_4(left,right,bottom,top,h_1);
    number_of_elements=4*N1*N2;
    matrix_size=[(N1+1)*(N2+1)+N1*N2 (N1+1)*(N2+1)+N1*N2];
    vector_size=(N1+1)*(N2+1)+N1*N2;
    [boundary_nodes,~]=generate_boundary_nodes_edges_mesh3(N1,N2,N1,N2);
end


%Guass quadrature's points and weights on the refenrece triangle and reference interval.
[G_coe_ref_tri,G_point_ref_tri]=generate_Gauss_reference_triangle(9);
%The following line is necessary if we need to deal with Neumann or Robin boundary condition. 
%Otherwise, we can comment this line to save time and memory.
%[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(4);

%Assemble the stiffness matrix.
A1=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,1,0);
A2=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,0,1,0,1);
A4=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,0,0);
A6=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,0,1);
A7=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,0,1,0,0);
A9=A6';
B1=assemble_matrix_SUPG(mesh_type,e,h_1,data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,1,0);
B2=assemble_matrix_SUPG(mesh_type,e,h_1,data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,0,1,0,1);
B6=assemble_matrix_SUPG(mesh_type,e,h_1,data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,0,1);
B9=B6';

%Assemble the load vector.
k1=assemble_vector(data.f,M,T,T,3,number_of_elements,vector_size,G_coe_ref_tri,G_point_ref_tri,0,0);
k2=assemble_vector(data.f,M,T,T,3,number_of_elements,vector_size,G_coe_ref_tri,G_point_ref_tri,1,0);
k3=assemble_vector(data.f,M,T,T,3,number_of_elements,vector_size,G_coe_ref_tri,G_point_ref_tri,0,1);
s2=assemble_vector_SUPG(mesh_type,e,h_1,data.f,M,T,T,3,number_of_elements,vector_size,G_coe_ref_tri,G_point_ref_tri,1,0);
s3=assemble_vector_SUPG(mesh_type,e,h_1,data.f,M,T,T,3,number_of_elements,vector_size,G_coe_ref_tri,G_point_ref_tri,0,1);

if mesh_type==1 || mesh_type==2
    A=(c*a_1^2+e)*A1+(c*a_2^2+e)*A2+c*a_1*a_2*(A6+A9)+a_1*A4+a_2*A7;
    k=k1+c*a_1*k2+c*a_2*k3;
elseif mesh_type==3 || mesh_type==4
    A=a_1^2*B1+e*A1+a_2^2*B2+e*A2+a_1*a_2*(B6+B9)+a_1*A4+a_2*A7;
    k=k1+a_1*s2+a_2*s3;
end


    


%Get the information matrices for boundary nodes and boundary edges.

%Deal with Neumann boundary condition. If we don't have Neumann boundary condition at all, we can comment the following line to save time and memory.
%b=treat_Neumann_boundary_triangle('function_q_tilde',b,boundary_edges,M,T,T,3,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,0);

%Deal with Robin boundary condition. If we don't have Robin boundary condition at all, we can comment the following line to save time and memory.
%[A,b]=treat_Robin_boundary_triangle('function_q','function_p',A,b,boundary_edges,M,T,T,T,3,3,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,0,basis_type,0,0);

%Deal with Dirichlet boundary condition.
[A,k]=treat_Dirichlet_boundary_triangle(data.g,A,k,boundary_nodes,M);

%Compute the numerical solution
r=A\k;






