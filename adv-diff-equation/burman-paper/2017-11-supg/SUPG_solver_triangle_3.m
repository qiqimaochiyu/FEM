function r=SUPG_solver_triangle_3(left,right,bottom,top,h_1,e,a_1,a_2)

data=example1;

N1=(right-left)/h_1(1);
N2=(top-bottom)/h_1(2);
%Mesh information for partition and finite element basis functions.
[M,T]=mesh_divise3(left,right,bottom,top,h_1);

%Guass quadrature's points and weights on the refenrece triangle and reference interval.
[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(9);
%The following line is necessary if we need to deal with Neumann or Robin boundary condition. 
%Otherwise, we can comment this line to save time and memory.
%[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(4);

%Assemble the stiffness matrix.
number_of_elements=4*N1*N2;
matrix_size=[(N1+1)*(N2+1)+N1*N2 (N1+1)*(N2+1)+N1*N2];

    number_of_trial_local_basis=3;
    number_of_test_local_basis=3;

A1=assemble_matrix(data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,1,0);
A2=assemble_matrix(data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1,0,1);
A4=assemble_matrix(data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,0);
% A6=assemble_matrix(data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1);
A7=assemble_matrix(data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1,0,0);
% A9=A6';
% A3=assemble_matrix(data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,0,0,0);
B1=assemble_matrix_SUPG(e,h_1,data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,1,0);
B2=assemble_matrix_SUPG(e,h_1,data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1,0,1);
B6=assemble_matrix_SUPG(e,h_1,data.a,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1);
B9=B6';
A=a_1^2*B1+e*A1+a_2^2*B2+e*A2+a_1*a_2*(B6+B9)+a_1*A4+a_2*A7;
%Assemble the load vector.
vector_size=(N1+1)*(N2+1)+N1*N2;
k1=assemble_vector(data.f,M,T,T,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,0);
k2=assemble_vector_SUPG(e,h_1,data.f,M,T,T,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0);
k3=assemble_vector_SUPG(e,h_1,data.f,M,T,T,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1);

k=k1+a_1*k2+a_2*k3;


%Get the information matrices for boundary nodes and boundary edges.
[boundary_nodes,~]=generate_boundary_nodes_edges(N1,N2,N1,N2);

%Deal with Neumann boundary condition. If we don't have Neumann boundary condition at all, we can comment the following line to save time and memory.
%b=treat_Neumann_boundary_triangle('function_q_tilde',b,boundary_edges,M,T,T,number_of_test_local_basis,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,0);

%Deal with Robin boundary condition. If we don't have Robin boundary condition at all, we can comment the following line to save time and memory.
%[A,b]=treat_Robin_boundary_triangle('function_q','function_p',A,b,boundary_edges,M,T,T,T,number_of_trial_local_basis,number_of_test_local_basis,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,0,basis_type,0,0);

%Deal with Dirichlet boundary condition.
[A,k]=treat_Dirichlet_boundary_triangle(data.g,A,k,boundary_nodes,M);

%Compute the numerical solution
r=A\k;






