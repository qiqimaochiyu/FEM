function r=B_P_solver_triangle(data,q,e,b,b_p,a_1,a_2)


global p v t

%Mesh information for partition and finite element basis functions.
M=p;T=t(1:3,:);
number_of_elements=length(t);l=length(p);
matrix_size=[l l];vector_size=l;

[boundary_nodes,boundary_edges]=generate_edge(v,T);

%Guass quadrature's points and weights on the refenrece triangle and reference interval.
[G_coe_ref_tri,G_point_ref_tri]=generate_Gauss_reference_triangle(9);
[G_coe_ref_1D,G_point_ref_1D]=generate_Gauss_reference_1D(8);
[Boundary,T_plus,T_index]=generate_in_edges(T,number_of_elements);
disp('#### 20%')
%Assemble the stiffness matrix.

A1=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,1,0);
A2=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,0,1,0,1);
A3=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,0,0,0,0);
A4=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,0,0);
A5=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,0,1,0,0);
%A6=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,0,0,1,0);
%A6=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,0,1,1,0);
%A7=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,0,1);

A_1=e*(A1+A2)+A3*b+a_1*A4+a_2*A5;

disp('######## 40%')
B1=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,matrix_size,G_coe_ref_1D,G_point_ref_1D,1,0,1,0);
B2=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,matrix_size,G_coe_ref_1D,G_point_ref_1D,0,1,0,1);
B3=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,matrix_size,G_coe_ref_1D,G_point_ref_1D,0,0,0,0);
B4=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,matrix_size,G_coe_ref_1D,G_point_ref_1D,1,0,0,0);
B5=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,matrix_size,G_coe_ref_1D,G_point_ref_1D,0,1,0,0);
B6=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,matrix_size,G_coe_ref_1D,G_point_ref_1D,1,0,0,1);
B=a_1^2*B1+a_2^2*B2+b^2*B3+a_1*a_2*(B6+B6')+a_1*b*(B4+B4')+a_2*b*(B5+B5');
%B=a_1^2*A1+a_2^2*A2+a_1*a_2*(A6+A6');
disp('############ 60%')
E=assemble_in_line_matrix_plus(q,data.a,Boundary,T_index,M,T,T_plus,3,3,matrix_size,G_coe_ref_1D,G_point_ref_1D,a_1,a_2,e);

A=A_1+E+b_p*B;
% %Assemble the load vector.


k1=assemble_vector(data.f,M,T,T,3,number_of_elements,vector_size,G_coe_ref_tri,G_point_ref_tri,0,0);
k2=assemble_out_line_vector(data.f,M,T,boundary_edges,T,3,vector_size,G_coe_ref_1D,G_point_ref_1D,0,0);
k3=assemble_out_line_vector(data.f,M,T,boundary_edges,T,3,vector_size,G_coe_ref_1D,G_point_ref_1D,1,0);
k4=assemble_out_line_vector(data.f,M,T,boundary_edges,T,3,vector_size,G_coe_ref_1D,G_point_ref_1D,0,1);
k_1=a_1*k3+a_2*k4+b*k2;


disp('################ 80%')
k=k1+k_1*b_p;
%Deal with Dirichlet boundary condition.
[A,k]=treat_Dirichlet_boundary_triangle(data.g,A,k,boundary_nodes,M);

%Compute the numerical solution
r=A\k;

end

