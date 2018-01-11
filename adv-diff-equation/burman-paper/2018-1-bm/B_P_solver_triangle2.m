function r=B_P_solver_triangle2(data,mesh_type,q,left,right,bottom,top,h_1,e,b,p,p_1,a_1,a_2)


N1=(right-left)/h_1(1);
N2=(top-bottom)/h_1(2);

%h=h_1(1)/sin(55/180*pi);
h=h_1(1);
c2=h^3/(2*b*h+h*sqrt(a_1^2+a_2^2)+e);

%Mesh information for partition and finite element basis functions.
if mesh_type==1
    [M,T]=mesh_divise_1(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
    matrix_size=[(N1+1)*(N2+1) (N1+1)*(N2+1)];
    vector_size=(N1+1)*(N2+1);
    [boundary_nodes,boundary_edges]=generate_boundary_nodes_edges(N1,N2,N1,N2);
elseif mesh_type==2
    [M,T]=mesh_divise_2(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
    matrix_size=[(N1+1)*(N2+1) (N1+1)*(N2+1)];
    vector_size=(N1+1)*(N2+1);
    [boundary_nodes,boundary_edges]=generate_boundary_nodes_edges_mesh2(T,N1,N2,N1,N2);
elseif mesh_type==3
    [M,T]=mesh_divise_3(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
    matrix_size=[(N1+1)*(N2+1) (N1+1)*(N2+1)];
    vector_size=(N1+1)*(N2+1);
    [boundary_nodes,boundary_edges]=generate_boundary_nodes_edges_mesh2(T,N1,N2,N1,N2);

elseif mesh_type==4
    [M,T]=mesh_divise_4(left,right,bottom,top,h_1);
    number_of_elements=4*N1*N2;
    matrix_size=[(N1+1)*(N2+1)+N1*N2 (N1+1)*(N2+1)+N1*N2];
    vector_size=(N1+1)*(N2+1)+N1*N2;
    [boundary_nodes,boundary_edges]=generate_boundary_nodes_edges_mesh3(N1,N2,N1,N2);
end

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
%A6=assemble_matrix(data.a,M,T,T,T,3,3,number_of_elements,matrix_size,G_coe_ref_tri,G_point_ref_tri,1,0,1,0);

A_1=e*(A1+A2)+A3*b+a_1*A4+a_2*A5;

disp('######## 40%')
B1=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,N1,N2,matrix_size,G_coe_ref_1D,G_point_ref_1D,1,0,1,0);
B2=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,N1,N2,matrix_size,G_coe_ref_1D,G_point_ref_1D,0,1,0,1);
B3=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,N1,N2,matrix_size,G_coe_ref_1D,G_point_ref_1D,0,0,0,0);
B4=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,N1,N2,matrix_size,G_coe_ref_1D,G_point_ref_1D,1,0,0,0);
B5=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,N1,N2,matrix_size,G_coe_ref_1D,G_point_ref_1D,0,1,0,0);
B6=assemble_out_line_matrix(data.a,M,T,boundary_edges,T,T,3,3,N1,N2,matrix_size,G_coe_ref_1D,G_point_ref_1D,1,0,0,1);
B=a_1^2*B1+a_2^2*B2+b^2*B3+a_1*a_2*(B6+B6')+a_1*b*(B4+B4')+a_2*b*(B5+B5');
BB=a_2^2*B1+a_1^2*B2-a_1*a_2*(B6'+B6);
%B=a_1^2*A1+a_2^2*A2+a_1*a_2*(A6+A6');
disp('############ 60%')

E=assemble_in_line_matrix_plus(q,data.a,Boundary,T_index,M,T,T_plus,3,3,matrix_size,G_coe_ref_1D,G_point_ref_1D,a_1,a_2,e);

A=A_1+E+p*c2*B+p_1*c2*BB;
% %Assemble the load vector.


k1=assemble_vector(data.f,M,T,T,3,number_of_elements,vector_size,G_coe_ref_tri,G_point_ref_tri,0,0);
k2=assemble_out_line_vector(data.f,M,T,boundary_edges,T,3,N1,N2,vector_size,G_coe_ref_1D,G_point_ref_1D,0,0);
k3=assemble_out_line_vector(data.f,M,T,boundary_edges,T,3,N1,N2,vector_size,G_coe_ref_1D,G_point_ref_1D,1,0);
k4=assemble_out_line_vector(data.f,M,T,boundary_edges,T,3,N1,N2,vector_size,G_coe_ref_1D,G_point_ref_1D,0,1);
k_1=a_1*k3+a_2*k4+b*k2;
disp('################ 80%')

k=k1+k_1*p*c2;
%Deal with Dirichlet boundary condition.
[A,k]=treat_Dirichlet_boundary_triangle(data.g,A,k,boundary_nodes,M);

%Compute the numerical solution
r=A\k;

end
