function r=SUPG_solver_triangle_f(A,uh,F1,F2,F3,data,mesh_type,left,right,bottom,top,h_1,e,a_1,a_2)

global b

N1=(right-left)/h_1(1);
N2=(top-bottom)/h_1(2); 

if mesh_type==1
    [M,T]=mesh_divise_1(left,right,bottom,top,h_1);
    number_of_elements=2*N1*N2;
    matrix_size=[(N1+1)*(N2+1) (N1+1)*(N2+1)];
    vector_size=(N1+1)*(N2+1);
    [boundary_nodes,~]=generate_boundary_nodes_edges(N1,N2,N1,N2);
    h=h_1(1)*sqrt(2);
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


k=(F1+c*a_1*F2+c*a_2*F3)*uh;
%Deal with Dirichlet boundary condition.
[Ag,k]=treat_Dirichlet_boundary_triangle(data.g,A,k,boundary_nodes,M);

%Compute the numerical solution
r=Ag\k;

end