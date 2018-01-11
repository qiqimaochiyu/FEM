function r=poisson_solver_triangle(left,right,bottom,top,h_1,e,b,a_1,a_2)


h=sqrt(4/3)*h_1(1);

if 1/(b*h)<=1
    xi=0;
else
    if 1/(b*h)>1 && h/e<1
        xi=1;
    else
        if 1/(b*h)>1 && h/e>=1
            xi=7*e/h;
        end
    end
end
c=h^2/(b*h^2+h*xi+6*e);

N1_partition=(right-left)/h_1(1);
N2_partition=(top-bottom)/h_1(2);

    N1_basis=N1_partition;
    N2_basis=N2_partition;   

%Mesh information for partition and finite element basis functions.
[M_partition,T_partition]=mesh_divise(left,right,bottom,top,h_1);

   M_basis=M_partition;
   T_basis=T_partition;

%Guass quadrature's points and weights on the refenrece triangle and reference interval.
[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(9);
%The following line is necessary if we need to deal with Neumann or Robin boundary condition. 
%Otherwise, we can comment this line to save time and memory.
%[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(4);

%Assemble the stiffness matrix.
number_of_elements=2*N1_partition*N2_partition;
matrix_size=[(N1_basis+1)*(N2_basis+1) (N1_basis+1)*(N2_basis+1)];

    number_of_trial_local_basis=3;
    number_of_test_local_basis=3;

A1=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,1,0);
A2=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1,0,1);
A3=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,0,0,0);
A4=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,0);
A5=A4';
% A5=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,0,1,0);
A6=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1);
A7=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1,0,0);
A8=A7';
% A8=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,0,0,1);
A9=A6';
A=(c*xi*a_1^2+e)*A1+(c*xi*a_2^2+e)*A2+c*xi*a_1*a_2*(A6+A9)+(a_1-c*a_1*b)*A4+(a_2-c*a_2*b)*A7+c*xi*b*a_1*A5+c*xi*b*a_2*A8+(b-c*b^2)*A3;
%Assemble the load vector.
vector_size=(N1_basis+1)*(N2_basis+1);
k1=assemble_vector('function_f',M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,0,a_1,a_2,e,b);
k2=assemble_vector('function_f',M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,a_1,a_2,e,b);
k3=assemble_vector('function_f',M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1,a_1,a_2,e,b);

k=(1-b*c)*k1+c*xi*a_1*k2+c*xi*a_2*k3;


%Get the information matrices for boundary nodes and boundary edges.
[boundary_nodes,boundary_edges]=generate_boundary_nodes_edges(N1_basis,N2_basis,N1_partition,N2_partition);

%Deal with Neumann boundary condition. If we don't have Neumann boundary condition at all, we can comment the following line to save time and memory.
%b=treat_Neumann_boundary_triangle('function_q_tilde',b,boundary_edges,M_partition,T_partition,T_basis,number_of_test_local_basis,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,0);

%Deal with Robin boundary condition. If we don't have Robin boundary condition at all, we can comment the following line to save time and memory.
%[A,b]=treat_Robin_boundary_triangle('function_q','function_p',A,b,boundary_edges,M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,0,basis_type,0,0);

%Deal with Dirichlet boundary condition.
[A,k]=treat_Dirichlet_boundary_triangle('function_g',A,k,boundary_nodes,M_basis);

%Compute the numerical solution
r=A\k;






