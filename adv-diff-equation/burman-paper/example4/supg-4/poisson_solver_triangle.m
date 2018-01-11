

load data2
e=10^(-4);
h=0.32;
c=h*(coth(h/(2*e))-2*e/h)/2;
a_1=1;a_2=0;

 M=p;
 T=t(1:3,:);
 M_partition=M;
 M_basis=M_partition;
 T_partition=T;
 T_basis=T_partition;
%Guass quadrature's points and weights on the refenrece triangle and reference interval.
[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(9);

%Assemble the stiffness matrix.
number_of_elements=length(t);
l=length(p);
matrix_size=[l l];

    number_of_trial_local_basis=3;
    number_of_test_local_basis=3;

A1=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,1,0);
A2=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1,0,1);
A4=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,0);
A6=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1);
A7=assemble_matrix('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1,0,0);
A9=A6';

A=(c*a_1^2+e)*A1+(c*a_2^2+e)*A2+c*a_1*a_2*(A6+A9)+a_1*A4+a_2*A7;
%Assemble the load vector.
vector_size=l;
k1=assemble_vector('function_f',M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,0);
k2=assemble_vector('function_f',M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0);
k3=assemble_vector('function_f',M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,0,1);

k=k1+c*a_1*k2+c*a_2*k3;


%Get the information matrices for boundary nodes and boundary edges.
% [boundary_nodes,boundary_edges]=generate_boundary_nodes_edges(N1_basis,N2_basis,N1_partition,N2_partition);
boundary_nodes=[v(5,:);v(1,:)];
boundary_edges=[v(5,:);v(1,:);v(2,:)];
% Deal with Neumann boundary condition. If we don't have Neumann boundary condition at all, we can comment the following line to save time and memory.

[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(8);

% k=treat_Neumann_boundary_triangle('function_q',b,boundary_edges,M_partition,T_partition,T_basis,number_of_test_local_basis,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,1,0,0);

%Deal with Robin boundary condition. If we don't have Robin boundary condition at all, we can comment the following line to save time and memory.
%[A,b]=treat_Robin_boundary_triangle('function_q','function_p',A,b,boundary_edges,M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,0,basis_type,0,0);


%Deal with Dirichlet boundary condition.
[A,k]=treat_Dirichlet_boundary_triangle('function_g',A,k,boundary_nodes,M_basis);

%Compute the numerical solution
r=A\k;
x1=p(1,:);
x2=p(2,:);



xlin = linspace(-3,9,2401);
ylin = linspace(-3,3,1201);
[X,Y]=meshgrid(xlin,ylin);
s = X.^2+Y.^2;
flag = find( s<1 );

Z = griddata(x1,x2,r',X,Y,'cubic');
Z(flag)=nan;

colormap('jet');
mesh(X,Y,Z);
colorbar;
figure
y11=Z(801,:);
plot(xlin,y11,'b');

figure
yl=linspace(0,3,601);
x11=Z(601:1201,1301)';
plot(yl,x11,'b');


save supg-result-data4


 


