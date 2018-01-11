function r=Gauss_quadrature_for_line_integral_trial_test_triangle_plus(T_index,lt,vertices_1,vertices_2,coefficient_function_name,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2,trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_index,test_derivative_degree_x,test_derivative_degree_y)


Gpn=length(Gauss_coefficient_reference_1D);
r=0;

if end_point_1(2)==end_point_2(2)
%The line is horizontal.

    lower_bound=min(end_point_1(1),end_point_2(1));
    upper_bound=max(end_point_1(1),end_point_2(1));
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    for i=1:Gpn
         r=r+Gauss_coefficient_local_1D(i)*feval(coefficient_function_name,Gauss_point_local_1D(i),end_point_1(2))*Tri_local_basis_plus(T_index,lt,vertices_1,vertices_2,Gauss_point_local_1D(i),end_point_1(2),trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y)*Tri_local_basis_plus(T_index,lt,vertices_1,vertices_2,Gauss_point_local_1D(i),end_point_1(2),test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
    end    

elseif end_point_1(1)==end_point_2(1)
%The line is vertical.

    lower_bound=min(end_point_1(2),end_point_2(2));
    upper_bound=max(end_point_1(2),end_point_2(2));
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    for i=1:Gpn
         r=r+Gauss_coefficient_local_1D(i)*feval(coefficient_function_name,end_point_1(1),Gauss_point_local_1D(i))*Tri_local_basis_plus(T_index,lt,vertices_1,vertices_2,end_point_1(1),Gauss_point_local_1D(i),trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y)*Tri_local_basis_plus(T_index,lt,vertices_1,vertices_2,end_point_1(1),Gauss_point_local_1D(i),test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
    end
    
else
%The slope of the edge is in (0,infinity).

    lower_bound=min(end_point_1(1),end_point_2(1));
    upper_bound=max(end_point_1(1),end_point_2(1));
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    slope=(end_point_2(2)-end_point_1(2))/(end_point_2(1)-end_point_1(1));
    Jacobi=sqrt(1+slope^2);
    for i=1:Gpn
         x=Gauss_point_local_1D(i);
         y=slope*(x-end_point_1(1))+end_point_1(2);
         r=r+Gauss_coefficient_local_1D(i)*Jacobi*feval(coefficient_function_name,x,y)*Tri_local_basis_plus(T_index,lt,vertices_1,vertices_2,x,y,trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y)*Tri_local_basis_plus(T_index,lt,vertices_1,vertices_2,x,y,test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
    end

end