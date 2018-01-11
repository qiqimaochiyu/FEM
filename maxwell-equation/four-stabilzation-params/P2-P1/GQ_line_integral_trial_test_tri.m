function r=GQ_line_integral_trial_test_tri(coefficient_function_name,end_point_1,end_point_2,trial_vertices,trial_basis_type,trial_basis_index,trial_degree_x,trial_degree_y,test_vertices,test_basis_type,test_basis_index,test_degree_x,test_degree_y)
% Yu wei 2017-05
%Use Gauss quadrature to compute a line integral on an edge for a matrix.
%Gpn: the number of the Gauss points of the Gauss formula we are using.
%Gauss_coefficient_local_1D,Gauss_point_local_1D:the Gauss coefficients and Gauss points on the local interval.
global Gauss_coefficient_reference_1D
Gpn=length(Gauss_coefficient_reference_1D);
r=0;

if end_point_1(2)==end_point_2(2)
%The line is horizontal.

    lower_bound=min(end_point_1(1),end_point_2(1));
    upper_bound=max(end_point_1(1),end_point_2(1));
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(lower_bound,upper_bound);
    for i=1:Gpn
         r=r+Gauss_coefficient_local_1D(i)*feval(coefficient_function_name,Gauss_point_local_1D(i),end_point_1(2))*tri_local_basis(Gauss_point_local_1D(i),end_point_1(2),trial_vertices,trial_basis_index,trial_degree_x,trial_degree_y,trial_basis_type)*tri_local_basis(Gauss_point_local_1D(i),end_point_1(2),test_vertices,test_basis_index,test_degree_x,test_degree_y,test_basis_type);
    end    

elseif end_point_1(1)==end_point_2(1)
%The line is vertical.

    lower_bound=min(end_point_1(2),end_point_2(2));
    upper_bound=max(end_point_1(2),end_point_2(2));
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(lower_bound,upper_bound);
    for i=1:Gpn
         r=r+Gauss_coefficient_local_1D(i)*feval(coefficient_function_name,end_point_1(1),Gauss_point_local_1D(i))*tri_local_basis(end_point_1(1),Gauss_point_local_1D(i),trial_vertices,trial_basis_index,trial_degree_x,trial_degree_y,trial_basis_type)*tri_local_basis(end_point_1(1),Gauss_point_local_1D(i),test_vertices,test_basis_index,test_degree_x,test_degree_y,test_basis_type);
    end
    
else
%The slope of the edge is in (0,infinity).

    lower_bound=min(end_point_1(1),end_point_2(1));
    upper_bound=max(end_point_1(1),end_point_2(1));
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(lower_bound,upper_bound);
    slope=(end_point_2(2)-end_point_1(2))/(end_point_2(1)-end_point_1(1));
    Jacobi=sqrt(1+slope^2);
    for i=1:Gpn
         x=Gauss_point_local_1D(i);
         y=slope*(x-end_point_1(1))+end_point_1(2);
         r=r+Gauss_coefficient_local_1D(i)*Jacobi*feval(coefficient_function_name,x,y)*tri_local_basis(x,y,trial_vertices,trial_basis_index,trial_degree_x,trial_degree_y,trial_basis_type)*tri_local_basis(x,y,test_vertices,test_basis_index,test_degree_x,test_degree_y,test_basis_type);
    end

end