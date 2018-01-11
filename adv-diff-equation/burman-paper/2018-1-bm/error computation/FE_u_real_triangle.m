function r=FE_u_real_triangle(accurate_function,left,right,bottom,top,h_1,Gauss_point_number)

global mesh_type

N1_partition=(right-left)/h_1(1);
N2_partition=(top-bottom)/h_1(2);
if mesh_type==2    
    number_of_elements=2*N1_partition*N2_partition;

    [M_partition,T_partition]=mesh_divise2(left,right,bottom,top,h_1);
elseif mesh_type==3 
    number_of_elements=4*N1_partition*N2_partition;

    [M_partition,T_partition]=mesh_divise3(left,right,bottom,top,h_1);
end

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

r=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    vertices=M_partition(:,T_partition(:,n));
   
    r=r+u_real_triangle(accurate_function,vertices,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
end
r=sqrt(r);