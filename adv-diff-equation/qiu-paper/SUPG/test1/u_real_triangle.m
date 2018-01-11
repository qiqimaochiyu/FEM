

function r=u_real_triangle(accurate_function,vertices,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)



Gpn=length(Gauss_coefficient_reference_triangle);

r=0;
[Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);

for i=1:Gpn
    r=r+Gauss_coefficient_local_triangle(i)*(feval(accurate_function,Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2)))^2;
end