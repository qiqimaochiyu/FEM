function r=FE_solution_triangle(x,y,uh_local,vertices,derivative_degree_x,derivative_degree_y)

r=0;
number_of_local_basis=length(uh_local);
for i=1:number_of_local_basis
    r=r+uh_local(i)*tri_local_basis(x,y,vertices,i,derivative_degree_x,derivative_degree_y);
end