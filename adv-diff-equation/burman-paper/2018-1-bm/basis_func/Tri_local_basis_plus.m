function r=Tri_local_basis_plus(T_index,lt,vertices_1,vertices_2,x,y,index,derivative_degree_x,derivative_degree_y)


if index==1
    r=tri_local_basis(x,y,vertices_1,T_index(1,lt),derivative_degree_x,derivative_degree_y)-tri_local_basis(x,y,vertices_2,T_index(2,lt),derivative_degree_x,derivative_degree_y);
elseif index==2
    r=tri_local_basis(x,y,vertices_1,T_index(3,lt),derivative_degree_x,derivative_degree_y)-tri_local_basis(x,y,vertices_2,T_index(4,lt),derivative_degree_x,derivative_degree_y);
elseif index==3
    r=tri_local_basis(x,y,vertices_1,T_index(5,lt),derivative_degree_x,derivative_degree_y);
elseif index==4
    r=-tri_local_basis(x,y,vertices_2,T_index(6,lt),derivative_degree_x,derivative_degree_y);
end
