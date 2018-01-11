function [A,b]=treat_Dirichlet_boundary_triangle(Dirichlet_boundary_function_name,A,b,boundary_nodes,M)


nbn=size(boundary_nodes,2);

%Check all boundary nodes of FE.
for k=1:nbn

%If the ith FE node X_i is a Dirichlet boundary node,then we reset the ith equation in the linear sysytem to be 1*u_i=u(X_i).
%Here u_i is the unknown associated with the node X_i and u(X_i)=Dirichelet_boundry_function_name(X_i).
%Also i=boundary_nodes(2,k) if boundary_nodes(1,k)==-1.

    if boundary_nodes(1,k)==-1 
        i=boundary_nodes(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i,1)=feval(Dirichlet_boundary_function_name,M(1,i),M(2,i));
    end

end