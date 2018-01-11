function [boundary_nodes,boundary_edges]=generate_boundary_nodes_edges(N1_basis,N2_basis,N1_partition,N2_partition)




%Information matrix for boundary nodes. It uses the index of FE, not the index of partition.

nbn=2*(N1_basis+N2_basis);
boundary_nodes=zeros(2,nbn); 

%The following boundary condition may change for different problems.
%All Dirichlet boundary nodes.
boundary_nodes(1,:)=-1;

%The index in the following is associated with the index in "generate_M_T_triangle.m".
%bottom boundary nodes.
for k=1:N1_basis
    boundary_nodes(2,k)=(k-1)*(N2_basis+1)+1;
end

%right boundary nodes.
for k=N1_basis+1:N1_basis+N2_basis
    boundary_nodes(2,k)=N1_basis*(N2_basis+1)+k-N1_basis;
end

%top boundary nodes.
for k=N1_basis+N2_basis+1:2*N1_basis+N2_basis
    boundary_nodes(2,k)=(2*N1_basis+N2_basis+2-k)*(N2_basis+1);
end

%left boundary nodes.
for k=2*N1_basis+N2_basis+1:nbn
    boundary_nodes(2,k)=2*N1_basis+2*N2_basis+2-k;
end





%Information matrix for boundary edges. It uses the index of partition, not the index of FE.

nbe=2*(N1_partition+N2_partition);
boundary_edges=zeros(4,nbe);


%The following boundary condition may change for different problems.
%All Dirichlet boundary edges.
boundary_edges(1,:)=-1;

%The index in the following is associated with the index in "generate_M_T_triangle.m".
%bottom boundary edges.
for k=1:N1_partition
    boundary_edges(2,k)=(k-1)*2*N2_partition+1;
    boundary_edges(3,k)=(k-1)*(N2_partition+1)+1;
    boundary_edges(4,k)=k*(N2_partition+1)+1;
end

%right boundary edges.
for k=N1_partition+1:N1_partition+N2_partition
    boundary_edges(2,k)=(N1_partition-1)*2*N2_partition+2*(k-N1_partition)-1;
    boundary_edges(3,k)=N1_partition*(N2_partition+1)+k-N1_partition;
    boundary_edges(4,k)=N1_partition*(N2_partition+1)+k-N1_partition+1;
end

%top boundary edges.
for k=N1_partition+N2_partition+1:2*N1_partition+N2_partition
    boundary_edges(2,k)=(2*N1_partition+N2_partition+1-k)*2*N2_partition;
    boundary_edges(3,k)=(2*N1_partition+N2_partition+2-k)*(N2_partition+1);
    boundary_edges(4,k)=(2*N1_partition+N2_partition+1-k)*(N2_partition+1);
end

%left boundary edges.
for k=2*N1_partition+N2_partition+1:nbe
    boundary_edges(2,k)=2*(2*N1_partition+2*N2_partition+1-k);
    boundary_edges(3,k)=2*N1_partition+2*N2_partition+2-k;
    boundary_edges(4,k)=2*N1_partition+2*N2_partition+1-k;
end



