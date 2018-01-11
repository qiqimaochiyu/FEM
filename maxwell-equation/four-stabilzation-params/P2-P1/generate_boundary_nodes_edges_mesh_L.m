function boundary_edges=generate_boundary_nodes_edges_mesh_L(N1,N2,T)


%Information matrix for boundary nodes. It uses the index of FE, not the index of partition.

nbn=2*(N1+N2);
boundary_edges=zeros(4,nbn); 
boundary_edges(1,:)=-1;

for k=1:N1/2
    boundary_edges(2,k)=(k-1)*2*N2+1;
    boundary_edges(3,k)=T(1,boundary_edges(2,k));
    boundary_edges(4,k)=T(2,boundary_edges(2,k));
end

for k=(N1/2+1):(N1/2+N2/2)
    boundary_edges(2,k)=2*(k-N1/2-1)+(N1/2-1)*2*N2+1;
    boundary_edges(3,k)=T(2,boundary_edges(2,k));
    boundary_edges(4,k)=T(3,boundary_edges(2,k));
end

for k=(N1/2+N2/2+1):(N1+N2/2)
    boundary_edges(2,k)=N1*N2+1+(k-N2/2-N1/2-1)*N2;
    boundary_edges(3,k)=T(1,boundary_edges(2,k));
    boundary_edges(4,k)=T(2,boundary_edges(2,k));
end

for k=(N1+N2/2+1):(N1+N2)
    boundary_edges(2,k)=N1*N2+1+(N1/2-1)*N2+2*(k-N2/2-N1-1);
     boundary_edges(3,k)=T(2,boundary_edges(2,k));
    boundary_edges(4,k)=T(3,boundary_edges(2,k));
end

for k=(N1+N2+1):(3/2*N1+N2)
    boundary_edges(2,k)=3/2*N1*N2-N2*(k-N1-N2-1);
    boundary_edges(3,k)=T(2,boundary_edges(2,k));
    boundary_edges(4,k)=T(3,boundary_edges(2,k));
end

for k=(3/2*N1+N2+1):(2*N1+N2)
    boundary_edges(2,k)=N1*N2-2*N2*(k-3/2*N1-N2-1);
    boundary_edges(3,k)=T(2,boundary_edges(2,k));
    boundary_edges(4,k)=T(3,boundary_edges(2,k));
end

for k=(2*N1+N2+1):(2*N1+2*N2)
    boundary_edges(2,k)=2*N2-2*(k-2*N1-N2-1);
    boundary_edges(3,k)=T(3,boundary_edges(2,k));
    boundary_edges(4,k)=T(1,boundary_edges(2,k));
end
end