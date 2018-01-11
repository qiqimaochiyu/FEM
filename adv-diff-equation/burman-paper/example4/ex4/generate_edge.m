function [boundary_nodes,boundary_edges]=generate_edge(v,T)

temp=length(v);
boundary_nodes=zeros(2,temp);
boundary_edges=zeros(4,temp);

boundary_nodes(1,:)=v(5,:);
boundary_nodes(2,:)=v(1,:);
boundary_edges(1,:)=v(5,:);
boundary_edges(3,:)=v(1,:);
boundary_edges(4,:)=v(2,:);

for i=1:temp
    [~,a]=find(T==v(1,i));
    [~,b]=find(T==v(2,i));
    boundary_edges(2,i)=intersect(a,b);
end



