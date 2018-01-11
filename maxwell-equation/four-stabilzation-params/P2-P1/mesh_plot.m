function result_plot

left=-1;right=1;bottom=-1;top=1;hf=[1/4 1/4];degree_P=1;
[P,T]=generate_tri_mesh_normal(left,right,bottom,top,hf,degree_P);
for i=1:size(T,2)
    A=P(:,T(:,i));
    xx=A(1,:);
    yy=A(2,:);
    plot(xx,yy,'k');
    hold on;
end
end



