function [uh1,uh2]=curl_solver_triangle(left,right,bottom,top,h)
% Yu wei 2017-05
% We use P2--P1 element to solve this problem;
% N1,N2: x y 轴上网格划分数；
% N1_2.N2_2 x y 轴上自由度个数-1；
% tnp:P1元的节点个数；tnp_2:P2元的节点个数

global w

hf=h(1);
N1=(right-left)/h(1);N2=(top-bottom)/h(2);
N1_2=3*N1;N2_2=3*N2;

%L型区域上网格剖分

[M,T]=generate_tri_mesh_Lshape(left,right,bottom,top,h,3);
[~,T1]=generate_tri_mesh_Lshape(left,right,bottom,top,h,1);

tnp=(N1/2+1)*(N2+1)+(N2/2+1)*N2/2;
tnp_2=(N1_2/2+1)*(N2_2+1)+(N2_2/2+1)*N2_2/2;
b_edges=generate_boundary_nodes_edges_mesh_L(N1,N2,T(1:3,:));
number_of_elements=3*N1*N2/2;% 三角单位元个数
number_of_edges=length(b_edges);% 外部边界树
%u v 属于P2；其余属于P1；
uv=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,3,0,0,3,0,0);
uxvx=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,3,1,0,3,1,0);
uyvy=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,3,0,1,3,0,1);
uxvy=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,3,1,0,3,0,1);
uyvx=uxvy';
vw1=assemble_out_line_matrix(1,'a',M,T,T1,tnp,tnp_2,number_of_edges,b_edges,1,0,0,3,0,0);
vw2=assemble_out_line_matrix(2,'a',M,T,T1,tnp,tnp_2,number_of_edges,b_edges,1,0,0,3,0,0);
wxzx=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,1,0,1,1,0);
wyzy=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,0,1,1,0,1);
wz=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,0,0,1,0,0);
uz1=assemble_out_line_matrix(1,'a',M,T,T1,tnp,tnp_2,number_of_edges,b_edges,3,0,0,1,0,0);
uz2=assemble_out_line_matrix(2,'a',M,T,T1,tnp,tnp_2,number_of_edges,b_edges,3,0,0,1,0,0);
vp1=vw1;vp2=vw2;
uv1=assemble_out_line_matrix(3,'a',M,T,T1,tnp,tnp_2,number_of_edges,b_edges,3,0,0,3,0,0);
uv2=assemble_out_line_matrix(4,'a',M,T,T1,tnp,tnp_2,number_of_edges,b_edges,3,0,0,3,0,0);

A11=hf*uv1+uxvx+uyvy-w^2*uv;A12=uyvx-uxvy;A13=vw1;A14=-vp1;
A21=uxvy-uyvx;A22=uxvx+uyvy+hf*uv2-w^2*uv;A23=vw2;A24=-vp2;
A31=-uz1;A32=-uz2;A33=wxzx+wyzy+wz;A34=sparse(tnp,tnp);
A41=sparse(tnp,tnp_2);A42=A41;A43=A34;A44=A33;

A=[A11 A12 A13 A14;
    A21 A22 A23 A24;
    A31 A32 A33 A34
    A41 A42 A43 A44];

f1v1=assemble_vector('f1',M,T,T1,tnp,tnp_2,number_of_elements,3,0,0);
f2v2=assemble_vector('f2',M,T,T1,tnp,tnp_2,number_of_elements,3,0,0);

pv1n2=assemble_out_line_vector(2,'p',M,T,T1,b_edges,tnp,tnp_2,number_of_edges,3,0,0);
pv2n1=assemble_out_line_vector(1,'p',M,T,T1,b_edges,tnp,tnp_2,number_of_edges,3,0,0);
qv1n1=assemble_out_line_vector(1,'q',M,T,T1,b_edges,tnp,tnp_2,number_of_edges,3,0,0);
qv2n2=assemble_out_line_vector(2,'q',M,T,T1,b_edges,tnp,tnp_2,number_of_edges,3,0,0);

qq=assemble_out_vector('q',M,T,T1,b_edges,tnp,tnp_2,number_of_edges,1,0,0);

k1=f1v1-pv1n2+hf*qv1n1;k2=f2v2+pv2n1+hf*qv2n2;
k3=zeros(tnp,1);k4=qq;

k=[k1;k2;k3;k4];


r=A\k;

uh1=r(1:tnp_2);
uh2=r(tnp_2+1:2*tnp_2);
end

