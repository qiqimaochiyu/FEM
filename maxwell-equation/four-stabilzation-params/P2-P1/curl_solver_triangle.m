function [uh1,uh2]=curl_solver_triangle(left,right,bottom,top,h)
% Yu wei 2017-05
% We use P2--P1 element to solve this problem;
% N1,N2: x y 轴上网格划分数；
% N1_2.N2_2 x y 轴上自由度个数-1；
% tnp:P1元的节点个数；tnp_2:P2元的节点个数

global w domain_type degree_P data


hf=h(1);
N1=(right-left)/h(1);N2=(top-bottom)/h(2);

if domain_type== 1%方形区域
    tnp=(N1+1)*(N2+1);
    [M,T]=generate_tri_mesh_normal(left,right,bottom,top,h,degree_P);
    [b_nodes,b_edges]=generate_boundary_nodes_edges(N1,N2,N1,N2);
    number_of_elements=2*N1*N2;
elseif domain_type== 2%L型区域
    tnp=(N1/2+1)*(N2+1)+(N2/2+1)*N2/2;
    [M,T]=generate_tri_mesh_Lshape(left,right,bottom,top,h,degree_P);
    b_edges=generate_boundary_nodes_edges_mesh_L(N1,N2,T);
    number_of_elements=3*N1*N2/2;
    
end
number_of_edges=length(b_edges);% 外部边界树
[Boundary,T_plus,T_index]=generate_in_edges(T,tnp);
T1=T;tnp_2=tnp;
uv=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,0,0,1,0,0);
uxvx=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,1,0,1,1,0);
uyvy=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,0,1,1,0,1);
uxvy=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,1,0,1,0,1);
uyvx=uxvy';
uxv=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,1,0,1,0,0);
uyv=assemble_matrix(M,T,T1,number_of_elements,tnp,tnp_2,1,0,1,1,0,0);
uvx=uxv';uvy=uyv';

uxvxF=assemble_in_line_matrix_plus('a',M,T,T_index,T_plus,tnp,tnp_2,number_of_edges,Boundary,1,1,0,1,1,0);
uyvyF=assemble_in_line_matrix_plus('a',M,T,T_index,T_plus,tnp,tnp_2,number_of_edges,Boundary,1,0,1,1,0,1);
uxvyF=assemble_in_line_matrix_plus('a',M,T,T_index,T_plus,tnp,tnp_2,number_of_edges,Boundary,1,1,0,1,0,1);
uyvxF=assemble_in_line_matrix_plus('a',M,T,T_index,T_plus,tnp,tnp_2,number_of_edges,Boundary,1,0,1,1,1,0);

uvt1=assemble_out_line_matrix(3,'a',M,T,T1,tnp,tnp_2,number_of_edges,b_edges,1,0,0,1,0,0);
uvt2=assemble_out_line_matrix(4,'a',M,T,T1,tnp,tnp_2,number_of_edges,b_edges,1,0,0,1,0,0);

a11=hf*uvt1+hf^2*w^2*uv+hf^2*uxvx+uyvyF;a12=-uxvyF+hf^2*uyvx;a13=uyvy-w*uv;
a14=-uxvy;a15=-uyvy+w*uv;a16=uxvy;a17=-uxv;a18=-uvy-uyv;a19=uxv;a110=uvy+uyv;
 
a21=-uyvxF+hf^2*uxvy;a22=hf*uvt2+uxvxF+hf^2*uyvy+hf^2*w^2*uv;a23=-uyvx;
a24=uxvx-w*uv;a25=uyvx;a26=-uxvx+w*uv;a27=-uyv;a28=uvx+uxv;a29=uyv;a210=-uvx-uxv;
a31=-uyvy+w*uv;a32=uxvy;a33=uv+uxvx+uyvy;
a41=uyvx;a42=-uxvx+w*uv;a44=uv+uyvy+uxvx;
a55=uv+uxvx+uyvy;a66=uv+uyvy+uxvx;
a71=uvx;a72=uvy;a77=uv+uxvx+uyvy;
a81=uyv+uvy;a82=-uxv-uvx;a88=uv+uxvx+uyvy;
a99=a88;a1010=a99;a0=sparse(tnp,tnp);
A=[a11 a12 a13 a14 a15 a16 a17 a18 a19 a110;
   a21 a22 a23 a24 a25 a26 a27 a28 a29 a210;
   a31 a32 a33 a0 a0 a0 a0 a0 a0 a0;
   a41 a42 a0 a44 a0 a0 a0 a0 a0 a0;
   a0 a0 a0 a0 a55 a0 a0 a0 a0 a0;
   a0 a0 a0 a0 a0 a66 a0 a0 a0 a0;
   a71 a72 a0 a0 a0 a0 a77 a0 a0 a0;
   a81 a82 a0 a0 a0 a0 a0 a88 a0 a0;
   a0 a0 a0 a0 a0 a0 a0 a0 a99 a0;
   a0 a0 a0 a0 a0 a0 a0 a0 a0 a1010;
];


f1v=assemble_vector(data.f1,M,T,T1,tnp,tnp_2,number_of_elements,1,0,0);
f2v=assemble_vector(data.f2,M,T,T1,tnp,tnp_2,number_of_elements,1,0,0);
gvx=assemble_vector(data.g,M,T,T1,tnp,tnp_2,number_of_elements,1,1,0);
gvy=assemble_vector(data.g,M,T,T1,tnp,tnp_2,number_of_elements,1,0,1);
Xvt1=assemble_out_line_vector(1,data.fai,M,T,T1,b_edges,tnp,tnp_2,number_of_edges,1,0,0);
Xvt2=assemble_out_line_vector(2,data.fai,M,T,T1,b_edges,tnp,tnp_2,number_of_edges,1,0,0);
gv=assemble_vector(data.g,M,T,T1,tnp,tnp_2,number_of_elements,1,0,0);
Xv=assemble_out_vector(data.fai,M,T,T1,b_edges,tnp,tnp_2,number_of_edges,1,0,0);


k1=-hf^2*w*f1v+hf^2*gvx+hf*Xvt1;
k2=-hf^2*w*f2v+hf^2*gvy+hf*Xvt2;
k3=sparse(tnp,1);k4=k3;
k5=f1v;k6=f2v;k7=k4;k8=k4;
k9=gv;k10=Xv;
k=[k1;k2;k3;k4;k5;k6;k7;k8;k9;k10];

[A,k]=treat_Dirichlet_boundary_triangle(tnp,data.g,A,k,b_nodes,M);

r=A\k;
uh1=r(1:tnp);
uh2=r(tnp+1:2*tnp);
end



