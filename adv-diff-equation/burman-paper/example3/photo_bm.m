format short e

global b h_1
% disp('----------BM-example2-QIU----------')

b=0;data=example_q;
left=0;right=1;bottom=0;top=1;
a_1=1/2;a_2=sqrt(3/4);e=10^(-6);
p=0.53;q=0.0125;
Gauss_point_number=9;
mesh_type=1;
    
h_1=[1/16 1/16]*(1/2);
uh=B_P_solver_triangle(data,mesh_type,q,left,right,bottom,top,h_1,e,b,p,a_1,a_2);
 
result_mesh;  