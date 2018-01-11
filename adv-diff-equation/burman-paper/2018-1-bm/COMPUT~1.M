function [L2,H1,L_max]=computer_error_L2_H1_max(mesh_type,data,uh)

global h_1
    left=0;right=1;bottom=0;top=1;
    a_1=1;a_2=0;e=10^(-5); 
    Gauss_point_number=9;
    L2=FE_solution_error_triangle(mesh_type,uh,data.u,left,right,bottom,top,...
        h_1,0,0,Gauss_point_number);
    H1_error_x=FE_solution_error_triangle(mesh_type,uh,data.u_x,left,right,bottom,top,...
        h_1,1,0,Gauss_point_number);
    H1_error_y=FE_solution_error_triangle(mesh_type,uh,data.u_y,left,right,bottom,top,...
        h_1,0,1,Gauss_point_number);
    H1=sqrt(H1_error_x^2+H1_error_y^2+L2^2);
    L_max=FE_solution_error_infinity_norm_triangle(mesh_type,uh,data.u,left,right,...
        bottom,top,h_1,0,0,Gauss_point_number);
    end
