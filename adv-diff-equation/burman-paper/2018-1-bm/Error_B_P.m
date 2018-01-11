
format short e

global b h_1

loop=4;loop_b=1;
L2=zeros(loop_b,loop);
H1=L2;L_max=L2;
% 
% disp('----------BM-example1----------')
% for j=1:loop_b
%     for t=1:4
%     for i=1:loop
%     b=1;data=example1;
%     disp('begin loop: ');disp(i);
%     left=0;right=1;bottom=0;top=1;
%     a_1=1;a_2=0;p=0; q=0.0125;e=10^(-5); 
%     Gauss_point_number=9;
%     mesh_type=t;
%     
%     h_1=[1/16 1/16]*(1/2)^i;
%     uh=B_P_solver_triangle(data,mesh_type,q,left,right,bottom,top,h_1,e,b,p,a_1,a_2);
% 
% 
%     [L2(j,i),H1(j,i),L_max(j,i)]=computer_error_L2_H1_max(mesh_type,data,uh);
% 
%     end
%     figure;MESH;
%     disp('b='); disp(b)
%     disp('mesh_type='); disp(mesh_type)
%     result_output(L2(j,:),H1(j,:),L_max(j,:))
%     end
% end

disp('----------BM-example2----------')


for j=1:loop_b
    for t=1:4
    for i=1:loop
    b=1;data=example2;
    disp('begin loop: ');disp(i);
    left=0;right=1;bottom=0;top=1;
    a_1=1;a_2=0;p=0; q=0.0125;e=10^(-5); 
    Gauss_point_number=9;
    mesh_type=t;
   %  figure;MESH;

    h_1=[1/16,1/16]*(1/2)^i;
    uh=B_P_solver_triangle(data,mesh_type,q,left,right,bottom,top,h_1,e,b,p,a_1,a_2);
    [L2(j,i),H1(j,i),L_max(j,i)]=computer_error_L2_H1_max(mesh_type,data,uh);

    end
    disp('b='); disp(b)
    disp('mesh_type='); disp(mesh_type)
    result_output(L2(j,:),H1(j,:),L_max(j,:))
    end

end

disp('-----------example2----------')

for j=1:loop_b
    for t=1:4
    for i=1:loop
    b=1;data=example2;
    disp('begin loop: ');disp(i);
    left=0;right=1;bottom=0;top=1;
    a_1=1;a_2=0;
    p=0.54; q=0.0125;e=10^(-5); 
    Gauss_point_number=9;
    mesh_type=t;
    %  figure;MESH;

    h_1=[1/16,1/16]*(1/2)^i;
    uh=B_P_solver_triangle(data,mesh_type,q,left,right,bottom,top,h_1,e,b,p,a_1,a_2);
    [L2(j,i),H1(j,i),L_max(j,i)]=computer_error_L2_H1_max(mesh_type,data,uh);

    end
    disp('b='); disp(b)
    disp('mesh_type='); disp(mesh_type)
    result_output(L2(j,:),H1(j,:),L_max(j,:))
    end
end

% Error_B_P2;
% plot_f;




