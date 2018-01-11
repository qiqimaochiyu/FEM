format short e

global b h_1

loop=2;loop_b=1;
L2=zeros(loop_b,loop);
H1=L2;L_max=L2;
disp('----------SUPG-example1----------')
for j=1:loop_b
    for t=1
    for i=1:loop
    b=1;data=example1;
    disp('begin loop: ');disp(i);
    left=0;right=1;bottom=0;top=1;
    a_1=1;a_2=0;e=10^(-5); 
    Gauss_point_number=9;
    mesh_type=t;
    
    h_1=[1/16 1/16]*(1/2)^i;
    if t<4
        uh=SUPG_solver_triangle(data,mesh_type,left,right,bottom,top,h_1,e,a_1,a_2);
    else 
        uh=SUPG_solver_triangle2(data,mesh_type,left,right,bottom,top,h_1,e,a_1,a_2);
    end
    [L2(j,i),H1(j,i),L_max(j,i)]=computer_error_L2_H1_max(mesh_type,data,uh);

    end
    % figure;MESH;
    disp('b='); disp(b)
    disp('mesh_type='); disp(mesh_type)
    result_output(L2(j,:),H1(j,:),L_max(j,:))
    end
end

disp('----------SUPG-example2----------')
for j=1:loop_b
    for t=1:4
    for i=1:loop
    b=1;data=example2;
    disp('begin loop: ');disp(i);
    left=0;right=1;bottom=0;top=1;
    a_1=1;a_2=0;e=10^(-5); 
    Gauss_point_number=9;
    mesh_type=t;
    
    h_1=[1/16 1/16]*(1/2)^i;
    if t<4
        uh=SUPG_solver_triangle(data,mesh_type,left,right,bottom,top,h_1,e,a_1,a_2);
    else 
        uh=SUPG_solver_triangle2(data,left,right,bottom,top,h_1,e,a_1,a_2);
    end
    [L2(j,i),H1(j,i),L_max(j,i)]=computer_error_L2_H1_max(mesh_type,data,uh);

    end
    % figure;MESH;
    disp('b='); disp(b)
    disp('mesh_type='); disp(mesh_type)
    result_output(L2(j,:),H1(j,:),L_max(j,:))
    end
end