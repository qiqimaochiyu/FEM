function [L2,Le]=computer_error_L2_Le(uh1,uh2,M,T)
%Yu wei 2017-05

 %% computer the error
   u_abs1=zeros(length(uh1),1);
   u_abs2=zeros(length(uh2),1);
   L2_error_u1=FE_solution_error_triangle(uh1,'u1',M,T,0,0);

   L2_error_u2=FE_solution_error_triangle(uh2,'u2',M,T,0,0);
   L2_error=sqrt(L2_error_u1+L2_error_u2);
   L2_abs_u1=FE_solution_error_triangle(u_abs1,'u1',M,T,0,0);
   L2_abs_u2=FE_solution_error_triangle(u_abs2,'u2',M,T,0,0);
   L2_abs=sqrt(L2_abs_u1+L2_abs_u2);
   L2=L2_error/L2_abs;
%ÄÜÁ¿·¶Êý   
   u1u1=FE_solution_error_triangle_2(uh1,'u1',M,T,0,0,uh1,'u1',0,0);
   u2u2=FE_solution_error_triangle_2(uh2,'u2',M,T,0,0,uh2,'u2',0,0);
   u1xu1x=FE_solution_error_triangle_2(uh1,'u1x',M,T,1,0,uh1,'u1x',1,0);
   u1xu2y=FE_solution_error_triangle_2(uh1,'u1x',M,T,1,0,uh2,'u2y',0,1);
   u1yu2x=FE_solution_error_triangle_2(uh1,'u1y',M,T,0,1,uh2,'u2x',1,0);
   u1yu1y=FE_solution_error_triangle_2(uh1,'u1y',M,T,0,1,uh1,'u1y',0,1);
   u2yu1x=FE_solution_error_triangle_2(uh2,'u2y',M,T,0,1,uh1,'u1x',1,0);
   u2yu2y=FE_solution_error_triangle_2(uh2,'u2y',M,T,0,1,uh2,'u2y',0,1);
   u2xu1y=FE_solution_error_triangle_2(uh2,'u2x',M,T,1,0,uh1,'u1y',0,1);
   u2xu2x=FE_solution_error_triangle_2(uh2,'u2x',M,T,1,0,uh2,'u2x',1,0);

   u1u1_abs=FE_solution_error_triangle_2(u_abs1,'u1',M,T,0,0,u_abs1,'u1',0,0);
   u2u2_abs=FE_solution_error_triangle_2(u_abs2,'u2',M,T,0,0,u_abs2,'u2',0,0);
   u1xu1x_abs=FE_solution_error_triangle_2(u_abs1,'u1x',M,T,1,0,u_abs1,'u1x',1,0);
   u1xu2y_abs=FE_solution_error_triangle_2(u_abs1,'u1x',M,T,1,0,u_abs2,'u2y',0,1);
   u1yu2x_abs=FE_solution_error_triangle_2(u_abs1,'u1y',M,T,0,1,u_abs2,'u2x',1,0);
   u1yu1y_abs=FE_solution_error_triangle_2(u_abs1,'u1y',M,T,0,1,u_abs1,'u1y',0,1);
   u2yu1x_abs=FE_solution_error_triangle_2(u_abs2,'u2y',M,T,0,1,u_abs1,'u1x',1,0);
   u2yu2y_abs=FE_solution_error_triangle_2(u_abs2,'u2y',M,T,0,1,u_abs2,'u2y',0,1);
   u2xu1y_abs=FE_solution_error_triangle_2(u_abs2,'u2x',M,T,1,0,u_abs1,'u1y',0,1);
   u2xu2x_abs=FE_solution_error_triangle_2(u_abs2,'u2x',M,T,1,0,u_abs2,'u2x',1,0);

   E=sqrt(u2xu2x-u1yu2x-u2xu1y+u1yu1y+u1xu1x+u1xu2y+u2yu2y+u2yu1x+L2_error_u1+L2_error_u2);
   E_abs=sqrt(u2xu2x_abs-u1yu2x_abs-u2xu1y_abs+u1yu1y_abs+...
   u1xu1x_abs+u1xu2y_abs+u2yu2y_abs+u2yu1x_abs+L2_abs_u1+L2_abs_u2);
   if E_abs==0
        Le=E;
   else
        Le=E/E_abs;
   end
   
 