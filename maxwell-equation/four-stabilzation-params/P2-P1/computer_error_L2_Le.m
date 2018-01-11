function [L2,Le]=computer_error_L2_Le(uh1,uh2,M,T)
%Yu wei 2017-05

 %% computer the error
   data=maxwell_example2;
   u_abs1=zeros(length(uh1),1);
   u_abs2=zeros(length(uh2),1);
   L2_error_u1=FE_solution_error_triangle(uh1,data.exact_u1,M,T,0,0);

   L2_error_u2=FE_solution_error_triangle(uh2,data.exact_u2,M,T,0,0);
   L2_error=sqrt(L2_error_u1+L2_error_u2);
   L2_abs_u1=FE_solution_error_triangle(u_abs1,data.exact_u1,M,T,0,0);
   L2_abs_u2=FE_solution_error_triangle(u_abs2,data.exact_u2,M,T,0,0);
   L2_abs=sqrt(L2_abs_u1+L2_abs_u2);
   L2=L2_error/L2_abs;
%curl·¶Êý   
   u1y=FE_solution_error_triangle(uh1,data.exact_u1y,M,T,0,1);
   u2x=FE_solution_error_triangle(uh2,data.exact_u2x,M,T,1,0);
   u1y_abs=FE_solution_error_triangle(u_abs1,data.exact_u1y,M,T,0,1);
   u2x_abs=FE_solution_error_triangle(u_abs2,data.exact_u2x,M,T,1,0);
   E=sqrt(u2x)-sqrt(u1y);
   E_abs=sqrt(u2x_abs)-sqrt(u1y_abs);
   if E_abs==0
        Le=E;
   else
        Le=E/E_abs;
   end
   
 