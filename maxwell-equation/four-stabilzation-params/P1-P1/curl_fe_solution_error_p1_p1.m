format short e
%P1-P1 Yu wei 2017-05
global Gauss_coefficient_reference_1D Gauss_point_reference_1D... 
Gauss_coefficient_reference_triangle Gauss_point_reference_triangle w...
domain_type degree_P data

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(9);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(8);
degree_P=1;

data=maxwell_example1;
num=6;
L2=zeros(1,num);
Le=zeros(1,num);
% domain_type='square';fprintf 'the domian is square \n\n'
domain_type=1;fprintf 'the domian is square \n\n'
for j=4
for i=1:num
  fprintf 'start loop:' 
  disp(i)
  W=[1 15 30 45];
  w=W(j);
  h=[1/2 1/2]*(1/2)^i;
  fprintf 'the size of the mesh:' 
  disp(1/h(1))
  left=0;right=1;
  bottom=0;top=1;
  
  [M,T]=generate_tri_mesh_normal(left,right,bottom,top,h,1);
  [uh1,uh2]=curl_solver_triangle(left,right,bottom,top,h);
  [L2(i),Le(i)]=computer_error_L2_Le(uh1,uh2,M,T);

end
fprintf 'w='
disp(w)
result_output(L2,Le);
fe_u_plot_nonsmooth(1,uh1,uh2,h)
fprintf '-------------Finished------------\n'
end
