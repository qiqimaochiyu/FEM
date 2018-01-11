function [M,T]=generate_tri_mesh_Lshape(left,right,bottom,top,h,degree_P)
% Yu wei 2017-05
% domain:[-1 1]*[-1 1]/[0 1]*[-1 0] 
% divide the domain into two parts [-1 0]*[-1 1] and [0 1]*[0 1]
% use the fuction generate_tri_mesh_mormal

%% first part
[M1,T1]=generate_tri_mesh_normal(left,0,bottom,top,h,degree_P);
m=size(M1);
%% second part
[M2,T2]=generate_tri_mesh_normal(0,right,0,top,h,degree_P);
n=size(T2);
T2=T2+ones(n(1),n(2))*(m(2)-1-top/h(2)*degree_P);%加上左侧编号
M2(:,1:1+top/h(2)*degree_P)=[];%去掉重合的部分
M=[M1 M2];T=[T1 T2];
%% plot
% for i=1:size(T,2)
%     A=M(:,T(:,i));
%     xx=A(1,:);
%     yy=A(2,:);
%     plot(xx,yy,'k');
%     hold on;
%     set(gca,'XTick',-1:1:1)
%     set(gca,'YTick',-1:1:1)
% end
end


