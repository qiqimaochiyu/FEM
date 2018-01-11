% h_1=[1/32,1/32];
% left=0;
% right=1;
% bottom=0;
% top=1;
% e=0.000001;
% b=32;
% a_1=sqrt(1/2);
% a_2=sqrt(1/2);
% 
% [uh,A,K1,K2,K3]=poisson_solver_triangle(left,right,bottom,top,h_1,e,b,a_1,a_2);
% uh=uh*b;
% 

 load test2
nt=192;
a=zeros(1089,nt);
a(:,1)=uh;
for i=2:nt
a(:,i)=poisson_solver_triangle1(a(:,i-1),A,K1,K2,K3,left,right,bottom,top,h_1,e,b,a_1,a_2)*b;

end

for n=1:10
set(gcf, 'doublebuffer', 'on')
for i = 1 : nt
  Z=reshape(a(:,i),33,33);
x=linspace(0,1,33);
y=linspace(0,1,33);
[X,Y]=meshgrid(x,y);
colormap(jet)
mesh(Y,X,Z)
view([135,45])

    drawnow
end
end