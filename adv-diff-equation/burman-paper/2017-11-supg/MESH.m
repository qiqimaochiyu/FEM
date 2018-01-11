%对网格剖分进行画图
global mesh_type
left=0;right=1;bottom=0;top=1;h_1=[1/8,1/8];
N1_basis=(right-left)/h_1(1);
N2_basis=(top-bottom)/h_1(2);
if mesh_type==0
    tnp=2*N1_basis*N2_basis;
    [M,T]=mesh_divise(left,right,bottom,top,h_1);
elseif mesh_type==1
    tnp=2*N1_basis*N2_basis;
    [M,T]=mesh_divise1(left,right,bottom,top,h_1);
elseif mesh_type==2
    tnp=2*N1_basis*N2_basis;
    [M,T]=mesh_divise2(left,right,bottom,top,h_1);
elseif mesh_type==3
    tnp=4*N1_basis*N2_basis;
    [M,T]=mesh_divise3(left,right,bottom,top,h_1); 
end
close all
for i=1:tnp
    plot([M(1,T(1,i)) M(1,T(2,i))], [M(2,T(1,i)) M(2,T(2,i))], 'k')
    hold on
    plot([M(1,T(1,i)) M(1,T(3,i))], [M(2,T(1,i)) M(2,T(3,i))], 'k')
    hold on
    plot([M(1,T(3,i)) M(1,T(2,i))], [M(2,T(3,i)) M(2,T(2,i))], 'k')
    hold on

    set(gca,'xtick',[]);    set(gca,'ytick',[]); %去掉途中坐标轴
end

