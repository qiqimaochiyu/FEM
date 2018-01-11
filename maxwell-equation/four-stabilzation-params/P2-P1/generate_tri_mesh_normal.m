function [M,T]=generate_tri_mesh_normal(left,right,bottom,top,h,degree_P)
% Yu Wei 2017-05
% Generate the matrix M and T for a tiangular mesh of a rectangular domain
% [left,right]*[bottom,top];
% h_f:the step size of the partition.
% We first form a rectangular mesh with element size h_f(1)*h_f(2).
% The tri_mesh is formed by cutting each rec into tri along the diagonal.
% degree_P:the type of the FE.
% degree_P=1:1D linear FE; 2: Lagrange quadratic FE.
% M(i,j): ith coordinate of the jth node.
% T(i,j): store the global index of the ith node in the jth element.
% 2 - - - 4
% | 2  /  |
% |  /  1 |  T(:,1)=[1 3 4] T(:,2)=[1 4 2]
% 1 - - - 3
% N1: the number of the sub-interval of the partition in x-direction.
% N2: the number of the sub-interval of the partition in y-direction.
% tnp: total number of the nodes.

%% P1

if degree_P==1
    
    N1=(right-left)/h(1);
    N2=(top-bottom)/h(2);
    tnp=(N1+1)*(N2+1);
    M=zeros(2,tnp);
    T=zeros(3,2*N1*N2);
    Q=zeros(N1+1,N2+1); 
    
    for j=1:tnp
        if mod(j,N2+1)==0
            M(1,j)=left+(j/(N2+1)-1)*h(1);
            M(2,j)=top;
        else
            M(1,j)=left+fix(j/(N2+1))*h(1);
            M(2,j)=bottom+(mod(j,N2+1)-1)*h(2);
        end
    end
    
    for i=1:N1+1
        for j=1:N2+1
            Q(i,j)=(i-1)*(N2+1)+j;
        end
    end
    
% go through all rec in the partition.
% for the nth rec, store the information of its two tri elements.

    for n=1:N1*N2
        
        if mod(n,N2)==0
            row=N2;
            column=n/N2;
        else
            row=mod(n,N2);
            column=fix(n/N2)+1;
        end
        
        T(1,2*n-1)=Q(column,row);
        T(2,2*n-1)=Q(column+1,row);
        T(3,2*n-1)=Q(column+1,row+1);
        
        T(1,2*n)=Q(column,row);
        T(2,2*n)=Q(column+1,row+1);
        T(3,2*n)=Q(column,row+1);
    end
 
elseif degree_P==2
%% P2
    
   N1=(right-left)/h(1);
   N2=(top-bottom)/h(2);
   dh=h/2;
   dN1=N1*2;
   dN2=N2*2;
   tnp=(dN1+1)*(dN2+1);
   M=zeros(2,tnp);
   T=zeros(6,2*N1*N2);
   Q=zeros(dN1+1,dN2+1);

   for j=1:tnp
      if mod(j,dN2+1)==0
         M(1,j)=left+(j/(dN2+1)-1)*dh(1);
         M(2,j)=top;
      else
         M(1,j)=left+fix(j/(dN2+1))*dh(1);
         M(2,j)=bottom+(mod(j,dN2+1)-1)*dh(2);
      end
   end

   for i=1:dN1+1
      for j=1:dN2+1
         Q(i,j)=(i-1)*(dN2+1)+j;
      end
   end

%Go through all rec in the partition. 
%For the nth rec, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,2*n-1)=Q(2*column-1,2*row-1);
      T(2,2*n-1)=Q(2*column+1,2*row-1); 
      T(3,2*n-1)=Q(2*column+1,2*row+1);
      T(4,2*n-1)=Q(2*column,2*row-1);
      T(5,2*n-1)=Q(2*column+1,2*row);
      T(6,2*n-1)=Q(2*column,2*row);


      T(1,2*n)=Q(2*column-1,2*row-1);
      T(2,2*n)=Q(2*column+1,2*row+1);
      T(3,2*n)=Q(2*column-1,2*row+1);
      T(4,2*n)=Q(2*column,2*row);
      T(5,2*n)=Q(2*column,2*row+1);
      T(6,2*n)=Q(2*column-1,2*row); 
      
  
   end
end

end

