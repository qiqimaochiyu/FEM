function [M,T]=mesh_divise_3(left,right,bottom,top,h_1)



N1=(right-left)/h_1(1);
   N2=(top-bottom)/h_1(2);
   tnp=(N1+1)*(N2+1);
   M=zeros(2,tnp);
   T=zeros(3,2*N1*N2);
   Q=zeros(N1+1,N2+1);

   for j=1:tnp
      if mod(j,N2+1)==0
         M(1,j)=left+(j/(N2+1)-1)*h_1(1);
         M(2,j)=top;
      else
         M(1,j)=left+fix(j/(N2+1))*h_1(1);
         M(2,j)=bottom+(mod(j,N2+1)-1)*h_1(2);
      end
   end

   for i=1:N1+1
      for j=1:N2+1
         Q(i,j)=(i-1)*(N2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
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
      if randsrc(1,1,randperm(2))==1
      % if mod(n-1,N2)~=0 && randsrc(1,1,randperm(2))==1

         T(3,2*n-1)=Q(column,row+1);  
         T(4,2*n-1)=-1;
         T(1,2*n)=Q(column,row+1);
         T(2,2*n)=Q(column+1,row);
         T(3,2*n)=Q(column+1,row+1);  
         T(4,2*n)=-1;
      
      else 
              
        T(3,2*n-1)=Q(column+1,row+1);  
        T(4,2*n-1)=-2;
        T(1,2*n)=Q(column,row);
        T(2,2*n)=Q(column+1,row+1);
        T(3,2*n)=Q(column,row+1); 
        T(4,2*n)=-2;
      end
  %T(4,n)=-1 表示剖分区域对角线斜率为－1； －2表示斜率为1；即向右上角。    
   end