function [M,T]=mesh_divise_4(left,right,bottom,top,h_1)



N1=(right-left)/h_1(1);
   N2=(top-bottom)/h_1(2);
   tnp=(N1+1)*(N2+1);
   tmp=N1*N2;
   M=zeros(2,tnp+tmp);
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
  for j=1:tmp
      if mod(j,N2)==0
          M(1,j+tnp)=left+1/2*h_1(1)+(j/N2-1)*h_1(1);
          M(2,j+tnp)=top-1/2*h_1(1);
      else
          M(1,j+tnp)=left+fix(j/N2)*h_1(1)+1/2*h_1(1);
          M(2,j+tnp)=bottom+(mod(j,N2)-1)*h_1(2)+1/2*h_1(1);
      end
  end
  
          
   
  for i=1:N1+1
      for j=1:N2+1
         Q(i,j)=(i-1)*(N2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:tmp 
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,4*n-3)=Q(column,row);
      T(2,4*n-3)=Q(column+1,row);
      T(3,4*n-3)=n+tnp;  
      
      T(1,4*n-2)=Q(column,row);
      T(2,4*n-2)=n+tnp;
      T(3,4*n-2)=Q(column,row+1);  
      
      T(1,4*n-1)=Q(column,row+1);
      T(2,4*n-1)=n+tnp;
      T(3,4*n-1)=Q(column+1,row+1);  
      
      T(1,4*n)=Q(column+1,row+1);
      T(2,4*n)=n+tnp;
      T(3,4*n)=Q(column+1,row);  
      
      
    
   end