function r=tri_reference_basis(x,y,basis_index,degree_x,degree_y,degree_P)
% Yu wei 2017-05

%This is the reference FE basis function on triangle ABC where A=(0,0), B=(1,0) and C=(0,1).
%x,y: the coordinates of the point where we want to evaluate the reference FE basis function.
%basis_type: the type of the FE.
%degree_P=1:2D linear FE.
%degree_P=2:2D Lagrange quadratic FE.
%basis_index: the index of FE basis function to specify which FE basis function we want to use.
%degree_x:the derivative degree of the FE basis function with respect to x.
%degree_y:the derivative degree of the FE basis function with respect to y.

%%
if degree_P==1
    r=tri_reference_basis_P1(x,y,basis_index,degree_x,degree_y);
elseif degree_P==2
    r=tri_reference_basis_P2(x,y,basis_index,degree_x,degree_y);
elseif degree_P==3    
    r=tri_reference_basis_P3(x,y,basis_index,degree_x,degree_y);

end

end
%%
function r=tri_reference_basis_P1(x,y,basis_index,degree_x,degree_y)

    if degree_x==0&&degree_y==0
        
        if basis_index==1
            r=1-x-y;
        elseif basis_index==2
            r=x;
        elseif basis_index==3
            r=y;
        end

    elseif degree_x==1&&degree_y==0
        
        if basis_index==1
            r=-1;
        elseif basis_index==2
            r=1;
        elseif basis_index==3
            r=0;
        end

    elseif degree_x==0&&degree_y==1
        
        if basis_index==1
            r=-1;
        elseif basis_index==2
            r=0;
        elseif basis_index==3
            r=1;
        end
        
    end
end

%%
function r=tri_reference_basis_P2(x,y,basis_index,degree_x,degree_y)
    if degree_x==0&&degree_y==0
        
        if basis_index==1
            r=1-3*x-3*y+2*x.^2+2*y.^2+4*x.*y;
        elseif basis_index==2
            r=2*x.^2-x;
        elseif basis_index==3
            r=2*y.^2-y;
        elseif basis_index==4
            r=4*x-4*x.^2-4*x.*y;
        elseif basis_index==5
            r=4*x.*y;
        elseif basis_index==6
            r=4*y-4*y.^2-4*x.*y;
        end
             
    elseif degree_x==1&&degree_y==0
 
        if basis_index==1
            r=-3+4*x+4*y;
        elseif basis_index==2
            r=4*x-1;
        elseif basis_index==3
            r=0;
        elseif basis_index==4
            r=4-8*x-4*y;
        elseif basis_index==5
            r=4*y;
        elseif basis_index==6
            r=-4*y;
        end           

                      
    elseif degree_x==0&&degree_y==1
            
        if basis_index==1
            r=-3+4*y+4*x;
        elseif basis_index==2
            r=0;
        elseif basis_index==3
            r=4*y-1;
        elseif basis_index==4
            r=-4*x;
        elseif basis_index==5
            r=4*x;
        elseif basis_index==6
            r=4-8*y-4*x;
        end
      
    elseif degree_x==2&&degree_y==0  
        
        if basis_index==1
            r=4;
        elseif basis_index==2
            r=4;
        elseif basis_index==3
            r=0;
        elseif basis_index==4
            r=-8;
        elseif basis_index==5
            r=0;
        elseif basis_index==6
            r=0;
        end

    elseif degree_x==0&&degree_y==2 

        if basis_index==1
            r=4;
        elseif basis_index==2
            r=0;
        elseif basis_index==3
            r=4;
        elseif basis_index==4
            r=0;
        elseif basis_index==5
            r=0;
        elseif basis_index==6
            r=-8;
        end

    elseif degree_x==1&&degree_y==1 

        if basis_index==1
            r=4;
        elseif basis_index==2
            r=0;
        elseif basis_index==3
            r=0;
        elseif basis_index==4
            r=-4;
        elseif basis_index==5
            r=4;
        elseif basis_index==6
            r=-4;
        end 
      
    end
    
end
function r=tri_reference_basis_P3(x,y,basis_index,degree_x,degree_y)
    if degree_x==0&&degree_y==0
        if basis_index==1
            r=-((3*x + 3*y - 1)*(3*x + 3*y - 2)*(x + y - 1))/2;
        elseif basis_index==2
            r=(x*(3*x - 1)*(3*x - 2))/2;
        elseif basis_index==3
            r=(y*(3*y - 1)*(3*y - 2))/2;
        elseif basis_index==4
            r=x*(3*x + 3*y - 2)*((9*x)/2 + (9*y)/2 - 9/2);
        elseif basis_index==5
            r=-(9*x*(3*x - 1)*(x + y - 1))/2;
        elseif basis_index==6
            r=(9*x*y*(3*x - 1))/2;
        elseif basis_index==7
            r=(9*x*y*(3*y - 1))/2;
        elseif basis_index==8
            r=-(9*y*(3*y - 1)*(x + y - 1))/2;
        elseif basis_index==9
            r=y*(3*x + 3*y - 2)*((9*x)/2 + (9*y)/2 - 9/2);
        elseif basis_index==10
            r=-x*y*(27*x + 27*y - 27);
        end

    elseif degree_x==1&&degree_y==0
        if basis_index==1
            r=- ((3*x + 3*y - 1)*(3*x + 3*y - 2))/2 - (3*(3*x + 3*y - 1)*(x + y - 1))/2 - (3*(3*x + 3*y - 2)*(x + y - 1))/2;
        elseif basis_index==2
            r=(3*x*(3*x - 1))/2 + (3*x*(3*x - 2))/2 + ((3*x - 1)*(3*x - 2))/2;
        elseif basis_index==3
            r=0;
        elseif basis_index==4
            r=(3*x + 3*y - 2)*((9*x)/2 + (9*y)/2 - 9/2) + (9*x*(3*x + 3*y - 2))/2 + 3*x*((9*x)/2 + (9*y)/2 - 9/2);
        elseif basis_index==5
            r=- (9*x*(3*x - 1))/2 - (27*x*(x + y - 1))/2 - (9*(3*x - 1)*(x + y - 1))/2;
        elseif basis_index==6
            r=(9*y*(3*x - 1))/2 + (27*x*y)/2;
        elseif basis_index==7
            r=(9*y*(3*y - 1))/2;
        elseif basis_index==8
            r=-(9*y*(3*y - 1))/2;
        elseif basis_index==9
            r=(9*y*(3*x + 3*y - 2))/2 + 3*y*((9*x)/2 + (9*y)/2 - 9/2);
        elseif basis_index==10
            r=- 27*x*y - y*(27*x + 27*y - 27);
        end

    elseif degree_x==0&&degree_y==1
        
        if basis_index==1
            r=- ((3*x + 3*y - 1)*(3*x + 3*y - 2))/2 - (3*(3*x + 3*y - 1)*(x + y - 1))/2 - (3*(3*x + 3*y - 2)*(x + y - 1))/2;
        elseif basis_index==2
            r=0;
        elseif basis_index==3
            r=(3*y*(3*y - 1))/2 + (3*y*(3*y - 2))/2 + ((3*y - 1)*(3*y - 2))/2;
        elseif basis_index==4
            r=(9*x*(3*x + 3*y - 2))/2 + 3*x*((9*x)/2 + (9*y)/2 - 9/2);
        elseif basis_index==5
            r=-(9*x*(3*x - 1))/2;

        elseif basis_index==6
            r=(9*x*(3*x - 1))/2;
        elseif basis_index==7
            r=(9*x*(3*y - 1))/2 + (27*x*y)/2;
        elseif basis_index==8
            r=- (9*y*(3*y - 1))/2 - (27*y*(x + y - 1))/2 - (9*(3*y - 1)*(x + y - 1))/2;
        elseif basis_index==9
            r=(3*x + 3*y - 2)*((9*x)/2 + (9*y)/2 - 9/2) + (9*y*(3*x + 3*y - 2))/2 + 3*y*((9*x)/2 + (9*y)/2 - 9/2);
        elseif basis_index==10
            r=- 27*x*y - x*(27*x + 27*y - 27);
        end
        
    end

end










































%  
% lamda1=tri_reference_basis_P1(x,y,1,0,0);
% lamda2=tri_reference_basis_P1(x,y,2,0,0);
% lamda3=tri_reference_basis_P1(x,y,3,0,0);
% lamda1dx=tri_reference_basis_P1(x,y,1,1,0);
% lamda2dx=tri_reference_basis_P1(x,y,2,1,0);
% lamda3dx=tri_reference_basis_P1(x,y,3,1,0);
% lamda1dy=tri_reference_basis_P1(x,y,1,0,1);
% lamda2dy=tri_reference_basis_P1(x,y,2,0,1);
% lamda3dy=tri_reference_basis_P1(x,y,3,0,1);
% %P2 element¡ª¡ª>dof=6,number of basis=6
% if degree_dx==0&&degree_dy==0
%    if basis_index==1
%        r=lamda1.*(2.*lamda1-1);
%    elseif basis_index==2
%        r=lamda2.*(2*lamda2-1);
%    elseif basis_index==3
%        r=lamda3.*(2*lamda3-1);
%    elseif basis_index==4
%        r=4*lamda1.*lamda2;
%    elseif basis_index==5   
%        r=4*lamda2.*lamda3;
%    elseif basis_index==6
%        r=4*lamda1.*lamda3;
%    else
%        fprintf('basis_index is more than the total number of basis functions')
%        return
%    end
% elseif degree_dx==1&&degree_dy==0
%    if basis_index==1
%        r=lamda1*4.*lamda1dx-lamda1dx;
%    elseif basis_index==2
%        r=lamda2*4.*lamda2dx-lamda2dx;
%    elseif basis_index==3
%        r=lamda3*4.*lamda3dx-lamda3dx;
%    elseif basis_index==4
%        r=4*lamda1dx.*lamda2+4*lamda1*lamda2dx;
%    elseif basis_index==5    
%        r=4*lamda2dx.*lamda3+4*lamda2*lamda3dx;
%    elseif basis_index==6
%        r=4*lamda1dx.*lamda3+4*lamda1*lamda3dx;
%    else
%        error('basis_index is more than the total number of basis functions')
%    end
% elseif degree_dx==0&&degree_dy==1
%    if basis_index==1
%        r=lamda1*4.*lamda1dy-lamda1dy;
%    elseif basis_index==2
%        r=lamda2*4.*lamda2dy-lamda2dy;
%    elseif basis_index==3
%        r=lamda3*4.*lamda3dy-lamda3dy;
%    elseif basis_index==4
%        r=4*lamda1dy.*lamda2+4*lamda1*lamda2dy;
%    elseif basis_index==5   
%        r=4*lamda2dy.*lamda3+4*lamda2*lamda3dy;
%    elseif basis_index==6
%        r=4*lamda1dy.*lamda3+4*lamda1*lamda3dy;
%    else
%        fprintf('basis_index is more than the total number of basis functions')
%        return
%    end
% elseif degree_dx==0&&degree_dy==2
%    if basis_index==1
%        r=4*lamda1dy.^2;
%    elseif basis_index==2
%        r=4*lamda2dy.^2;
%    elseif basis_index==3
%        r=4*lamda3dy.^2;
%    elseif basis_index==4
%        r=8*lamda1dy.*lamda2dy;
%    elseif basis_index==5   
%        r=8*lamda2dy.*lamda3dy;
%    elseif basis_index==6
%        r=8*lamda1dy.*lamda3dy;
%    else
%        fprintf('basis_index is more than the total number of basis functions')
%        return
%    end
% elseif degree_dx==2&&degree_dy==0
%    if basis_index==1
%        r=4*lamda1dx.^2;
%    elseif basis_index==2
%        r=4*lamda2dx.^2;
%    elseif basis_index==3
%        r=4*lamda3dx.^2;
%    elseif basis_index==4
%        r=8*lamda1dx.*lamda2dx;
%    elseif basis_index==5   
%        r=8*lamda2dx.*lamda3dx;
%    elseif basis_index==6
%        r=8*lamda1dx.*lamda3dx;
%    else
%        fprintf('basis_index is more than the total number of basis functions')
%        return
%    end
% elseif degree_dx==1&&degree_dy==1
%    if basis_index==1
%        r=lamda1dx*4.*lamda1dy;
%    elseif basis_index==2
%        r=lamda2dx*4.*lamda2dy;
%    elseif basis_index==3
%        r=lamda3dx*4.*lamda3dy;
%    elseif basis_index==4
%        r=4*lamda1dx.*lamda2dy+4*lamda1dy.*lamda2dx;
%    elseif basis_index==5   
%        r=4*lamda2dx.*lamda3dy+4*lamda2dy.*lamda3dx;
%    elseif basis_index==6
%        r=4*lamda1dx.*lamda3dy+4*lamda1dy.*lamda3dx;
%    else
%        fprintf('basis_index is more than the total number of basis functions')
%        return
%    end
% else
%     r=0;
% end
% 
