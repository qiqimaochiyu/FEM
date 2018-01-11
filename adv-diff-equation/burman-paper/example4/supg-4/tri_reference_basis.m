function r=tri_reference_basis(x,y,basis_index,derivative_degree_x,derivative_degree_y)



   if derivative_degree_x==0&&derivative_degree_y==0
        
        if basis_index==1
            r=1-x-y;
        elseif basis_index==2
            r=x;
        elseif basis_index==3
            r=y;
        end

    elseif derivative_degree_x==1&&derivative_degree_y==0
        
        if basis_index==1
            r=-1;
        elseif basis_index==2
            r=1;
        elseif basis_index==3
            r=0;
        end

    elseif derivative_degree_x==0&&derivative_degree_y==1
        
        if basis_index==1
            r=-1;
        elseif basis_index==2
            r=0;
        elseif basis_index==3
            r=1;
        end
        
    end