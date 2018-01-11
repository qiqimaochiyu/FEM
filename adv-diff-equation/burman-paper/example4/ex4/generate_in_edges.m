function [Boundary,T_plus,T_index]=generate_in_edges(T,tnp)


t=0;
for j=1:tnp
    [~,n]=find(T==T(1,j));
    [~,q]=find(T==T(2,j));
    [~,l]=find(T==T(3,j));
    a=intersect(n,q);
    b=intersect(q,l);
    c=intersect(l,n);
    if length(a)==2 && max(a)~=j
        t=t+1;
        Boundary(1,t)=j;
        Boundary(2,t)=max(a);
        Boundary(3,t)=T(1,j);
        Boundary(4,t)=T(2,j);
        T_plus(1,t)=T(1,j);
        T_plus(2,t)=T(2,j);
        T_plus(3,t)=T(3,j);
        T_plus(4,t)=T(setdiff(setdiff([1 2 3],find(T(:,max(a))==T(1,j))),find(T(:,max(a))==T(2,j))),max(a));
        T_index(1,t)=1;
        T_index(2,t)=find(T(:,max(a))==T(1,j));
        T_index(3,t)=2;
        T_index(4,t)=find(T(:,max(a))==T(2,j));
        T_index(5,t)=3;
        T_index(6,t)=setdiff(setdiff([1 2 3],find(T(:,max(a))==T(1,j))),find(T(:,max(a))==T(2,j)));
    end
    
    if length(b)==2 && max(b)~=j
        t=t+1;
        Boundary(1,t)=j;
        Boundary(2,t)=max(b);
        Boundary(3,t)=T(2,j);
        Boundary(4,t)=T(3,j);
        T_plus(1,t)=T(2,j);
        T_plus(2,t)=T(3,j);
        T_plus(3,t)=T(1,j);
        T_plus(4,t)=T(setdiff(setdiff([1 2 3],find(T(:,max(b))==T(2,j))),find(T(:,max(b))==T(3,j))),max(b));
        T_index(1,t)=2;
        T_index(2,t)=find(T(:,max(b))==T(2,j));
        T_index(3,t)=3;
        T_index(4,t)=find(T(:,max(b))==T(3,j));
        T_index(5,t)=1;
        T_index(6,t)=setdiff(setdiff([1 2 3],find(T(:,max(b))==T(2,j))),find(T(:,max(b))==T(3,j)));
        
    end
    
    if length(c)==2 && max(c)~=j
        t=t+1;
        Boundary(1,t)=j;
        Boundary(2,t)=max(c);
        Boundary(3,t)=T(3,j);
        Boundary(4,t)=T(1,j);
        T_plus(1,t)=T(3,j);
        T_plus(2,t)=T(1,j);
        T_plus(3,t)=T(2,j);
        T_plus(4,t)=T(setdiff(setdiff([1 2 3],find(T(:,max(c))==T(3,j))),find(T(:,max(c))==T(1,j))),max(c));
        T_index(1,t)=3;
        T_index(2,t)=find(T(:,max(c))==T(3,j));
        T_index(3,t)=1;
        T_index(4,t)=find(T(:,max(c))==T(1,j));
        T_index(5,t)=2;
        T_index(6,t)=setdiff(setdiff([1 2 3],find(T(:,max(c))==T(3,j))),find(T(:,max(c))==T(1,j)));
    end
        
        
end

     
    