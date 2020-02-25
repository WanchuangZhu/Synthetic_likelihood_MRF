function vec=findstructure(x,y,Z,indicator)
vec=[];
if indicator==1
    
    % vec=zeros(1,G_start);
    if x-1>0 && Z(x-1,y)~=0
        % vec(Z(x-1,y))=vec(Z(x-1,y))+1;
        vec=[vec,Z(x-1,y)];
    end
    if x+1<=size(Z,1) && Z(x+1,y)~=0
        % vec(Z(x+1,y))=vec(Z(x+1,y))+1;
        vec=[vec,Z(x+1,y)];
    end
    if y-1>0 && Z(x,y-1)~=0
        % vec(Z(x,y-1))=vec(Z(x,y-1))+1;
        vec=[vec,Z( x,y-1)];
    end
    if y+1<=size(Z,2) && Z(x,y+1)~=0
        % vec(Z(x,y+1))=vec(Z(x,y+1))+1;
        vec=[vec,Z(x,y+1)];
    end
elseif indicator==0
   
    % vec=zeros(1,G_start);
    if x-1>0 && y-1>0 && Z(x-1,y-1)~=0
        % vec(Z(x-1,y-1))=vec(Z(x-1,y-1))+1;
        vec=[vec,Z( x-1,y-1)];
    end
    if x-1>0 && y+1<=size(Z,2) && Z(x-1,y+1)~=0
        % vec(Z(x-1,y+1))=vec(Z(x-1,y+1))+1;
        vec=[vec,Z( x-1,y+1)];
    end
    if x+1<=size(Z,1) && y-1>0 && Z(x+1,y-1)~=0
        % vec(Z(x+1,y-1))=vec(Z(x+1,y-1))+1;
        vec=[vec,Z( x+1,y-1)];
    end
    if x+1<=size(Z,1) && y+1<=size(Z,2) && Z(x+1,y+1)~=0
        % vec(Z(x+1,y+1))=vec(Z(x+1,y+1))+1;
        vec=[vec,Z( x+1,y+1)];
    end
end
end