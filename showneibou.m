function vec=showneibou(x,y,Z,G_start,indicator)
% this function is used to show the constitution of a specific pixel's neighbourhood.
% x,y is the coordinates in the field. Z is the field. G_start is number of the components
% for each pixel, this funtion returns a vec which shows the constitution of its neibourghhood
% indicator shows which calculation we should use

if indicator==1
    vec=zeros(1,G_start);
    if x-1>0 && Z(x-1,y)~=0
        vec(Z(x-1,y))=vec(Z(x-1,y))+1;
    end
    if x+1<=size(Z,1) && Z(x+1,y)~=0
        vec(Z(x+1,y))=vec(Z(x+1,y))+1;
    end
    if y-1>0 && Z(x,y-1)~=0
        vec(Z(x,y-1))=vec(Z(x,y-1))+1;
    end
    if y+1<=size(Z,2) && Z(x,y+1)~=0
        vec(Z(x,y+1))=vec(Z(x,y+1))+1;
    end
elseif indicator==0
    vec=zeros(1,G_start);
    if x-1>0 && y-1>0 && Z(x-1,y-1)~=0
        vec(Z(x-1,y-1))=vec(Z(x-1,y-1))+1;
    end
    if x-1>0 && y+1<=size(Z,2) && Z(x-1,y+1)~=0
        vec(Z(x-1,y+1))=vec(Z(x-1,y+1))+1;
    end
    if x+1<=size(Z,1) && y-1>0 && Z(x+1,y-1)~=0
        vec(Z(x+1,y-1))=vec(Z(x+1,y-1))+1;
    end
    if x+1<=size(Z,1) && y+1<=size(Z,2) && Z(x+1,y+1)~=0
        vec(Z(x+1,y+1))=vec(Z(x+1,y+1))+1;
    end
end
end
        