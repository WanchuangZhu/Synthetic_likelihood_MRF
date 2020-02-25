function [neibcell,Z_res] = neibstructure(Z,times)
%composedecom returns the matrix block for each split using RCoDA.
N=size(Z);
if mod(N(1),2)==1
    Z=[Z;zeros(1,N(2))];
end
if mod(N(2),2)==1
    Z=[Z,zeros(N(1),1)];
end
chess=repmat(eye(2),N/2); %% creat the chess board
chess=chess(1:N(1),1:N(2));
mask=(Z>0);

for t=1:times
    n=0;v=[];
    if mod(t,2)==1
        for i=1:N(1)
            for j=1:N(2)
                if chess(i,j)==1 && mask(i,j)==1
                    
                    n=n+1;
                    
                    nei=findstructure(i,j,Z,1);
                    v(n,1:length(nei)+1)=[Z(i,j),nei];
                end
            end
        end
        chess(:,1:2:(N(2)-1))=1;
    elseif mod(t,2)==0
        for i=1:N(1)
            for j=1:N(2)
                if chess(i,j)==0 && mask(i,j)==1
                    n=n+1;
                    nei=findstructure(i,j,Z,0);
                    v(n,1:length(nei)+1)=[Z(i,j),nei];
                end
            end
        end
        Z=Z(2:2:N(1),1:2:(N(2)-1));
        N=size(Z);
        if mod(N(1),2)==1
            Z=[Z;zeros(1,N(2))];
            N=size(Z);
        end
        if mod(N(2),2)==1
            Z=[Z,zeros(N(1),1)];
            N=size(Z);
        end
        mask=(Z>0);
        if mod(N(1),2)==0 && mod(N(2),2)==0
            chess=repmat(eye(2),N/2);
        end
    end
    neibcell{t}=v;
end
Z_res=Z;
end

