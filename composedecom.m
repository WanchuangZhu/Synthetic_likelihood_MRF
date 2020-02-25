function [y,Z_res] = composedecom(Z,times,G_start )
%composedecom returns the matrix block for each split using RCoDA.
N=size(Z);
chess=repmat(eye(2),N/2);
mask=(Z>0);
for t=1:times
    n=0;v=[];
    if mod(t,2)==1
        for i=1:N(1)
            for j=1:N(2)
                if chess(i,j)==1 && mask(i,j)==1
                    n=n+1;
                    nei=showneibou(i,j,Z,G_start,1);
                    v(n,:)=[Z(i,j),nei];
                end                
            end
        end
        chess(:,1:2:(N(2)-1))=1;
    else
        for i=1:N(1)
            for j=1:N(2)
                if chess(i,j)==0 && mask(i,j)==1
                    n=n+1;
                    nei=showneibou(i,j,Z,G_start,0);
                    v(n,:)=[Z(i,j),nei];
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
    y{t}=v;
end
Z_res=Z;
end

