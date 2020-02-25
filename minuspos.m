function tem=minuspos(x,N)
%% x is row number or column number
%% N is corresponding total row number or total column number
tem=x-1;
if x==1
    tem=N;
end
