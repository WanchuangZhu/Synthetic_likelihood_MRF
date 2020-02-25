function num=findneighbour(x)
% x is allocation variable matrix
num=(x==circshift(x, [ 0 1])) + ... % compare with  left neighbours
    (x==circshift(x, [ 0 -1])) + ...%compare with  right neighbours
    (x==circshift(x, [ 1 0])) + ... %compare with  backward neighbours
    (x==circshift(x, [ -1 0]));
num=0;
d=(x==circshift(x, [0 1]));
d(:,1)=0;
num=num+d;
d=(x==circshift(x, [0 -1]));
d(:,size(x,2))=0;
num=num+d;
d=(x==circshift(x, [1 0]));
d(1,:)=0;
num=num+d;
d=(x==circshift(x, [-1 0]));
d(size(x,1),:)=0;
num=num+d;
end