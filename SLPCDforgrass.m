clear
kappa_grid=0:0.002:2;
Z=imread('grass.tiff');
Z=im2double(Z);
Z=Z(1:256,1:256);
observation=Z;
m=512/2;n=512/2;
N=[m,n];
q=2;
%chess=zeros(m,n);
%chess=conschess(chess);
G_start=q;
%% priors and starting points
mu_0=mean(mean(observation)); %prior of mean parameter is normal distribution
mu=normrnd(mu_0,0.2,[1,q]); %starting points
sigma=(unifrnd(0,0.1,[1,q])).^2;
iteration=6000;
burnin=2000;
alpha=1;tau=0.01;
data=observation;
clear('observation');
Z_start=randi(G_start,N);
% pi_start=sample_dirichlet(ones(1,G_start),1);
% pi_start=sample_dirichlet(ones(1,G_start),1);
kappa_start=0.3;
R=max(data(:))-min(data(:));
invp1=0.001;invp2=0.001;
mu_sig=100^2;
%% Chess
threedeye=eye(2);
chess = repmat(threedeye,N/2);
kappa_upper=4;
Zshell=Z;
Zshell(:)=1:m*n;

%% (1 item)
p1=dlmread(strcat('p1tableCT=1','q=',num2str(q),'.txt'));
%% Conditional type(2 items)
nCT=min(2,q);
p2=[];
for i=1:nCT
    p2(:,:,i)=dlmread(strcat('p2tableCT=',num2str(i),'q=',num2str(q),'.txt'));
end

p1(:,2:end)=log(p1(:,2:end));
p2(:,2:end,:)=log(p2(:,2:end,:));
p1(p1==-Inf)=0;
p2(p2==-Inf)=0;
iteration=6000;
burnin=2000;
%% MCMC
for iter=1:iteration
    chess=1-chess;
    %% update Z
    Z_new=(1-chess).*Z_start + chess.*prop_new_potts(Z_start,q); %propose new value
    neighbors_new = findneighbour(Z_new);
    neighbors_start =findneighbour(Z_start);
    likelihood_start=log(normpdf(data,mu(Z_start),sqrt(sigma(Z_start))));
    likelihood_new=log(normpdf(data,mu(Z_new),sqrt(sigma(Z_new))));
    Z_prob=likelihood_new-likelihood_start+...
        (neighbors_new-neighbors_start)*kappa_start;
    transitions = (log(rand(N)) < Z_prob ).* chess .* (Z_new - Z_start);
    Z_start=Z_start+transitions;
    chess=1-chess;
    %% update Z
    Z_new=((1-chess).*Z_start+chess.*prop_new_potts(Z_start,q)); %propose new value
    neighbors_new = findneighbour(Z_new);
    neighbors_start =findneighbour(Z_start);
    likelihood_start=log(normpdf(data,mu(Z_start),sqrt(sigma(Z_start))));
    likelihood_new=log(normpdf(data,mu(Z_new),sqrt(sigma(Z_new))));
    Z_prob=likelihood_new-likelihood_start+...
        (neighbors_new-neighbors_start)*kappa_start;
    transitions = (log(rand(N)) < Z_prob ).* chess .* (Z_new - Z_start);
    Z_start=Z_start+transitions;
    
    %% ============ update kappa
    d=(Z_start==circshift(Z_start, [1 0]));
    cnt1=d(2:end,1)+0;
    d=(Z_start==circshift(Z_start, [0 1]));
    cnt2=d(1,2:end)+0;
    tem=tabulate([cnt1;cnt2']);
    counts1=tem(:,2)';
    d=(Z_start==circshift(Z_start, [1 0]))+0;
    d=d+(Z_start==circshift(Z_start, [0 1]));
    d(1,:)=[];
    d(:,1)=[];
    nCT=min(2,q); CTmat=zeros(size(Z_start)-1);
    for ii=2:size(Z_start,1)
        for jj=2:size(Z_start,2)
            neib=[Z_start(minuspos(ii,m),jj),Z_start(ii,minuspos(jj,n))];
            CTmat(ii-1,jj-1)= length(unique(neib));
        end
    end
    
    p2tableCT=zeros(nCT,2+1);
    for ii=1:nCT
        tem=tabulate(d(CTmat==ii));
        p2tableCT(ii,tem(:,1)+1)=p2tableCT(ii,tem(:,1)+1)+tem(:,2)';
    end
    
    kappa_new=kappa_start+0.02*randn(1);
    if kappa_new >0 && kappa_new<kappa_upper
        kappa_start_r=abs(kappa_start-kappa_grid);
        kappa_new_r=abs(kappa_new-kappa_grid);
        rnk_start=find(kappa_start_r==min(kappa_start_r),1);
        rnk_new=find(kappa_new_r==min(kappa_new_r),1);
        like_p1start=sum(counts1.*p1(rnk_start,2:end));
        like_p1new=sum(counts1.*p1(rnk_new,2:end));
        like_p2start=0;
        for ii=1:nCT
            like_p2start=like_p2start+sum(p2tableCT(ii,:).*p2(rnk_start,2:4,ii));
        end
        like_p2new=0;
        for ii=1:nCT
            like_p2new=like_p2new+sum(p2tableCT(ii,:).*p2(rnk_new,2:4,ii));
        end
        
        like_new=like_p1new+like_p2new;
        like_start=like_p1start+like_p2start;
        if like_new-like_start>log(rand(1))
            kappa_start=kappa_new;
        end
        
    end

    %% update mu
    for i=1:G_start
        a=(Z_start==i);
        aa=sum(a(:));
        block=data(a);
        mu_new=normrnd(mu(i),0.002);
        mu_prob=sum(log(normpdf(block,mu_new,sqrt(sigma(i))))-...
            log(normpdf(block,mu(i),sqrt(sigma(i)))))...
            +log(normpdf(mu_new,mu_0,sqrt(mu_sig)))-log(normpdf(mu(i),mu_0,sqrt(mu_sig)));
        if mu_prob>log(rand(1))
            mu(i)=mu_new;
        end
%         ave=(sum(block)+tau*mu_0)/(aa+tau);
%         sig=sigma(i)/(aa+tau);
%         mu(i)=normrnd(ave,sqrt(sig));
    end
    %% update of mean parameter is over
    for i=1:G_start
        a=(Z_start==i);
        block=data(a);
        sigma_new=(normrnd(sqrt(sigma(i)),0.001))^2;
        sigma_prob=sum(log(normpdf(block,mu(i),sqrt(sigma_new)))-log(normpdf(block,mu(i),sqrt(sigma(i)))))...
            +inversegampdf(sigma_new,invp1,invp2)-inversegampdf(sigma(i),invp1,invp2);
        if sigma_prob>log(rand(1))
            sigma(i)=sigma_new;
        end
    end
%     sigma=1./phi;
%% update of sigma is over. need to update \kappa
mu_mat(iter,:)=mu;
sigma_mat(iter,:)=sigma;
kappa_mat(iter)=kappa_start;
    if iter>burnin
        Z_cell{iter-burnin}=Z_start;
    end
end
save(strcat('SLPCDforgrass','q=',num2str(q),'.mat'));
exit;

plot(kappa_mat(2001:end))

mean(mu_mat(2001:end,:))
std(mu_mat(2001:end,:))
mean(sigma_mat(2001:end,:))
std(sigma_mat(2001:end,:))
mean(kappa_mat(2001:end))
std(kappa_mat(2001:end))

MSE=0;
for i=1:4000
    mu_sim=mu_mat(i+2000,:);
    sigma_sim=sqrt(sigma_mat(i+2000,:));
    Zsim=normrnd(mu_sim(Z_cell{i}),sigma(Z_cell{i}));
    MSE=MSE+mean((Zsim(:)-Z(:)).^2);
end
MSE=MSE/4000

imshow(Z_start-1)
set(gca,'position',[0,0,1,1])

set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'grass_SLPCD', 'jpg') %Save figure

