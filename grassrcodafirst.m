clear
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

times=12;
alpha_start=0.3;

R=max(data(:))-min(data(:));
invp1=0.001;invp2=0.001;
mu_sig=100^2;
% gammap1=4; gammap2=R^2/100; %priors for phi
%% Chess
threedeye=eye(2);
chess = repmat(threedeye,N/2);
kappa_upper=4;
% parpool('local', 2)
%% creat the neighbourhood structure. It returns index of each pixel and related neighbours' indices.
% this can be used to update kappa and alpha
Zshell=Z;
Zshell(:)=1:m*n;

[neibindex,Z_resind] = neibstructure(Zshell,times);
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
    
    [neibcell,~]=composedecom(Z_start,times,q); % find the matrix blocks using RCoDA
    %% ============ update kappa
    Z_vec=Z_start(:);
    Nhat=N(1)/(2^(times/2));
    file=strcat('NCwhole_',num2str(Nhat),'times',num2str(Nhat),num2str(q),'.txt');
    NCserial=importdata(file);
    NCserial=NCserial.data;
    
    kappa_new=kappa_start+0.02*randn(1);
    % [~,Z_res]=composelike(kappa_start,Z_start,times,G_start,alpha_start);
    Z_res=Z_start(Z_resind);
    pairs=0;
    d=(Z_res==circshift(Z_res, [0 1]));
    d(:,1)=0;
    pairs=pairs+sum(sum(d));
    d=(Z_res==circshift(Z_res, [0 -1]));
    d(:,size(Z_res,2))=0;
    pairs=pairs+sum(sum(d));
    d=(Z_res==circshift(Z_res, [1 0]));
    d(1,:)=0;
    pairs=pairs+sum(sum(d));
    d=(Z_res==circshift(Z_res, [-1 0]));
    d(size(Z_res,1),:)=0;
    pairs=pairs+sum(sum(d));
    pairs=pairs/2;
    if kappa_new >0 && kappa_new<kappa_upper
        kappa_serial_start=kappa_start*(alpha_start.^(0:(times-1)));
        kappa_serial_new=kappa_new*(alpha_start.^(0:(times-1)));
        like_start=RCoDAlike(neibcell,kappa_serial_start,times,q);
        % like_start=RCoDAlikestructure(neibindex,kappa_serial_start,times,q,Z_vec);
        like_new=RCoDAlike(neibcell,kappa_serial_new,times,q);
        kappa_prob=like_new-like_start+pairs*(kappa_new-kappa_start)*alpha_start^(times)...
            +ncintegnew(alpha_start^(times)*kappa_start,NCserial) -...
            ncintegnew(alpha_start^(times)*kappa_new,NCserial);
        if kappa_prob>log(rand(1))
            kappa_start=kappa_new;
        end
    end
    %% update alpha
    alpha_new=normrnd(alpha_start,0.01);
    if alpha_new>0 && alpha_new<1
        kappa_serial_start=kappa_start*(alpha_start.^(0:(times-1)));
        kappa_serial_new=kappa_start*alpha_new.^(0:(times-1));
        like_alpha_start=RCoDAlike(neibcell,kappa_serial_start,times,q);
        like_alpha_new=RCoDAlike(neibcell,kappa_serial_new,times,q);
        
        alpha_prob=like_alpha_new-like_alpha_start + ...
            pairs*kappa_start*(alpha_new^(times)-alpha_start^(times))+...
            ncintegnew(alpha_start^(times)*kappa_start,NCserial) - ...
            ncintegnew(alpha_new^(times)*kappa_start,NCserial);
        if alpha_prob>log(rand(1))
            alpha_start=alpha_new;
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
%         shape=(gammap1+length(a))/2;
%         scale=2/(gammap2+sum(((data(Z_start==i)-mu(i))).^2));
%         phi(i)=gamrnd(shape,scale);
    end
%     sigma=1./phi;
%% update of sigma is over. need to update \kappa
mu_mat(iter,:)=mu;
sigma_mat(iter,:)=sigma;
kappa_mat(iter)=kappa_start;
alpha_mat(iter)=alpha_start;
    if iter>burnin
        Z_cell{iter-burnin}=Z_start;
    end
end
save(strcat('grassrcodafirst','q=',num2str(q),'.mat'));
exit;
%%%%%%%%%%%%%%%%%%%%%%%%%%%   label switching
for i=1:size(mu_mat,1)
    [b,arr]=sort(mu_mat(i,:));
    mu_mat(i,:)=b;
    sigma_mat(i,:)=sigma_mat(i,arr);
end
mu_est=mean(mu_mat(2000:end,:));
sigma_est=mean(sigma_mat(2000:end,:));
mean(kappa_mat(2000:end))
sqrt(var(mu_mat(1:end,:)))
sqrt(var(sigma_mat(1:end,:)))
sqrt(var(kappa_mat(2:end)))
%% MAP of Z
postmean=zeros(N);
for i=1:size(data,1)
    for j=1:size(data,2)
           for k=1:G_start
               v(k)=normpdf(data(i,j),mu_est(k),sqrt(sigma_est(k)));
            alloc=find(v==max(v));
           end
        postmean(i,j)=alloc;
    end
end
imshow(3-postmean,[1 2])
set(gca,'position',[0 0.02 1 0.96])
%% posterior mean of Z

%% posterior prediction
postpredict=zeros(N);
for i=1:size(Z_cell,2)
    for j=1:N(1)
        for k=1:N(2)
            alloc=Z_cell{i}(j,k);
            postsim(j,k)=normrnd(mu_mat(i,alloc),sqrt(sigma_mat(i,alloc)));
        end
    end
    postpredict=postpredict+postsim;
end
postpredict=postpredict/size(Z_cell,2);
% distance between posterior prediction and data
distant=postpredict-data;
imshow(abs(distant))
sum(sum(abs(distant))) % =
%% posterior mode
mu_est=mean(mu_mat);
sigma_est=mean(sigma_mat);
tem=Z_cell{1};
for i=2:3000
    tem(:,:,i)=Z_cell{i};  
end
Z_t=zeros(256,256);
for i=1:256
    for j=1:256
        t=tem(i,j,:);
        
        w_mat1(i,j)=sum(t==1)/3000;
        w_mat2(i,j)=sum(t==2)/3000;
%         if prob(1)*normpdf(data(i,j),mu_est(1),sqrt(
%             Z_t(i,j)=1;
%         else Z_t(i,j)=2;
%         end
    end
end
for i=1:4000
    mu_point=mu_mat(i,:);
    sigma_point=sigma_mat(i,:);
    Z_sim(:,:,i)=w_mat1.*normrnd(mu_point(1),sqrt(sigma_point(1)))+...
        w_mat2.*normrnd(mu_point(2),sqrt(sigma_point(2)));
end
for i=1:256
    for j=1:256
        if Z(i,j)<=quantile(Z_sim(i,j,:),0.975) && Z(i,j)>= quantile(Z_sim(i,j,:),0.025)
            Z_in(i,j)=1;
        else Z_in(i,j)=0;
        end
    end
end
sum(Z_in(:))/256/256  %0.9850
imshow(Z_t)
set(gca,'position',[0 0.02 1 0.96])
%% compare with truth
load('grassrcodafirstq=2.mat')
mu_mat=mu_mat(2001:end,:);
sigma_mat=sigma_mat(2001:end,:);
for i=1:4000
    ii=i;
    z_tem=Z_cell{ii};
    mu=mu_mat(ii,:);
    sigma=sigma_mat(ii,:);
    obs=normrnd(mu(z_tem),sqrt(sigma(z_tem)));
    Z_com(:,:,i)=obs;
end
Z_est=mean(Z_com,3);
sum(sum(abs(Z_est-data)))/256/256

for i=1:4000
    ii=i;
    z_tem=Z_cell{ii};
    mu=mu_mat(ii,:);
    sigma=sigma_mat(ii,:);
    obs=normrnd(mu(z_tem),sqrt(sigma(z_tem)));
    Z_com(:,:,i)=(obs-data).^2;
end
Z_est=mean(Z_com,3);
sum(sum((Z_est)))/256/256 %0.0328