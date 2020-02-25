clear
%% simulate dataset
kappa=0.1:0.1:0.3;
kappa_grid=0:0.002:0.5;
q=2;
jobid=getenv('PBS_ARRAYID');
jobid=str2double(jobid);
lattsize=[32 128 256];
m=lattsize(jobid);n=m;
%% read in the pretable
%% Conditional type(1 item)
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

%% Conditional type(3 items)
nCT=min(3,q);
p3=[];
for i=1:nCT
    p3(:,:,i)=dlmread(strcat('p3tableCT=',num2str(i),'q=',num2str(q),'.txt'));
end
p3(:,2:end,:)=log(p3(:,2:end,:));
p3(p3==-Inf)=0;
%% Conditional type(4 items)
nCT=3;
p4=[];
for i=1:nCT
    p4(:,:,i)=dlmread(strcat('p4tableCT=',num2str(i),'q=',num2str(q),'.txt'));
end
p4(:,2:end,:)=log(p4(:,2:end,:));
p4(p4==-Inf)=0;
% 
% walksteps=[0.1,0.1,0.095,0.085,0.085,0.07,0.06,0.045;...
%     0.04,0.04,0.035,0.035,0.035,0.02,0.01,0.01;...
%     0.01,0.01,0.015,0.01,0.01,0.0025,0.0025,0.0025;];
parpool('local',3)
spmd
    kappa_true=kappa(labindex);
    iteration=6000;
    reps=200;
    for rep=1:reps
        fadd=fullfile('snddata',strcat(num2str(m),num2str(q)),...
            strcat('Znd',num2str(kappa_true),num2str(rep),'.txt'));
        Z=dlmread(fadd);
        %% starting point
        kappa_start=0.15;
        %% count the pairs (1 item)
        d=(Z==circshift(Z, [1 0]));
        cnt1=d(2:end,1)+0;

        tem=tabulate(cnt1);
        counts1=tem(:,2)';
        %% count the pairs (2 item)
        d=(Z==circshift(Z, [0 1]));
        d=d+(Z==circshift(Z, [-1 1]));
        cnt2=d(1,2:end)+0;
        
        nCT=min(2,q); CTmat=zeros(size(Z));
        for ii=1:1
            for jj=1:size(Z,2)
                neib=[Z(addpos(ii,m),minuspos(jj,n)),Z(ii,minuspos(jj,n))];
                CTmat(ii,jj)= length(unique(neib));
            end
        end
        p2tableCT=zeros(nCT,2+1);
        for ii=1:nCT
            tem=tabulate(cnt2(CTmat(1,2:end)==ii));
            p2tableCT(ii,tem(:,1)+1)=p2tableCT(ii,tem(:,1)+1)+tem(:,2)';
        end
        %% count the pairs (3 item)
        d=(Z==circshift(Z, [0 1]));
        d=d+(Z==circshift(Z, [1 1]));
        d=d+(Z==circshift(Z, [1 0]));
        cnt3=d(end,2:end)+0;
        
        nCT=min(3,q); CTmat=zeros(size(Z));
        for ii=1:size(Z,1)
            for jj=size(Z,2):size(Z,2)
                neib=[Z(minuspos(ii,m),minuspos(jj,n)),Z(ii,minuspos(jj,n)),Z(minuspos(ii,m),jj)];
                CTmat(ii,jj)= length(unique(neib));
            end
        end
        p3tableCT=zeros(nCT,3+1);
        for ii=1:nCT
            tem=tabulate(cnt3(CTmat(end,2:end)==ii));
            if isempty(tem)
                continue;
            end
            p3tableCT(ii,tem(:,1)+1)=p3tableCT(ii,tem(:,1)+1)+tem(:,2)';
        end
        %% count the pairs (4 item)
        d=(Z==circshift(Z, [0 1]));
        d=d+(Z==circshift(Z, [1 1]));
        d=d+(Z==circshift(Z, [1 0]));
        d=d+(Z==circshift(Z, [-1 1]));
        cnt4=d(2:end-1,2:end-1)+0;
        
        nCT=3; CTmat=zeros(size(Z));
        for ii=1:size(Z,1)
            for jj=1:size(Z,2)
                neib=[Z(minuspos(ii,m),minuspos(jj,n)),Z(ii,minuspos(jj,n)),...
                    Z(minuspos(ii,m),jj),Z(addpos(ii,m),minuspos(jj,n))];
                CTmat(ii,jj)= length(unique(neib));
                if sum(neib==2)==2
                    CTmat(ii,jj)=3;
                end
            end
        end
        p4tableCT=zeros(nCT,4+1);
        for ii=1:nCT
            tem=tabulate(cnt4(CTmat(2:end-1,2:end-1)==ii));
            if isempty(tem)
                continue;
            end
            p4tableCT(ii,tem(:,1)+1)=p4tableCT(ii,tem(:,1)+1)+tem(:,2)';
        end
        %% find the optimal sigma of proposal distribution
        sigma_scale_start=0.1; % simga for proposal distribution
        iMax=2000; % if in 100 iterations restart is not happening we find the sigma
        pstar=0.44; % optimal acceptance rate
        n0=round(5/(pstar*(1-pstar)));
        i=1;
        sigma_vec(1)=sigma_scale_start;
        sigma_scale=sigma_scale_start;
        
         while i<iMax
            kappa_new=normrnd(kappa_start,sigma_scale);
            if kappa_new>0  && kappa_new<0.5
                kappa_start_r=abs(kappa_start-kappa_grid);
                kappa_new_r=abs(kappa_new-kappa_grid);
                rnk_start=find(kappa_start_r==min(kappa_start_r),1);
                rnk_new=find(kappa_new_r==min(kappa_new_r),1);
                %% 1 item
                like_p1start=sum(counts1.*p1(rnk_start,2:end));
                like_p1new=sum(counts1.*p1(rnk_new,2:end));
                %% 2 items
                nCT=min(2,q);
                like_p2start=0;
                for ii=1:nCT
                    like_p2start=like_p2start+sum(p2tableCT(ii,:).*p2(rnk_start,2:end,ii));
                end
                like_p2new=0;
                for ii=1:nCT
                    like_p2new=like_p2new+sum(p2tableCT(ii,:).*p2(rnk_new,2:end,ii));
                end
                %% 3 items
                nCT=min(3,q);
                like_p3start=0;
                for ii=1:nCT
                    like_p3start=like_p3start+sum(p3tableCT(ii,:).*p3(rnk_start,2:end,ii));
                end
                like_p3new=0;
                for ii=1:nCT
                    like_p3new=like_p3new+sum(p3tableCT(ii,:).*p3(rnk_new,2:end,ii));
                end
                %% 4 items
                nCT=3;
                like_p4start=0;
                for ii=1:nCT
                    like_p4start=like_p4start+sum(p4tableCT(ii,:).*p4(rnk_start,2:end,ii));
                end
                like_p4new=0;
                for ii=1:nCT
                    like_p4new=like_p4new+sum(p4tableCT(ii,:).*p4(rnk_new,2:end,ii));
                end
                
                like_new=like_p1new+like_p2new+like_p3new+like_p4new;
                like_start=like_p1start+like_p2start+like_p3start+like_p4start;
                if like_new-like_start>log(rand(1))
                    kappa_start=kappa_new;
                    if i>n0
                        sigma_scale=sigma_scale + sigma_scale/(pstar*i);
                        sigma_vec=[sigma_vec,sigma_scale];
                    end
                else
                    if i>n0
                        sigma_scale=sigma_scale-sigma_scale/((1-pstar)*i);
                        sigma_vec=[sigma_vec,sigma_scale];
                    end
                end
            else continue;
            end
            
            if i==n0
                sigma_scale_start=sigma_scale;
            end
            
            if i>n0
                Toobig=sigma_scale>3*sigma_scale_start;
                Toosmall=sigma_scale < sigma_scale_start/3;
                if Toobig || Toosmall
                    i=1;sigma_vec=[];continue;
                end
            end
            i=i+1;
        end
        sigma_scale=mean(sigma_vec(end-400:end));
        %% MCMC
        for iter=1:iteration
            kappa_new=normrnd(kappa_start,sigma_scale);
            if kappa_new>0 && kappa_new<0.5
                kappa_start_r=abs(kappa_start-kappa_grid);
                kappa_new_r=abs(kappa_new-kappa_grid);
                rnk_start=find(kappa_start_r==min(kappa_start_r),1);
                rnk_new=find(kappa_new_r==min(kappa_new_r),1);
                %% 1 item
                like_p1start=sum(counts1.*p1(rnk_start,2:end));
                like_p1new=sum(counts1.*p1(rnk_new,2:end));
                %% 2 items
                nCT=min(2,q);
                like_p2start=0;
                for ii=1:nCT
                    like_p2start=like_p2start+sum(p2tableCT(ii,:).*p2(rnk_start,2:end,ii));
                end
                like_p2new=0;
                for ii=1:nCT
                    like_p2new=like_p2new+sum(p2tableCT(ii,:).*p2(rnk_new,2:end,ii));
                end
                %% 3 items
                nCT=min(3,q);
                like_p3start=0;
                for ii=1:nCT
                    like_p3start=like_p3start+sum(p3tableCT(ii,:).*p3(rnk_start,2:end,ii));
                end
                like_p3new=0;
                for ii=1:nCT
                    like_p3new=like_p3new+sum(p3tableCT(ii,:).*p3(rnk_new,2:end,ii));
                end
                %% 4 items
                nCT=3;
                like_p4start=0;
                for ii=1:nCT
                    like_p4start=like_p4start+sum(p4tableCT(ii,:).*p4(rnk_start,2:end,ii));
                end
                like_p4new=0;
                for ii=1:nCT
                    like_p4new=like_p4new+sum(p4tableCT(ii,:).*p4(rnk_new,2:end,ii));
                end
                
                like_new=like_p1new+like_p2new+like_p3new+like_p4new;
                like_start=like_p1start+like_p2start+like_p3start+like_p4start;
                if like_new-like_start>log(rand(1))
                    kappa_start=kappa_new;
                end
            end
            kappa_vec(iter)=kappa_start;
        end
        kappa_mat(rep,:)=kappa_vec;
    end
    % kappa_cell{kappanum}=kappa_mat;
end
dsave(strcat('MCAPCD',num2str(m),'q=',num2str(q)','.mat'));
exit


%%%%%%
%%%%%%
for i=1:3
    kappa_mmat=kappa_mat{i};
    kappa_mmat=kappa_mmat(:,2001:end);
    sd(i)=mean(std(kappa_mmat'));
    kappa_point(i,:)=mean(kappa_mmat');
    esti(i)=mean(mean(kappa_mmat'));
end
for i=1:3
    kappa_rmse(i)=sqrt(mean((kappa_point(i,:)-0.1*i).^2));
end
latex(esti,'%.3f','nomath')
latex(kappa_rmse,'%.3f','nomath')

latex(esti,'%.3f','nomath')
latex(sd,'%.3f','nomath')