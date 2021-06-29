
%%%%%%%%%%%%%%%%%%%%%%%% Chaotic AOA (CAOA)

clear;
close all;
clc;



nPop=100;
Maxit=20;
nVar=2; %% number of optimization parameters
C1=2;C2=6; C3=1;C4=2; u=.9;l=.1;   %%% The paramters of the algorithm


%%%% Initialization

lb=-4.5*ones(1,nVar); %%% same search area for different parameters
ub=4.5*ones(1,nVar);

% lb=[0*ones(1,D-1),0.000001]; %% different search area for each parameter
% ub=[5*ones(1,D-1),0.5];





for i=1:nPop
    X(i,:)=lb+rand(1,nVar).*(ub-lb);
    Y(i)=cost(X(i,:));
    den(i,:)=rand(1,nVar); 
    vol(i,:)=rand(1,nVar);
    acc(i,:)=lb+rand(1,nVar).*(ub-lb);
end



[Scorebest, Score_index] = min(Y);
Xbest = X(Score_index,:);
den_best=den(Score_index,:);
vol_best=vol(Score_index,:);
acc_best=acc(Score_index,:);
acc_norm=acc;

chnum=1;  %%% selection of chaotic map  chnum=1,...,10
chstart=0.65;  %%% starting point for chaotic map
chmap=ch_map(chnum,chstart,Maxit);


for it = 1:Maxit
    
    TF(it)=exp(((it-Maxit)/(Maxit)))+(-0.2+0.4*chmap(it));  %%%% TF + chaotic_value
    
%     TF(t)=exp(((t-Maxit)/(Maxit))); %%%% Basic algorithm


    if TF(it)>1
        TF(it)=1;
    end
    d=exp((Maxit-it)/Maxit)-(it/Maxit); 
    acc=acc_norm;
    r=rand();
    for i=1:nPop
        den(i,:)=den(i,:)+r*(den_best-den(i,:));   
        vol(i,:)=vol(i,:)+r*(vol_best-vol(i,:));
        if TF(it) < 0.45     
            mr=randi(nPop);
            acc_temp(i,:)=(den(mr,:)+(vol(mr,:).*acc(mr,:)))./(rand*den(i,:).*vol(i,:));   
        else
            acc_temp(i,:)=(den_best+(vol_best.*acc_best))./(rand*den(i,:).*vol(i,:));   
        end
    end
    
    acc_norm=((u*(acc_temp-min(acc_temp(:))))./(max(acc_temp(:))-min(acc_temp(:))))+l;   
    
    for i=1:nPop
        if TF(it) < 0.4
            for j=1:size(X,2)
                mrand=randi(nPop);
                Xnew(i,j)=X(i,j)+C1*rand*acc_norm(i,j).*(X(mrand,j)-X(i,j))*d;  
            end
        else
            for j=1:size(X,2)
                p=2*rand-C4;  
                T=C3*TF(it);
                if T>1
                    T=1;
                end
                if p<.5
                    Xnew(i,j)=Xbest(j)+C2*rand*acc_norm(i,j).*(T*Xbest(j)-X(i,j))*d;  
                else
                    Xnew(i,j)=Xbest(j)-C2*rand*acc_norm(i,j).*(T*Xbest(j)-X(i,j))*d;
                end
            end
        end
    end
    
    
    for i=1:nPop
        for j=1:nVar        
            if Xnew(i,j) < lb(j)
                    Xnew(i,j)= lb(j);
            elseif Xnew(i,j) > ub(j)
                    Xnew(i,j)= ub(j);
            end    
        end    
    end
    
    
    
    
    for i=1:nPop
        v=cost( Xnew(i,:));
        if v<Y(i)
            X(i,:)=Xnew(i,:);
            Y (i)=v;
        end
        
    end
    [var_Ybest,var_index] = min(Y);
    Convergence_curve(it)=var_Ybest;
    if var_Ybest<Scorebest
        Scorebest=var_Ybest;
        Score_index=var_index;
        Xbest = X(var_index,:);
        den_best=den(Score_index,:);
        vol_best=vol(Score_index,:);
        acc_best=acc_norm(Score_index,:);
    end
    
end


figure
semilogy(Convergence_curve,'b')
