%gaussian shift-in-mean change detection

miu0=0; miu1=0.2;
change=50;
total=200;
K=10; B=7; %first 7 sensors are affected
X=zeros(K,total); %observation
L=zeros(K,total); %likelihood ratio
Z=zeros(K,total); %part SUM 1-t
W=zeros(K,total); %cusum with (largest positive drift)

for j=1:B  %affected
    X(j,1:change-1)=normrnd(miu0*ones(1,change-1),ones(1,change-1));
    X(j,change:total)=normrnd(miu1*ones(1,total-change+1),ones(1,total-change+1));
end
for j=B+1:K
    X(j,1:total)=normrnd(miu0*ones(1,total),ones(1,total));  %unaffected
end

L=X*(miu1-miu0)+(miu0^2-miu1^2)/2;
Z=cumsum(L,2); %cum sum in row
Z1=sum(Z(1:B,:),1);  %SUM-CUSUM (all affected)
Z2=sum(Z(1:K,:),1);  %SUM-CUSUM (partly affected)


for t=1:total
    W(:,t)=Z(:,t)-min(Z(:,1:t),[],2);
    W1(t)=Z1(t)-min(Z1(1:t),[],2);
    W2(t)=Z2(t)-min(Z2(1:t),[],2);
end


%% 
for j=1:B
    plot(W(j,:),'b'); hold on;
end
for j=B+1:K
   % plot(W(j,:),'r'); hold on;
end
%% 
plot(W1,'g'); hold on;