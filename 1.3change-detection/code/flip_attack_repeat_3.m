%% gaussian shift-in-mean change detection
% "reverse" attack 
%% generate part
miu0=0; miu1=0.5;
change=1000;
total=10000;
K=3; B=2; %first 2 sensors are safe
X=zeros(K,total); %observation
L=zeros(K,total); %likelihood ratio
Z=zeros(K,total); %part SUM 1-t
W=zeros(K,total); %cusum with (largest positive drift)

times=10;  %100 times 65 seconds
h=10:50:1000; %thresholds
T=zeros(3,length(h),times);


for h_index=1:times
    
    for j=1:B  %normally affected
        X(j,1:change-1)=normrnd(miu0*ones(1,change-1),ones(1,change-1));
        X(j,change:total)=normrnd(miu1*ones(1,total-change+1),ones(1,total-change+1));
    end
    j=B+1;
    X(j,1:change-1)=normrnd(miu1*ones(1,change-1),ones(1,change-1));
    X(j,change:total)=normrnd(miu0*ones(1,total-change+1),ones(1,total-change+1));  %reversed
    
    L=X*(miu1-miu0)+(miu0^2-miu1^2)/2;
    Z=cumsum(L,2); %cum sum in row
    
    %calculate cusum
    for j=1:K
        if  W(j,1)<0
            W(j,1)=0;
        else
            W(j,1)=L(j,1);
        end
        for t=2:total
            W(j,t)=W(j,t-1)+L(j,t);
            if  W(j,t)<0
                W(j,t)=0;
            end
        end
    end
    
    %detection part
  
    % Lth alarm L=1
    % voting rule L=1 they are the same: the first alarm
    for i=1:length(h)
        for t=1:total
            if max(W(:,t))>h(i)
                T(1,i,h_index)=t; %alarm time
                break;
            end
        end
    end
    % Low-sum h*L/N(L==N so h is still the same value sequence)
    for i=1:length(h)
        for t=1:total
            ws=sort(W(:,t),1);
            if sum(ws(1:2))>h(i)
                T(2,i,h_index)=t; %alarm time
                break;
            end
        end
    end
    % Top-sum
    for i=1:length(h)
        for t=1:total
            ws=sort(W(:,t),1);
            if sum(ws(2:3))>h(i)
                T(3,i,h_index)=t; %alarm time
                break;
            end
        end
    end

end

% calculate expectations
E=zeros(3,length(h)); %expectaions of delay
for i=1:3
    for h_index=1:length(h)
        E(i,h_index)=sum(T(i,h_index,:)-change)/times;
    end
end
%% plot delay in different alarm strategy
plot(h,E(1,:),'r'); hold on;
plot(h,E(2,:),'k'); hold on;
plot(h,E(3,:),'b');

%% 
% for j=1:B
%     plot(W(j,:),'b'); hold on;
% end
% plot(W(B+1,:),'r'); hold on;
% plot(W(B+2,:),'y'); hold on;
% plot(W(B+3,:),'k'); hold on;
% for j=B+1:K
%    plot(W(j,:),'r'); hold on;
% end
%% 
% plot(W1,'g'); hold on;