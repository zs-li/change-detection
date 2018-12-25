%gausian shift-in-mean change detection

miu0=0; miu1=0.5;
change=1000;
total=10000;
X=zeros(1,total); %observation
L=zeros(1,total); %likelihood ratio
Z=zeros(1,total); %part SUM 1-t
W1=zeros(1,total); %cusum with (stay above zero)
W0=zeros(1,total); %cusum with (largest positive drift)

X(1:change-1)=normrnd(miu0*ones(1,change-1),ones(1,change-1));
X(change:total)=normrnd(miu1*ones(1,total-change+1),ones(1,total-change+1));

L=X*(miu1-miu0)+(miu0^2-miu1^2)/2;
Z=cumsum(L);

% 
% for t=1:total
%     W0(t)=Z(t)-min(Z(1:t));
%    
% end

W1(1)=L(1);
if  W1(1)<0
    W1(1)=0;
end
for t=2:total
    W1(t)=W1(t-1)+L(t);
    if  W1(t)<0
        W1(t)=0;
    end
end
