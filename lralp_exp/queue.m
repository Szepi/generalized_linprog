clear all;
close all;
addpath ~/proj/MDPtoolbox
a=1e-1*2;
b=5;
N=1000;
discount=1-1/N;
epsilon=0.1;
p=0.2;
q=[0.2 0.4 0.6 0.8];
d=length(q);

for i=2:N-1
	for j=1:d
	P(i,i-1,j)=q(j);
	P(i,i,j)=1-q(j)-p;
	P(i,i+1,j)=p;
	g(i,j)=-((i/N)^2+q(j)^3);
	end;
end;

for j=1:d
g(N,j)=-((N/N)^2+(q(j))^3);
end;

for j=1:d
g(1,j)=-((1/N)+(q(j))^3);
end;

g=g;

for i=1:d
P(1,1,i)=1-p;
P(1,2,i)=p;
P(N,N,i)=1-q(i);
P(N,N-1,i)=q(i);
end;


%[V_LP,dum, policy, cpu_time,lambda] = mdp_LP(P, g, discount,eye(N),eye(N*d),ones(N,1));
%[V_LP,dum, policy,] = mdp_LP(P, g, discount,eye(N),eye(N*d),ones(N,1));
[V_PI, policy, iter, cpu_time] = mdp_policy_iteration(P, g, discount);
V_LP=V_PI;

P_pol=zeros(N,N);
g_pol=zeros(N,1);
for i=1:N
P_pol(i,:)=P(i,:,policy(i));
g_pol(i)=g(i,policy(i));
end;
V=inv(eye(N)-discount*P_pol)*g_pol;




c=ones(N,1);
c1=zeros(N,1);
zeta=discount;
for i=1:N
Phi(i,2)=i;
Phi(i,1)=1;
Phi(i,3)=i*i;
Phi(i,4)=i^3;
c(i)=(1-zeta)*zeta^(i);
end;
c=c/sum(c);
[V_ALP,r_ALP,policy_ALP]=mdp_LP(P,g,discount,Phi,eye(N*d),c);


m=5;

fact=N*d/m;
W=zeros(m,N*d);
fact=N/m;

for i=1:m
%W(i,(i-1)*fact+1:i*fact)=1;
ind=(i-1)*fact+1;
ind=lot(c,1);
W(i,ind)=1;
c1(ind)=1;
end;
%W(i,N)=1;
%c1(N)=1;
for i=2:d
W(:,(i-1)*N+1:i*N)=W(:,1:N);
end;



%W=rand(size(W));


[V_RLP,r_RLP,policy_RLP]=mdp_LP(P,g,discount,Phi,W,c);


P_pol_alp=zeros(N,N);
g_pol_alp=zeros(N,1);
for i=1:N
P_pol_alp(i,:)=P(i,:,policy_ALP(i));
g_pol_alp(i)=g(i,policy_ALP(i));
end;

V_pol_alp=inv(eye(N)-discount*P_pol_alp)*g_pol_alp;




P_pol_rlp=zeros(N,N);
g_pol_rlp=zeros(N,1);
for i=1:N
P_pol_rlp(i,:)=P(i,:,policy_RLP(i));
action=randi(4);
%P_pol_rlp(i,:)=P(i,:,action);
g_pol_rlp(i)=g(i,policy_RLP(i));
%g_pol_rlp(i)=g(i,action);
end;

V_pol_rlp=inv(eye(N)-discount*P_pol_rlp)*g_pol_rlp;


[sum(V_LP)  sum(V_pol_alp) sum(V_pol_rlp)]

figure(1)
plot(V_LP); hold on; plot(V_ALP,'r'); plot(V_RLP,'g');
legend('LP','ALP','RLP');
d=1000;

figure(2)
plot(V_LP(1:d)); hold on; plot(V_pol_alp(1:d),'r'); 
plot(V_pol_rlp(1:d),'g');
legend('LP-pol','ALP-pol','RLP-pol');

figure(3)
plot(policy(1:d)); hold on; plot(policy_ALP(1:d),'r'); 
plot(policy_RLP(1:d),'g');
legend('LP-pol','ALP-pol','RLP-pol');
