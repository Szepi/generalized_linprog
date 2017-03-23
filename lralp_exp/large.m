clear all;

addpath ~/proj/MDPtoolbox
N=1000;
d=10;
k=2;
discount=0.9;
[P,R]=mdp_example_rand(N,d);
PR=mdp_computePR(P,R);
g=PR*100;

[V_LP,dum, policy, cpu_time,lambda] = mdp_LP(P, g, discount,eye(N),eye(N*d),ones(N,1));

P_pol=zeros(N,N);
g_pol=zeros(N,1);
for i=1:N
P_pol(i,:)=P(i,:,policy(i));
g_pol(i)=g(i,policy(i));
end;
V=inv(eye(N)-discount*P_pol)*g_pol;

Phi=round(rand(N,k));
Phi(:,1)=ones(N,1);

%for i=1:N
%	for j=1:k
%	Phi(i,j)=i^j;
%	end;
%end;


[V_ALP,r_ALP,policy_ALP,cpu_time,lambda_ALP]=mdp_LP(P,g,discount,Phi,eye(N*d),ones(N,1));


P_pol_alp=zeros(N,N);
g_pol_alp=zeros(N,1);
for i=1:N
P_pol_alp(i,:)=P(i,:,policy_ALP(i));
g_pol_alp(i)=g(i,policy_ALP(i));
end;

V_pol_alp=inv(eye(N)-discount*P_pol_alp)*g_pol_alp;

m=50;
W=round(rand(m,N*d));

[V_RLP,r_RLP,policy_RLP,cpu_time,lambda_RLP]=mdp_LP(P,g,discount,Phi,W,ones(N,1));


P_pol_rlp=zeros(N,N);
g_pol_rlp=zeros(N,1);
for i=1:N
P_pol_rlp(i,:)=P(i,:,policy_RLP(i));
g_pol_rlp(i)=g(i,policy_RLP(i));
end;

V_pol_rlp=inv(eye(N)-discount*P_pol_rlp)*g_pol_rlp;



[sum(V_LP) sum(V_pol_alp) sum(V_pol_rlp)]
