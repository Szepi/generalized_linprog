clear all;
addpath ~/proj/MDPtoolbox

a=1e-1*2;
b=5;
N=10;
discount=0.98;
epsilon=0.1;
p=0.2;
q=[0.2 0.4];
d=length(q);
for i=2:N-1
	for j=1:d
	P(i,i-1,j)=q(j);
	P(i,i,j)=1-q(j)-p;
	P(i,i+1,j)=p;
	g(i,j)=-(i-1+60*q(j)^3);
	end;
end;

for j=1:d
g(N,j)=-(N-1+60*q(j)^3);
end;

for j=1:d
g(1,j)=-(0+60*q(j)^3);
end;

for i=1:d
P(1,1,i)=1-p;
P(1,2,i)=p;
P(N,N,i)=1-q(i);
P(N,N-1,i)=q(i);
end;
[V_PI, policy, iter] = mdp_policy_iteration(P, g, discount);
V_LP=V_PI;

P_pol=zeros(N,N);
g_pol=zeros(N,1);
for i=1:N
P_pol(i,:)=P(i,:,policy(i));
g_pol(i)=g(i,policy(i));
end;
V=inv(eye(N)-discount*P_pol)*g_pol;


for i=1:N
Phi(i,2)=i-1;
Phi(i,1)=1;
end;
[V_ALP,r_ALP,policy_ALP]=mdp_LP(P,g,discount,Phi,eye(N*d),ones(N,1));

m=5;
W=zeros(m,N*d);
fact=N/m;
for i=1:m
W(i,(i-1)*fact+1:i*fact)=1;
end;
for i=2:d
W(:,(i-1)*N+1:i*N)=W(:,1:N);
end;

%W=rand(size(W));

gj=gam(V_LP,Phi);
Hj=qbell(P,g,discount,V_LP);
ghj=gamt(Hj,Phi,W);


[V_RLP,r_RLP,policy_RLP]=mdp_LP(P,g,discount,Phi,W,ones(N,1));







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
g_pol_rlp(i)=g(i,policy_RLP(i));
end;
V_pol_rlp=inv(eye(N)-discount*P_pol_rlp)*g_pol_rlp;

[V V_pol_alp V_pol_rlp]

for i=1:N
P_pol=P(:,:,policy(i));
end;

d_pol=diag(P_pol^100);

sstates=lot(ones(N,1)/N,5);
W=zeros(m,N*d);
fact=N/m;
for i=1:m
W(i,sstates(i))=1;
end;
for i=2:d
W(:,(i-1)*N+1:i*N)=W(:,1:N);
end;

gj=gam(V_LP,Phi);
Hj=qbell(P,g,discount,V_LP);
ghj=gamt(Hj,Phi,W);

