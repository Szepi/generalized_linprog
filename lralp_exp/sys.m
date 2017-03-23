clear all;

addpath ~/proj/MDPtoolbox
a=1e-1*2;
b=5;
N=10;
d=2;
discount=0.9;
epsilon=0.1;
p=1/2-epsilon/2;
q=[1/4 1/2+epsilon/2];

for i=2:N-1
	for j=1:2
	P(i,i-1,j)=q(j);
	P(i,i,j)=1-q(j)-p;
	P(i,i+1,j)=p;
	g(i,j)=-(a*i*i+b*j);
	end;
end;
g(N,1)=-(a*N*N+b*1);
g(N,2)=-(a*N*N+b*2);
g(1,1)=-(a*1*1+b*1);
g(1,2)=-(a*1*1+b*2);
P(1,1,1)=1-p;
P(1,2,1)=p;
P(1,1,2)=P(1,1,1);
P(1,2,2)=P(1,2,1);
P(N,N,1)=1-q(1);
P(N,N-1,1)=q(1);
P(N,N,2)=1-q(2);
P(N,N-1,2)=q(2);
g=g*1e-1;
[V_LP,dum, policy, cpu_time,lambda] = mdp_LP(P, g, discount,eye(N),eye(N*d),ones(N,1));

for i=1:N
Phi(i,2)=i*i*0.1;
Phi(i,1)=1;
end;
[V_ALP,r_ALP,policy_ALP,cpu_time,lambda_ALP]=mdp_LP(P,g,discount,Phi,eye(N*d),ones(N,1));
%W(1,:)=[1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0];
%W(2,:)=[0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 1 0];
%W(3,:)=[0 1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0];
%W(4,:)=[0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 1];

%W(1,:)=[1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%W(2,:)=[0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0];
%W(3,:)=[0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0];
%W(4,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1];



W(1,:)=[1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0];
W(2,:)=[0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0];
W(3,:)=[0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0];
W(4,:)=[0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0];
W(5,:)=[0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1];

%W(1,:)=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%W(2,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];


[V_RLP,r_RLP,policy_RLP,cpu_time,lambda_RLP]=mdp_LP(P,g,discount,Phi,W,ones(N,1));
for i=1:N
e=zeros(N,1);
e(i)=1;
[V_A,r_A(:,i),policy_ALP,cpu_time,lambda_A]=mdp_LP(P,g,discount,Phi,eye(N*d),e);
end;
for i=1:N
e=zeros(N,1);
e(i)=1;
[V_R,r_R(:,i),policy_R,cpu_time,lambda_R]=mdp_LP(P,g,discount,Phi,W,e);
end;


iters=100;
V=zeros(N,1);
for i=1:iters
V=bell(P,g,discount,V);
end;
iters=100;
V_g=min((Phi*r_A)')';
for i=1:iters
TV=bell(P,g,discount,V_g);
V_g=gam(TV,Phi);
end;

iters=100;
V_t=min((Phi*r_R)')';
V_t=zeros(N,1);
for i=1:iters
TV_t=qbell(P,g,discount,V_t);
[V_t,r_t]=gamt(TV_t,Phi,W);
end;

jb=gam(V_LP,Phi)
Tjb=bell(P,g,discount,jb);
gjb=gam(Tjb,Phi);


jb=gam(V_LP,Phi)
Tjb=qbell(P,g,discount,jb);
gtjb=gamt(Tjb,Phi,W);

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



plt;

