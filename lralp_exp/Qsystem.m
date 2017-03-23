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

