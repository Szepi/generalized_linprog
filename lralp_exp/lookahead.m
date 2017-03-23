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





policy_CS=rlppol(P,g,discount,Phi,c);

policy_cone=conepol(P,g,discount,Phi,c);


P_pol_CS=zeros(N,N);
g_pol_CS=zeros(N,1);
for i=1:N
P_pol_CS(i,:)=P(i,:,policy_CS(i));
g_pol_CS(i)=g(i,policy_CS(i));
end;

V_pol_CS=inv(eye(N)-discount*P_pol_CS)*g_pol_CS;




P_pol_cone=zeros(N,N);
g_pol_cone=zeros(N,1);
for i=1:N
P_pol_cone(i,:)=P(i,:,policy_cone(i));
action=randi(4);
%P_pol_rlp(i,:)=P(i,:,action);
g_pol_cone(i)=g(i,policy_cone(i));
%g_pol_rlp(i)=g(i,action);
end;

V_pol_cone=inv(eye(N)-discount*P_pol_cone)*g_pol_cone;


[sum(V_LP)  sum(V_pol_CS) sum(V_pol_cone)]


figure(1)
plot(V_LP); hold on; plot(V_pol_CS,'r'); 
plot(V_pol_cone,'g');
legend('LP-pol','CS-pol','cone-pol');

figure(2)
plot(policy); hold on; plot(policy_CS,'r'); 
plot(policy_cone,'g');
legend('LP-pol','CS-pol','cone-pol');
