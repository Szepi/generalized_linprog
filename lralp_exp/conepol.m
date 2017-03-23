function [policy_RLP]=conepol(P,g,discount,Phi,c)
N=size(P,1);
policy_RLP=ones(N,1);
q=[0.2 0.4 0.6 0.8];
A=4;
Q=zeros(A,1);
for s=2:N-1
m=5;
W=zeros(m+1,N*A);
	fact=N/m;

	for i=1:m-1
	ind=(i-1)*fact+1;
	W(i,ind)=1;
	end;
  W(m,N)=1;
  W(m+1,s-1)=1;

	for i=2:A
	W(:,(i-1)*N+1:i*N)=W(:,1:N);
	end;
c_local=zeros(N,1);
c_local(s-1)=1;
[V_RLPp,r_RLP,policy_RLP]=mdp_LP(P,g,discount,Phi,W,c_local);
  W(m+1,s-1)=0;
    W(m+1,s+1)=1;
    c_local(s-1)=0;
    c_local(s+1)=1;
[V_RLPn,r_RLP,policy_RLP]=mdp_LP(P,g,discount,Phi,W,c_local);
	for a=1:A
		Q(a)=g(s,a)+discount*P(s,s-1,a)*V_RLPp(s-1)+discount*P(s,s+1,a)*V_RLPn(s+1);
	end;
	policy_RLP(s)=min(find(Q==max(Q)));
end;
