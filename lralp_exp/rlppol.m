function [policy_RLP]=rlppol(P,g,discount,Phi,c)
N=size(P,1);
policy_RLP=ones(N,1);
q=[0.2 0.4 0.6 0.8];
A=4;
Q=zeros(A,1);
for s=2:N-1
m=5;
[W,c_local]=generateW(m,N,A,discount,s-1);
[V_RLPp,r_RLP,policy_RLP]=mdp_LP(P,g,discount,Phi,W,c_local);
[W,c_local]=generateW(m,N,A,discount,s+1);
[V_RLPn,r_RLP,policy_RLP]=mdp_LP(P,g,discount,Phi,W,c_local);
	for a=1:A
		Q(a)=g(s,a)+discount*P(s,s-1,a)*V_RLPp(s-1)+discount*P(s,s+1,a)*V_RLPn(s+1);
	end;
	policy_RLP(s)=min(find(Q==max(Q)));
end;
