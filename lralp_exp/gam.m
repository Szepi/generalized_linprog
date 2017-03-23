function [GV]=gam(V,Phi);

for i=1:length(V)
e=zeros(length(V),1);
e(i)=1;
f=e'*Phi;
[r(:,i)]=glpk(f',-Phi,-V,-150*ones(size(Phi,2),1),150*ones(size(Phi,2),1), repmat(["U"],size(Phi,1),1),repmat(["C"],size(Phi,2),1),1);
end;

GV=min((Phi*r)')';
