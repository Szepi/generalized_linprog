function [GV,r]=gamt(V,Phi,W);

d=length(V)/size(Phi,1);

k=1;
%for i=1:length(Phi)
%	for j=1:d
%		A(k,:)=Phi(i,:);
%		k=k+1;
%	end;
%end;
A=[];
for a=1:d
A=[A;Phi];
end;

for i=1:size(Phi,1)
e=zeros(size(Phi,1),1);
e(i)=1;
f=e'*Phi;

[r(:,i)]=glpk(f',-W*A,-W*V,-150*ones(size(Phi,2),1),150*ones(size(Phi,2),1), repmat(["U"],size(W*A,1),1),repmat(["C"],size((f))),1);
end;
GV=min((Phi*r)')';
