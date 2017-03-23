function [TV]=bell(P,g,discount,V);

for i=1:size(P,1)
	for a=1:size(P,3)
		Q(a)=g(i,a)+discount*P(i,:,a)*V;
	end;
TV(i)=max(Q);
end;
TV=TV';
