function [TV]=qbell(P,g,discount,V);
k=1;
for a=1:size(P,3)
	for i=1:size(P,1)
		TV(k)=g(i,a)+discount*P(i,:,a)*V;
		k=k+1;
	end;
end;
TV=TV';
