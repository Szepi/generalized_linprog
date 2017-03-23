function [vec]=lot(A,num);
B=cumsum(A);
x=rand(1,num);
	for i=1:num
	vec(i)=min(find(B>x(i)));
	end;



