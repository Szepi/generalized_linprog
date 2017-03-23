clear all;
addpath ~/proj/MDPtoolbox
a=1e-1*2;
b=5;
N=10;
discount=0.98;
epsilon=0.1;
p=0.2;
q=[0.2 0.4];
d=length(q);
for i=2:N-1
	for j=1:d
	P(i,i-1,j)=q(j);
	P(i,i,j)=1-q(j)-p;
	P(i,i+1,j)=p;
	g(i,j)=-((i-1)+60*q(j)^3);
	end;
end;

for j=1:d
g(N,j)=-(N-1+60*q(j)^3);
end;

for j=1:d
g(1,j)=-(0+60*q(j)^3);
end;

for i=1:d
P(1,1,i)=1-p;
P(1,2,i)=p;
P(N,N,i)=1-q(i);
P(N,N-1,i)=q(i);
end;

g=g/(1-discount);

fptr=fopen('queue.tex','w');
fprintf(fptr,'\\begin{tikzpicture}\n');
fprintf(fptr,'\\begin{axis}\n')

x=1;
	c='red';
	for j=1:d
		fprintf(fptr,'\\addplot[domain=-100:100,smooth,variable=\\y,%s]  plot ({',c);
		fprintf(fptr,'%f+%f*\\y',g(x,j),discount/(1-discount)*((1-p)*fun(x)+p*fun(x+1))-fun(x)/(1-discount));
		fprintf(fptr,'},{\\y});\n');
		c='black';
	end;


for x=2:N-1
	c='red';
	for j=1:d
		fprintf(fptr,'\\addplot[domain=-60:60,smooth,variable=\\y,%s]  plot ({',c);
		fprintf(fptr,'%f+%f*\\y',g(x,j),discount/(1-discount)*(q(j)*fun(x-1)+(1-q(j)-p)*fun(x)+p*fun(x+1)) -fun(x)/(1-discount));
		fprintf(fptr,'},{\\y});\n');
		c='black';
	end;

end;


x=N;
	c='red';
	for j=1:d
		fprintf(fptr,'\\addplot[domain=-60:60,smooth,variable=\\y,%s]  plot ({',c);
		fprintf(fptr,'%f+%f*\\y',g(x,j),discount/(1-discount)*(q(j)*fun(x-1)+(1-q(j))*fun(x))-fun(x)/(1-discount));
		fprintf(fptr,'},{\\y});\n');
		c='black';
	end;


fprintf(fptr,'\\end{axis}\n');
fprintf(fptr,'\\end{tikzpicture}\n');
fclose(fptr);



A(1,1)=g(1,1);
A(1,2)=discount/(1-discount)*((1-p)*fun(1)+p*fun(1+1)) -fun(1)/(1-discount);


j=1;
	
	for x=2:N-1
		A(x,1)=g(x,j);
		A(x,2)=discount/(1-discount)*(q(j)*fun(x-1)+(1-q(j)-p)*fun(x)+p*fun(x+1)) -fun(x)/(1-discount);
	end;

A(N,1)=g(N,1);
A(N,2)=discount/(1-discount)*(q(1)*fun(N-1)+(1-q(1))*fun(N)) -fun(N)/(1-discount);

A(N+1,1)=g(1,2);
A(N+1,2)=discount/(1-discount)*((1-p)*fun(1)+p*fun(1+1)) -fun(1)/(1-discount);

j=2;

	for x=2:N-1
		A(x+N,1)=g(x,j);
		A(x+N,2)=discount/(1-discount)*(q(j)*fun(x-1)+(1-q(j)-p)*fun(x)+p*fun(x+1)) -fun(x)/(1-discount);
	end;

A(N+N,1)=g(N,2);
A(N+N,2)=discount/(1-discount)*(q(2)*fun(N-1)+(1-q(2))*fun(N)) -fun(N)/(1-discount);



m=5;
W=zeros(m,N*d);
fact=N/m;
for i=1:m
W(i,(i-1)*fact+1:i*fact)=1;
end;
for i=2:d
W(:,(i-1)*N+1:i*N)=W(:,1:N);
end;
W=W/5;
B=W*A;



fptr=fopen('queuered.tex','w');
fprintf(fptr,'\\begin{tikzpicture}\n');
fprintf(fptr,'\\begin{axis}\n')
for i=1:length(B)
		fprintf(fptr,'\\addplot[domain=-60:60,smooth,variable=\\y,black]  plot ({');
		fprintf(fptr,'%f+%f*\\y',B(i,1),B(i,2));
		fprintf(fptr,'},{\\y});\n');

end;
fprintf(fptr,'\\end{axis}\n');
fprintf(fptr,'\\end{tikzpicture}\n');
fclose(fptr);

