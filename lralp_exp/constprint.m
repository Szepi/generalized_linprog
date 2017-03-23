a=1e-1*2;
b=5;
N=10;
discount=0.9;
epsilon=0.1;
p=1/2-epsilon/2;
q=[1/4 1/2+epsilon/2];


for i=2:N-1
	for j=1:2
	P(i,i-1,j)=q(j);
	P(i,i,j)=1-q(j)-p;
	P(i,i+1,j)=p;
	g(i,j)=-(a*i*i+b*j);
	end;
end;
g(N,1)=-(a*N*N+b*1);
g(N,2)=-(a*N*N+b*2);
g(1,1)=-(a*1*1+b*1);
g(1,2)=-(a*1*1+b*2);
g=g*1e-1/(1-discount);

fptr=fopen('constraints.tex','w');
fprintf(fptr,'\\begin{tikzpicture}\n');
fprintf(fptr,'\\begin{axis}[xmax=10,xmin=-20,ymax=5,ymin=-10]\n')
%fprintf(fptr,'\\begin{axis}\n')
%fprintf(fptr,'\\draw[->] (-3,0) -- (4.2,0) node[right] {$x$};\n');
%fprintf(fptr,'\\draw[->] (0,-3) -- (0,4.2) node[above] {$y$};\n');

for x=2:N-1
	c='red';
	for j=1:2
		fprintf(fptr,'\\addplot[scale=1.0,domain=-5:5,smooth,variable=\\y,%s]  plot ({',c);
%		fprintf(fptr,'%f+%f*(\\y-1)+%f*\\y+%f*(\\y+1)-%f*\\y',(a*x+b*q(j)),discount*q(j),discount*(1-q(j)-p),discount*p,x);
		fprintf(fptr,'%f+%f*\\y',g(x,j),discount/(1-discount)*(q(j)*fun(x-1)+(1-q(j)-p)*fun(x)+p*fun(x+1)) -fun(x)/(1-discount));
		fprintf(fptr,'},{\\y});\n');
		c='black';
	end;

end;




fprintf(fptr,'\\addplot[scale=1.0,domain=-5:5,smooth,variable=\\y,red]  plot ({',c);
fprintf(fptr,'%f+%f*\\y',g(1,1) ,discount/(1-discount)*((1-p)*fun(1)+p*fun(1+1)) -fun(1)/(1-discount));
fprintf(fptr,'},{\\y});\n');


fprintf(fptr,'\\addplot[scale=1.0,domain=-5:5,smooth,variable=\\y,black]  plot ({',c);
fprintf(fptr,'%f+%f*\\y',g(1,2) ,discount/(1-discount)*((1-p)*fun(1)+p*fun(1+1)) -fun(1)/(1-discount));
fprintf(fptr,'},{\\y});\n');

fprintf(fptr,'\\addplot[scale=1.0,domain=-5:5,smooth,variable=\\y,red]  plot ({',c);
fprintf(fptr,'%f+%f*\\y',g(N,1) ,discount/(1-discount)*(q(1)*fun(N-1)+(1-q(1))*fun(N)) -fun(N)/(1-discount));
fprintf(fptr,'},{\\y});\n');


fprintf(fptr,'\\addplot[scale=1.0,domain=-5:5,smooth,variable=\\y,black]  plot ({',c);
fprintf(fptr,'%f+%f*\\y',g(N,2) ,discount/(1-discount)*(q(2)*fun(N-1)+(1-q(2))*fun(N)) -fun(N)/(1-discount));
fprintf(fptr,'},{\\y});\n');


fprintf(fptr,'\\end{axis}\n');

fprintf(fptr,'\\end{tikzpicture}\n');
fclose(fptr);
%({2*\y+3},{\y})



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

B=W*A;



fptr=fopen('reduced.tex','w');
fprintf(fptr,'\\begin{tikzpicture}\n');
fprintf(fptr,'\\begin{axis}[xmax=10,xmin=-20,ymax=5,ymin=-10]\n')
for i=1:length(B)
		fprintf(fptr,'\\addplot[scale=1.0,domain=-5:5,smooth,variable=\\y,black]  plot ({');
		fprintf(fptr,'%f+%f*\\y',B(i,1),B(i,2));
		fprintf(fptr,'},{\\y});\n');

end;
fprintf(fptr,'\\end{axis}\n');

fprintf(fptr,'\\end{tikzpicture}\n');

fclose(fptr);
