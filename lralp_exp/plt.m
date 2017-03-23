fptr=fopen('V','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-V(i));
end;
fclose(fptr);


fptr=fopen('V_pol_alp','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-V_pol_alp(i));
end;
fclose(fptr);

fptr=fopen('V_pol_rlp','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-V_pol_rlp(i));
end;
fclose(fptr);

fptr=fopen('V_alp','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-V_ALP(i));
end;
fclose(fptr);

fptr=fopen('V_rlp','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-V_RLP(i));
end;
fclose(fptr);

fptr=fopen('V_g','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-V_g(i));
end;
fclose(fptr);

fptr=fopen('V_t','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-V_t(i));
end;
fclose(fptr);


fptr=fopen('gj','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-gj(i));
end;
fclose(fptr);


fptr=fopen('ghj','w');
for i=1:N
fprintf(fptr,'%f\t%f\n',i,-ghj(i));
end;
fclose(fptr);


fptr=fopen('r_ALP','w');
fprintf(fptr,'%f\t%f\n',r_ALP(1),r_ALP(2));
fclose(fptr);


fptr=fopen('r_RLP','w');
fprintf(fptr,'%f\t%f\n',r_RLP(1),r_RLP(2));
fclose(fptr);

