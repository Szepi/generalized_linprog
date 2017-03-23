function [V,r, policy] = mdp_LP(P, R, discount,Phi,W,f)
% mdp_LP    Resolution of discounted MDP with linear programming
% Arguments --------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA) = transition matrix 
%              P could be an array with 3 dimensions or 
%              a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%              a cell array (1xA), each cell containing a sparse matrix (SxS) or
%              a 2D array(SxA) possibly sparse  
%   discount = discount rate, in ]0; 1[
% Evaluation -------------------------------------------------------------
%   V(S)   = optimal values
%   policy(S) = optimal policy
%   cpu_time = used CPU time

% MDPtoolbox: Markov Decision Processes Toolbox
% Copyright (C) 2009  INRA
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
%    * Redistributions of source code must retain the above copyright notice, 
%      this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice, 
%      this list of conditions and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the <ORGANIZATION> nor the names of its contributors 
%      may be used to endorse or promote products derived from this software 
%      without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
% OF THE POSSIBILITY OF SUCH DAMAGE.


cpu_time = cputime;


% check of arguments
if (discount <= 0 || discount >= 1)
    disp('--------------------------------------------------------')
    disp('MDP Toolbox ERROR: Discount rate must be in ]0; 1[')
    disp('--------------------------------------------------------')
else
    if iscell(P)
        S = size(P{1},1);
        A = length(P);
    else
        S = size(P,1);
        A = size(P,3);
    end

    PR = mdp_computePR(P,R);
    
    
    % The objective is to resolve : min V / V >= PR + discount*P*V
    % The function linprog of the optimisation Toolbox of Mathworks resolves :
    % min f'* x / M * x <= b
    % So the objective could be expressed as : min V / (discount*P-I) * V <= - PR
    % To avoid loop on states, the matrix M is structured following actions M(A*S,S)
    
    %f=ones(S,1);
    M = [];
    PR_new=[];	
    if iscell(P)
        for a=1:A; M=[M;discount*P{a}-speye(S)]; end
    else
        for a=1:A; M=[M;discount*P(:,:,a)-speye(S)]; PR_new=[PR_new;PR(:,a)];end

    end


%k=1;
%for i=1:size(P,1)
%PR_new(k)=R(i,1);
%k=k+1;
%PR_new(k)=R(i,2);
%k=k+1;
%end;
%PR_new=PR_new';

%[M ,PR_new]
B=1e4;
%PARAM.lpsolver=1;
%PARAM.msglev=3;
[r,FMIN,ERRNUM,EXTRA] = glpk((f'*Phi)',W*M*Phi,-W*PR_new,-B*ones(size(Phi,2),1),B*ones(size(Phi,2),1), repmat(["U"],size(W*M*Phi,1),1),repmat(["C"],size((f'*Phi))),1);
%FMIN
%ERRNUM
%EXTRA
cost=(f'*Phi)';
Amat=W*M*Phi;
bmat=-W*PR_new;
%[r,dum1,dum2,dum3,lambda] = linprog((f'*Phi)',M*Phi,-PR_new);
[V, policy] =  mdp_bellman_operator(P,PR,discount,Phi*r);
V=Phi*r;
end
cpu_time = cputime - cpu_time;
