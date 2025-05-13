function [R] = assignment_matrix(Z_g,K)
% Sampling the assignment probability matrix from the posterior 
%
% Version 1.0
% 10-March-2021
% Copyright (c) 2021, Lingbin Bian

[N,S]=size(Z_g);
R=zeros(N,K);

for i=1:N
   R(i,:)=assignment_vector(Z_g(i,:),S);
end


% Nested function -----------------------------------------------------–---
function [r] = assignment_vector(z_g,S)
% z_g is the row vector: the labels of the subjects in the group for a
% specific node

% r: a row vector of assignment probability matrix for a specific subject

alpha=ones(1,K);   % hyper parameters of Dirichlet prior

S_k=zeros(1,K);
 
for s=1:S
    S_k(z_g(s))=S_k(z_g(s))+1;
end 

alpha_p=S_k+alpha;

n=1;    % one trial of multinomial
r=drchrnd(alpha_p,n);   % Dirichrit prior      

% Dirichlet random variable
    function [ r ] = drchrnd( a,n )
        % Dirichlet random variable
        p=length(a);
        r=gamrnd(repmat(a,n,1),1,n,p);
        r=r./repmat(sum(r,2),1,p);
    end
end

end