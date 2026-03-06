function [S,Q] = multilayer_modularity(A,gamma,omega)
%   Input: A: Cell array of NxN adjacency matrices for each layer of an
%          unordered multilayer undirected network
%          gamma: intralayer resolution parameter
%          omega: interlayer coupling strength
%
%   Output: S: N_node x N_subj matrix of community labels
N=length(A{1,1});
T=length(A);
[B,twom] = multicat(A,gamma,omega);
%[S,Q]= iterated_genlouvain(B,'moverandw');
[S,Q]= iterated_genlouvain(B,'move');
S=postprocess_categorical_multilayer(S);

Q=Q/twom;
S=reshape(S,N,T);

% original_B=full(B);

end