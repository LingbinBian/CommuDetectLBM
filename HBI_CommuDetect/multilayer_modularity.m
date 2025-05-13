function [S,Q] = multilayer_modularity(A,gamma,omega)

N=length(A{1,1});
T=length(A);
[B,twom] = multicat(A,gamma,omega);
[S,Q]= iterated_genlouvain(B,'moverandw');
S=postprocess_categorical_multilayer(S);

Q=Q/twom;
S=reshape(S,N,T);

% original_B=full(B);

end