function [R_esti] = assign_esti(Z_g,K,N_simu)
% This function estimates the assignment probability matrix
%
% Version 1.0
% 11-March-2021
% Copyright (c) 2021, Lingbin Bian


[N,S_subj]=size(Z_g);

R_esti=zeros(N,K);
R_group=cell(N_simu,1);

for n=1:N_simu
    R_group{n,1}=assignment_matrix(Z_g,K);
end

R_mean=zeros(N_simu,1);

    for i=1:N
        for k=1:K
            for n=1:N_simu
                R_mean(n,1)=R_group{n,1}(i,k);
            end
            R_esti(i,k)=mean(R_mean);
        end
    end

end

