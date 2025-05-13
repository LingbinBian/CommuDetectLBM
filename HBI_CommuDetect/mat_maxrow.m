function [mat_max] = mat_maxrow(mat)
% maximum value of the row

[N,M]=size(mat);
mat_max=zeros(N,M);

[a b]=max(mat');

for i=1:N
   mat_max(i,b(i))=a(i);
end

end

