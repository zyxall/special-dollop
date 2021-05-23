%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_eq                                                %
% author: Nima Rahbar                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id,neq]=number_eq(idb,nnp,ndf)
%------------------------------------------------------------------------
% Purpose:
% number equations in global matrix
% id(i,N): equation number associated with dof i of node N
%
% Synopsis:
% [id,neq]=number_eq(idb,nnp,ndf)
%
% Variable Description:
% nnp - number of nodes
% ndf - number of equations per node
%------------------------------------------------------------------------

id = zeros(ndf,nnp); % initialize id
neq = 0; % number of equations
for N = 1:nnp
    for i = 1:ndf
        if (idb(i,N) == 0) % assign an equation number to all non prescribed nodes
            neq = neq+1;
            id(i,N) = neq;
        end;
    end;
end;