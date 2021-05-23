function [point,w] = gauss(Nintx,Ninty,nsd)

if nsd==1
    [point,w]=gauss_1d(Nintx);
elseif nsd==2
    [point,w]=gauss_2d(Nintx,Ninty);

end

