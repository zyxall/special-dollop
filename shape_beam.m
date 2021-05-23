function [N,Nx,NXi,DetJ]=shape_beam(Xi,xn,ien,nen)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


for a=1:nen
    x(a)=xn(1,ien(a));
end


[N,Nxi]=shape_2beam(Xi,nen);
[J,DetJ]=jacobian_1D(Nxi,x);



    
    
end


