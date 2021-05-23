function [Nx,Ny]=shape_truss(xbar,L)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here




N(1)=1-xbar/L;
N(2)=xbar/L;

Nx=-1/L;
Ny=1/L;

end

