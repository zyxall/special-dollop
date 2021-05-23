function [J,DetJ]=jacobian_1D(Nxi,x)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

J=zeros(4:4);


    J(1,1:4)=[Nxi(1)*x(1),Nxi(1)*x(2),Nxi(1)*x(1),Nxi(1)*x(2)];
    J(1,1:4)=[Nxi(2)*x(1),Nxi(2)*x(2),Nxi(2)*x(1),Nxi(2)*x(2)];
    J(1,1:4)=[Nxi(3)*x(1),Nxi(3)*x(2),Nxi(3)*x(1),Nxi(3)*x(2)];
    J(1,1:4)=[Nxi(4)*x(1),Nxi(4)*x(2),Nxi(4)*x(1),Nxi(4)*x(2)];
    

DetJ=det(J);







