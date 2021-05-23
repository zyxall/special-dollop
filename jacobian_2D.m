function [J,DetJ]=jacobian_2D(NXi,NEta,x,y,nen)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

J=zeros(2:2);

for a=1:nen
    J(1,1)=J(1,1)+NXi(a)*x(a);
    J(1,2)=J(1,2)+NXi(a)*y(a);
    J(2,1)=J(2,1)+NEta(a)*x(a);
    J(2,2)=J(2,2)+NEta(a)*y(a);
    
end

DetJ=det(J);







