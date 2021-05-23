% Author Yixiao Zhang
function [Ke,ke,Qe,L]=Ke_truss_3D_finite(E,A,xn,ien,nen,ndf,nsd)


ke=zeros(2,2);

v=xn(:,ien(2))-xn(:,ien(1));% calculate the coordinates of elmemnts
L=norm(v);
[B] = B_truss(L);
Qe=zeros(nen,ndf*nen); %construct 2x4 matrix
v=v/L; % calculate the basic matrix in accordance witn β
if (nsd ==2)          % 2D case
    Qe(1,1)=v(1);
    Qe(1,2)=v(2);
    Qe(2,3)=v(1);
    Qe(2,4)=v(2);
    
    Qee(1,1)=v(1);      Qee(1,2)=v(2); 
    Qee(2,1)=-v(2);     Qee(2,2)=v(1); 
    Qee(3,3)=v(1);      Qee(3,4)=v(2); 
    Qee(4,3)=-v(2);      Qee(4,4)=v(1);  
    
    
elseif (nsd == 3)
    Qe(1,1)=v(1);
    Qe(1,2)=v(2);
    Qe(1,3)=v(3);
    Qe(2,4)=v(1);
    Qe(2,5)=v(2);
    Qe(2,6)=v(3);
end
ke(:,:)=E*A*(B)'*(B)*L;


Ke(:,:)=Qe(:,:)'*ke(:,:)*Qe(:,:);
end





