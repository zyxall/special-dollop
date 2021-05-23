%Author Yixiao Zhang
function [Ke,ke,Qe]=Ke_beam3D_finite(E,A,Iy,Iz,G,J,psi,xn,ien,nen,ndf,nsd,nel,Nintx)
kte=zeros(2,2,nel);
ktt=zeros(2,2,nel);
kbe=zeros(4,4,nel);
kbe2=zeros(4,4,nel);

if nsd==2
    ke=zeros(6,6,nel);
elseif nsd==3
    ke=zeros(12,12,nel);
end

for e=1:nel

v=xn(:,ien(2,e))-xn(:,ien(1,e)); % two ponints' coordinate difference
L=norm(v); % Length of member

[point,w]=gauss(Nintx,0,0,nsd);
if nsd ==2 
    I(e)=Iz(e);
 

    for i=1:Nintx
        Xi=point(i,1);
        Wi(i,1)=w(i,1);
       [Bb,Bb2] = B_beam(Xi,L);   
       kbe(:,:,e)=kbe(:,:,e)+E(e)*I(e)*(Bb)'*Bb*Wi(i,1)*L/2;%Beam B matrix 
    end
 
[Bt]=B_truss(L);
kte(:,:,e)=E(e)*A(e)*(Bt)'*(Bt)*L; %truss B matrix 


ke(2:3,2:3,e)=kbe(1:2,1:2,e); %rearrange Ke 4x4 matrix into 6x6
ke(2:3,5:6,e)=kbe(1:2,3:4,e);
ke(5:6,2:3,e)=kbe(3:4,1:2,e);
ke(5:6,5:6,e)=kbe(3:4,3:4,e);

ke(1,1,e)=kte(1,1,e);
ke(1,4,e)=kte(1,2,e);
ke(4,1,e)=kte(2,1,e);
ke(4,4,e)=kte(2,2,e);

 %construct 2x4 matri
        v=v/L; % calculate the basic matrix in accordance witn β        
      
            % 2D case        
           Qe(1,1,e)=v(1);        
           Qe(1,2,e)=v(2); 
           Qe(2,1,e)=-v(2);
           Qe(2,2,e)=v(1);
           Qe(3,3,e)=1;
           Qe(4,4,e)=v(1);        
           Qe(4,5,e)=v(2);
           Qe(5,4,e)=-v(2);        
           Qe(5,5,e)=v(1);
           Qe(6,6,e)=1;
           
        
        elseif  (nsd ==3)

   for i=1:Nintx
        Xi=point(i,1);
        Wi(i,1)=w(i,1);
       [Bb,Bb2] = B_beam(Xi,L);   
       kbe(:,:,e)=kbe(:,:,e)+E(e)*Iz(e)*(Bb)'*Bb*Wi(i,1)*L/2;%Beam B matrix 
       kbe2(:,:,e)=kbe2(:,:,e)+E(e)*Iy(e)*(Bb2)'*Bb2*Wi(i,1)*L/2;
   end
    
[Bt]=B_truss(L);
kte(:,:,e)=E(e)*A(e)*(Bt)'*(Bt)*L; %truss B matrix 
ktt(:,:,e)=G*J*(Bt)'*(Bt)*L; %truss B matrix 

ke(1,1,e)=kte(1,1,e); ke(4,4,e)=ktt(1,1,e);
ke(1,7,e)=kte(1,2,e); ke(4,10,e)=ktt(1,2,e);
ke(7,1,e)=kte(2,1,e); ke(10,4,e)=ktt(2,1,e);
ke(7,7,e)=kte(2,2,e); ke(10,10,e)=ktt(2,2,e);

ke(2,:,e)=[  0 ,kbe(1,1,e),       0   ,  0  ,       0   ,kbe(1,2,e), 0  ,kbe(1,3,e),     0     ,  0  ,     0      ,kbe(1,4,e)];                              %rearrange Ke 4x4 matrix into 6x6
ke(3,:,e)=[  0 ,    0     ,kbe2(1,1,e),  0  ,kbe2(1,2,e),     0    , 0  ,    0     ,kbe2(1,3,e),  0  , kbe2(1,4,e),    0     ];
ke(5,:,e)=[  0 ,    0     ,kbe2(2,1,e),  0  ,kbe2(2,2,e),     0    , 0  ,    0     ,kbe2(2,3,e),  0  , kbe2(2,4,e),    0     ];
ke(6,:,e)=[  0 ,kbe(2,1,e),       0   ,  0  ,       0   ,kbe(2,2,e), 0  ,kbe(2,3,e),     0     ,  0  ,     0      ,kbe(2,4,e)];

ke(8,:,e)=[  0 ,kbe(3,1,e),       0   ,  0  ,       0   ,kbe(3,2,e), 0  ,kbe(3,3,e),     0     ,  0  ,     0      ,kbe(3,4,e)];                              %rearrange Ke 4x4 matrix into 6x6
ke(9,:,e)=[  0 ,    0     ,kbe2(3,1,e),  0  ,kbe2(3,2,e),     0    , 0  ,    0     ,kbe2(3,3,e),  0  , kbe2(3,4,e),    0     ];
ke(11,:,e)=[  0 ,    0     ,kbe2(4,1,e),  0  ,kbe2(4,2,e),     0    , 0  ,    0     ,kbe2(4,3,e),  0  , kbe2(4,4,e),    0     ];
ke(12,:,e)=[  0 ,kbe(4,1,e),       0   ,  0  ,       0   ,kbe(4,2,e), 0  ,kbe(4,3,e),     0     ,  0  ,     0      ,kbe(4,4,e)];




%         ke(:,:,e)=[(E(e)*A(e)/L),0          ,0          ,0    ,0          ,0          ,-E(e)*A(e)/L,0           ,0           ,0     ,0          ,0;       
%                   0     ,12*E(e)*Iz(e)/L^3,0          ,0    ,0          ,6*E(e)*Iz(e)/L^2 ,0     ,-12*E(e)*Iz(e)/L^3,0           ,0     ,0          ,6*E(e)*Iz(e)/L^2;
%                   0     ,0          ,12*E(e)*Iy(e)/L^3,0    ,-6*E(e)*Iy(e)/L^2,0          ,0     ,0           ,-12*E(e)*Iy(e)/L^3,0     ,-6*E(e)*Iy(e)/L^2,0;
%                   0     ,0          ,0          ,G*J/L,0          ,0          ,0     ,0           ,0           ,-G*J/L,0          ,0; 
%                   0     ,0          ,-6*E(e)*Iy(e)/L^2,0    ,4*E(e)*Iy(e)/L   ,0          ,0     ,0           ,6*E(e)*Iy(e)/L^2  ,0     ,2*E(e)*Iy(e)/L   ,0;
%                   0     ,6*E(e)*Iz(e)/L^2 ,0          ,0    ,0          ,4*E(e)*Iz(e)/L   ,0     ,-6*E(e)*Iz(e)/L^2 ,0           ,0     ,0          ,2*E(e)*Iz(e)/L;
%                 -(E(e)*A(e)/L),0          ,0          ,0    ,0          ,0          ,E(e)*A(e)/L ,0           ,0           ,0     ,0          ,0;       
%                   0     ,-12*E(e)*Iz(e)/L^3,0         ,0    ,0          ,-6*E(e)*Iz(e)/L^2,0     ,12*E(e)*Iz(e)/L^3 ,0           ,0     ,0          ,-6*E(e)*Iz(e)/L^2;
%                   0     ,0          ,-12*E(e)*Iy(e)/L^3,0   ,6*E(e)*Iy(e)/L^2 ,0          ,0     ,0           ,12*E(e)*Iy(e)/L^3 ,0     ,6*E(e)*Iy(e)/L^2 ,0;
%                   0     ,0          ,0          ,-G*J/L,0         ,0          ,0     ,0           ,0           ,G*J/L,0           ,0; 
%                   0     ,0          ,-6*E(e)*Iy(e)/L^2,0    ,2*E(e)*Iy(e)/L   ,0          ,0     ,0           ,6*E(e)*Iy(e)/L^2  ,0     ,4*E(e)*Iy(e)/L   ,0;
%                   0     ,6*E(e)*Iz(e)/L^2 ,0          ,0    ,0          ,2*E(e)*Iz(e)/L   ,0     ,-6*E(e)*Iz(e)/L^2 ,0           ,0     ,0          ,4*E(e)*Iz(e)/L];

 %construct 12x12 matrix
v=xn(:,ien(2,e))-xn(:,ien(1,e));% calculate the coordinates of elmemnts
v=v/norm(v); %r

    if v(1)==0 && v(3) ==0
            if v(2)>0
                r(:,:,e) = [0 1 0; -1 0 0; 0 0 1];
            else
                r(:,:,e) = [0 -1 0; 1 0 0; 0 0 1];
            end
    else
        % 3D case
        rx(1,:,e)=[v(1),v(2),v(3)];
        rx(2,:,e)=[-v(1)*v(2)/sqrt(v(1)^2+v(3)^2),sqrt(v(1)^2+v(3)^2),-v(2)*v(3)/sqrt(v(1)^2+v(3)^2)];
        rx(3,:,e)=[-v(3)/sqrt(v(1)^2+v(3)^2),0,v(1)/sqrt(v(1)^2+v(3)^2)];
        rpsi(1,:,e)=[1,0,0];
        rpsi(2,:,e)=[0,cos(psi),sin(psi)];
        rpsi(3,:,e)=[0,-sin(psi),cos(psi)];
        r(:,:,e)=rx(:,:,e)*rpsi(:,:,e);
    end
    Qe(1:3,1:3,e)=r(:,:,e);
    Qe(4:6,4:6,e)=r(:,:,e);
    Qe(7:9,7:9,e)=r(:,:,e);
    Qe(10:12,10:12,e)=r(:,:,e);
    
    
end

% The big 12x12 Ke matrix
Ke(:,:,e)=Qe(:,:,e)'*ke(:,:,e)*Qe(:,:,e);

end
end





