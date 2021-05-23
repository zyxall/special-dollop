%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vibration                                                %
% author: Yixiao Zhang                                     %
% Date:11/27/2020                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; % removes all variables from the workspace.
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%
% Mesh %
%%%%%%%%

nsd=3; % number of space dimension
ndf=6; % number of degree of freedom per node
nen=2;

xinc=1; % per length in x direction
% de=input('Degree from (0,60] (ex:60):  ')
yinc=xinc*sind(60); % per length in y direction
zinc=1; %per length in z direction
n=8; %number of groups (not less than 1 story)
% le=input('Length of cell consisted by integer under [3,+∞) (ex:5):  ')
nxd=n; %nubmer in one row of x direction****fixed*****
nyd=3; %nubmer in one row of y direction****fixed*****
nzd=2; %nubmer in one row of z direction****fixed*****
% hi=input('height of cell consisted by integer under [1,+∞] (ex: 2): ')
% 
nnn=n;

% change it to square model
%%%%%%%%%%%%%
% Geometric %
%%%%%%%%%%%%%


%%%%%%%%%%%%%
% Vibration %
%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%
% Nodal coordinates %
%%%%%%%%%%%%%%%%%%%%%
% xn(i,N):= coordinate i for node N
% N=1,...,nnp
% i=1,...,nsd

pn=0;

for k=1:nzd
 for j=1:2:(nyd+(n-1)*2)
    for i=1:nxd
        pn=pn+1;
        xn(1,pn)=(i-1)*xinc;
        xn(2,pn)=(j-1)*yinc;
        xn(3,pn)=(k-1)*zinc;
    end
    
 end  
  
end

for k=1:nzd
for j=2:2:(nyd+(n-1)*2)
    for i=1:nxd+1
        pn=pn+1;
        xn(1,pn)=(i-1)*xinc-xinc/2;
        xn(2,pn)=(j-1)*yinc;
        xn(3,pn)=(k-1)*zinc;
    end
end 
end

nnp=size(xn,2);
Length=norm(xn(:,1)-xn(:,nxd*n+1)); % Length of member
%%%%%%%%%%%%%%%%
% Connectivity %
%%%%%%%%%%%%%%%%
% ien(a,e) = N
% N: global node number - N=1,...,nnp
% e: element number - e=1,...,nel
% a: local node number - a=1,...,nen

rn=nnz(1:2:(nyd+(n-1)*2)); %amount of odd row

E1=1;
E2=2;
rn=nnz(1:2:(nyd+(n-1)*2)); %amount of odd row

    e=0;
    
    for k=1:nzd
    for j=1:n
        for i=1:nxd
            
            
            e=e+1;
            ien(1,e)=i+(j-1)*nxd+(k-1)*(n+1)*nxd;
            ien(2,e)=(n+1)*nxd*nzd+(n)*(nxd+1)*(k-1)+1+(i-1)+(j-1)*(nxd+1);
            E(e)=E1;
            e=e+1;
            ien(1,e)=(n+1)*nxd*nzd+1+(n)*(nxd+1)*(k-1)+1+(i-1)+(j-1)*(nxd+1);
            ien(2,e)=i+(j-1)*nxd+(k-1)*(n+1)*nxd;
            E(e)=E2;
            e=e+1;
            ien(1,e)=i+nxd+(j-1)*nxd+(k-1)*(n+1)*nxd;
            ien(2,e)=(n+1)*nxd*nzd+(n)*(nxd+1)*(k-1)+1+(i-1)+(j-1)*(nxd+1);
            E(e)=E1;
            e=e+1;
            ien(1,e)=(n+1)*nxd*nzd+1+(n)*(nxd+1)*(k-1)+1+(i-1)+(j-1)*(nxd+1);
            ien(2,e)=i+nxd+(j-1)*nxd+(k-1)*(n+1)*nxd;
            E(e)=E2;
            
        end
        
    end
    end
    
    for k=1:nzd
    for j=1:n+1
        for i=1:nxd-1
            
            e=e+1;
            ien(1,e)=i+(j-1)*nxd+(k-1)*(n+1)*nxd;
            ien(2,e)=i+1+(j-1)*nxd+(k-1)*(n+1)*nxd;
            E(e)=E1;
        end
    end
    end
    
    for k=1:nzd
        for j=1:n
            for i=2:nxd+1
                e=e+1;
                ien(1,e)=i-1+rn*nxd*nzd+(j-1)*(nxd+1)+(n)*(nxd+1)*(k-1);
                ien(2,e)=i+rn*nxd*nzd+(j-1)*(nxd+1)+(n)*(nxd+1)*(k-1);
                 E(e)=E2;
                
            end
        end
    end
    
   for k=2:nzd
    for j=1:n
        for i=1:nxd+1
            e=e+1;
            ien(1,e)=(n+1)*nxd*nzd+(n)*(nxd+1)*(k-2)+1+(i-1)+(j-1)*(nxd+1);
            ien(2,e)=(n+1)*nxd*nzd+(n)*(nxd+1)*(k-1)+1+(i-1)+(j-1)*(nxd+1);
            E(e)=E1;
        end
    end
   end 
    
   
    for k=2:nzd
    for j=1:n+1
        for i=1:nxd
            
            e=e+1;       
            ien(1,e)=i+(j-1)*nxd+(k-2)*(n+1)*nxd;
            ien(2,e)=i+(j-1)*nxd+(k-1)*(n+1)*nxd;
            E(e)=E1;
            
        end
    end
    end 

nel=size(ien,2);
   
   
%%%%%%%%%%%%%%%%%%%%%%%
% Material Properties %
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% Material %
%%%%%%%%%%%%

% E(1:nel)=1;
 for i=1:nzd
E(4*n+(i-1)*4*n^2:4*n:4*n^2+(i-1)*4*n^2)=1; %In beam case, simple model need to ensure E=I=1;

 E(4*n-2+(i-1)*4*n^2:4*n:4*n^2-2+(i-1)*4*n^2)=1;
 
 E((n-1)*(n+1)*nzd+1+(4*n^2+4*n^2)+(i-1)*n^2:n:(n-1)*(n+1)*nzd+1+(4*n^2+4*n^2)+n^2+(i-1)*n^2)=1;
 end

A(1:nel)=1e3;
EA=E(1)*A(1);

Iz(1:nel)=1;
Iy(1:nel)=1;
G=1;
J=1;
psi=0;
                       

%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%
% prescribed temperature (essential boundary condition)
%
% idb(i,N)=1 if the degree of freedom i of the node N is prescribed
%         =0 otherwise
%
% 1) initialize idb to 0


idb=zeros(ndf,nnp);
idbx=zeros(ndf,nnp);


% for i=1:nzd
% idb(1,1+(i-1)*(nxd+1)*nxd:nxd+(i-1)*(nxd+1)*nxd) = 1;
% idb(2,1+(i-1)*(nxd+1)*nxd:nxd+(i-1)*(nxd+1)*nxd) = 1;
% 
% idb(1,1+nxd*nxd+(i-1)*(nxd+1)*nxd:(nxd+1)*nxd+(i-1)*(nxd+1)*nxd) = 1;
% idb(2,1+nxd*nxd+(i-1)*(nxd+1)*nxd:(nxd+1)*nxd+(i-1)*(nxd+1)*nxd) = 1;
% end

%%%%%%%%%%%%%%%%%%%%%%%
% Up-down Boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nzd
idbx (1:6,nxd/2+(i-1)*(n+1)*nxd:nxd/2+(i-1)*(n+1)*nxd+1) = 1;    
idbx(4:6,1+(i-1)*(n+1)*nxd:nxd+(i-1)*(n+1)*nxd) = 1;
idbx(2,1+(i-1)*(n+1)*nxd:nxd+(i-1)*(n+1)*nxd) = 1;

idbx(1:6,1+n*nxd+(i-1)*(n+1)*nxd:(n+1)*nxd+(i-1)*(n+1)*nxd) = -1;

end

for i=1:nzd
idbx(3:6,1+nxd+(i-1)*(nxd+1)*nxd:nxd:1+nxd*(n-1)+(i-1)*(nxd+1)*nxd) = 1; %outer left
idbx(1,1+nxd+(i-1)*(nxd+1)*nxd:nxd:1+nxd*(n-1)+(i-1)*(nxd+1)*nxd) = 1; %outer left

idbx(3:6,nxd+nxd+(i-1)*(nxd+1)*nxd:nxd:(nxd)*nxd+(i-1)*(nxd+1)*nxd) = -1; %outer right
idbx(1,nxd+nxd+(i-1)*(nxd+1)*nxd:nxd:(nxd)*nxd+(i-1)*(nxd+1)*nxd) = -1; %outer right

idbx(3:6,(n+1)*nxd*nzd+1+(i-1)*(nxd+1)*nxd:nxd+1:(n+1)*nxd*nzd+1+nxd*n+(i-1)*(nxd+1)*nxd)=1; % Left side inner
idbx(1,(n+1)*nxd*nzd+1+(i-1)*(nxd+1)*nxd:nxd+1:(n+1)*nxd*nzd+1+nxd*n+(i-1)*(nxd+1)*nxd)=1; % Left side inner

idbx(3:6,(nxd)*(n+1)*nzd+nxd+1+(i-1)*(nxd+1)*nxd:nxd+1:(nxd+1)*(n)*nzd+nxd+(nxd)*(n)+(i-1)*(nxd+1)*nxd)=-1; % right side inner
idbx(1,(nxd)*(n+1)*nzd+nxd+1+(i-1)*(nxd+1)*nxd:nxd+1:(nxd+1)*(n)*nzd+nxd+(nxd)*(n)+(i-1)*(nxd+1)*nxd)=-1; % right side inner

idbx(1,1+nxd+(i-1)*(nxd+1)*nxd)=0; % bottom Left side outer
idbx(1,nxd+nxd+(i-1)*(nxd+1)*nxd)=0; % bottom right side outer
idbx(1,(n+1)*nxd*nzd+1+(i-1)*(nxd+1)*nxd)=0; % bottom Left side inner
idbx(1,(nxd)*(n+1)*nzd+nxd+1+(i-1)*(nxd+1)*nxd)=0; % bottom right side inner


end

idb(:,:)=abs(idbx(:,:));

% idb(1,nxd/2+nxd^2:nxd*(nxd+1):(nzd-1)*nxd*(nxd+1)+nxd/2+nxd*(nxd+1)-nxd)=1;  
% idb(1,nxd/2+nxd^2+1:nxd*(nxd+1):(nzd-1)*nxd*(nxd+1)+nxd/2+1+nxd*(nxd+1)-nxd)=1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Lateral Boundary conditions %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:nzd
% idb(1:6,1+(i-1)*(nxd+1)*nxd:nxd:1+nxd*nxd+(i-1)*(nxd+1)*nxd) = 1;
% 
% idb(1:6,(nxd+1)*nxd*nzd+1+(i-1)*(nxd-1)*nxd:nxd-1:(nxd-1)*(nxd-1)+(nxd+1)*nxd*nzd+1+(i-1)*(nxd-1)*nxd) = 1;
% 
% idb(1:6,nxd+(i-1)*(nxd+1)*nxd:nxd:(nxd+1)*nxd+(i-1)*(nxd+1)*nxd) = 1;
% 
% idb(1:6,(nxd+1)*nxd*nzd+nxd-1+(i-1)*(nxd-1)*nxd:nxd-1:(nxd-1)*(nxd-1)+(nxd+1)*nxd*nzd+nxd-1+(i-1)*(nxd-1)*nxd) = 1;
% 
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Up-Down Boundary conditions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  for i=1:nzd
%  idb(2:6,1+(i-1)*(nxd+1)*nxd:nxd+(i-1)*(nxd+1)*nxd) = 1;
%  end
% 
% 
% idb(1,nxd/2:nxd*(nxd+1):(nzd-1)*nxd*(nxd+1)+nxd/2)=1;  idb(1,nxd/2+1:nxd*(nxd+1):(nzd-1)*nxd*(nxd+1)+nxd/2+1)=1;

% 2) enter the flag for prescribed displacement boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed nodal temperature boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(i,N): prescribed displacement for the dof i of node N
% initialize g

g=zeros(ndf,nnp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Up-down displacement boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:nzd
% 
% 
% g(2,1+(i-1)*(nxd+1)*nxd:nxd+(i-1)*(nxd+1)*nxd) = -1;
% g(2,1+nxd*nxd+(i-1)*(nxd+1)*nxd:(nxd+1)*nxd+(i-1)*(nxd+1)*nxd) = 1;
% % g(1,1+(i-1)*(nxd+1)*nxd:nxd+(i-1)*(nxd+1)*nxd) = -1;
% % g(1,1+nxd*nxd+(i-1)*(nxd+1)*nxd:(nxd+1)*nxd+(i-1)*(nxd+1)*nxd) = 1;
% 
% end

 for i=1:nzd
% g(2,1+(i-1)*(nxd+1)*nxd:nxd+(i-1)*(nxd+1)*nxd) = -1;
g(2,1+n*nxd+(i-1)*(n+1)*nxd:(n+1)*nxd+(i-1)*(n+1)*nxd) = 1;
end

% enter the values

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed nodal forces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(i,N): prescribed force for the dof i of node N
% initialize f

fnn=zeros(ndf*nen,nel);
f = zeros(ndf,nnp);

%  for i=1:nzd
% % f(2,1+(i-1)*(nxd+1)*nxd:nxd+(i-1)*(nxd+1)*nxd) = -1;
% f(2,1+nxd*nxd+(i-1)*(nxd+1)*nxd:(nxd+1)*nxd+(i-1)*(nxd+1)*nxd) = 1;
% end

% f(2,nxd*n+1:nxd*n+nxd)=1;
% f(2,1:nxd)=-1;

% f(1,nxd+1:nxd:nxd*(n-2)+nxd+1)=1;
% f(1,nxd+nxd:nxd:nxd*(n-2)+nxd+nxd)=-1;


% f(1,nxd+1:nxd:nxd*(n-2)+nxd+1)=-1;
% f(1,nxd+nxd:nxd:nxd*(n-2)+nxd+nxd)=1;

% enter the values

%---------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number the equations; build the id table %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[id,neq] = number_eq(idb,nnp,ndf);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the elemental quantities in the elemental coordinate system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nintx=4;
[Ke,ke,Qe]=Ke_beam3D_finite(E,A,Iy,Iz,G,J,psi,xn,ien,nen,ndf,nsd,nel,Nintx);


% Contribution of the prescribed displacements to the elemental force vector
% fe=fe-Ke*Ue

fe=zeros(ndf*nen,nel); % fe may be non zero in general
Ue=zeros(ndf*nen,nel);
for e=1:nel
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e) = g(i,ien(n,e));
        end
    end
    fe(:,e)=fe(:,e)-Ke(:,:,e)*Ue(:,e);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly operation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------
% build K and F
%----------------

K=zeros(neq,neq); 
F=zeros(neq,1);

% input the prescribed nodal forces in F

for N=1:nnp
    for i=1:ndf
        if (id(i,N) > 0)
            P=id(i,N);
            F(P)=f(i,N);
        end
    end
end

% compute global K and F

if (neq > 0)
    for e=1:nel
        K = addstiff(K,id,Ke(:,:,e),ien(:,e),nen,ndf);
        F = addforce(F,id,fe(:,e),ien(:,e),nen,ndf);
    end
end

% Solve the system

if (neq > 0)
    U=K\F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% complete U %
%%%%%%%%%%%%%%

Ucomp=zeros(ndf,nnp);
for N=1:nnp
    for i=1:ndf
        if (id(i,N) == 0)
            Ucomp(i,N)=g(i,N);
        else
            P=id(i,N);
            Ucomp(i,N)=U(P);
        end
    end
end

% print results
% disp('Nodal Displacements:')
% disp(' node        d1              d2            ')
% for N=1:nnp 
%     disp(sprintf('%5d   %1.3e       %1.3e        %1.3e        %1.3e         %1.3e      %1.3e',N,Ucomp(:,N)))
% end
% disp(' ')

%%%%%%%%%%%%%
% REACTIONS %
%%%%%%%%%%%%%
% build the idb table; overwrite original idb table
% idb(i,N): equation number associated with dof i of node N

ineq=0; % number of equations
for i=1:ndf
    for N=1:nnp
        if (idb(i,N) > 0) % assign an equation number to all prescribed nodes
            ineq=ineq+1;
            idb(i,N)=ineq;
        end
    end
end

% Contribution of the displacement to the elemental force vector
% fe=Ke*Ue

for e=1:nel
    Ue(:,e)=zeros(ndf*nen,1);
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e)=Ucomp(i,ien(n,e));
        end
    end
    fe(:,e)=Ke(:,:,e)*Ue(:,e);
end

% compute reactions R %

R=zeros(ineq,1);
for e=1:nel
    R = addforce(R,idb,fe(:,e),ien(:,e),nen,ndf);
end

% Collect reactions

Rcomp=zeros(ndf,nnp);
for N=1:nnp
    for i=1:ndf
        if (idb(i,N) > 0)
            Rcomp(i,N)=R(idb(i,N));
        end
    end
end

% print results

% disp('Nodal Reactions')  
% disp(' node         R1              R2           ')
% for N=1:nnp
%     disp(sprintf('%5d    %1.1e       %1.1e  ',N,Rcomp(:,N)))
% end 
% disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%
% AXIAL FORCES/STRESSES %
%%%%%%%%%%%%%%%%%%%%%%%%%

for e=1:nel
    Ue(:,e)=zeros(ndf*nen,1);
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e)=Ucomp(i,ien(n,e));
        end
    end
   
    v=xn(:,ien(2,e))-xn(:,ien(1,e));% calculate the coordinates of elmemnts
    L=norm(v);
    [B] = B_truss(L);
    

    epsn(e)=B*[Qe(1,1:2,e),Qe(1,4:5,e);...
        Qe(4,1:2,e),Qe(4,4:5,e);]*[Ue(1:2,e);Ue(7:8,e)];
    
    str(e)=epsn(e)*E(e);
    
    fe(:,e) = Ke(:,:,e)*Ue(:,e)-fnn(:,e); %Global force vector (**no relation with area**)
    force(:,e)=Qe(:,:,e)*fe(:,e)-fnn(:,e); %Local force vector (**no relation with area**)
    
    
% [point,w]=gauss(Nintx,0,0,nsd);
%      for i=1:Nintx
%         Xi=point(i,1);
%         Wi(i,1)=w(i,1);
%         [Bb,Bb2] = B_beam(Xi,L);  
%         epstz(e)=[Bb(1,2),Bb(1,4)]*[Qe(3,3,e),0;0,Qe(6,6,e)]*[Ue(6,e);Ue(12,e)];
%         epsty(e)=[Bb2(1,2),Bb2(1,4)]*[Qe(9,9,e),0;0,Qe(12,12,e)]*[Ue(5,e);Ue(11,e)];
%         strtz(e)=epstz(e)*E(e);
%         strty(e)=epsty(e)*E(e);
%         
%      end
%       


end
ii=1;
fracture_flag = zeros(2,nel);

nelv(ii)=nel;
xn_x(:,:,ii)=xn(:,:);
ien_x{ii}=ien;
E_x{ii}=E;
Ucomp_x(:,:,ii) = Ucomp;  
force_x{ii}=fe;
Qe_x{ii}=Qe;

g_main = g;

for lambda = 0.1:0.1:1%0.1*Length*2 % the maximum of all pinned boundary is 128
    lambda
    

    g = lambda *g_main;

    Ke=zeros(12,12,nel);
    ke=zeros(12,12,nel);
    Qe=zeros(12,12,nel);
         

Nintx=4;
[Ke,ke,Qe]=Ke_beam3D_finite(E_x{ii},A,Iy,Iz,G,J,psi,xn,ien_x{ii},nen,ndf,nsd,nelv(ii),Nintx);

% Contribution of the prescribed displacements to the elemental force vector
% fe=fe-Ke*Ue

fe=zeros(ndf*nen,nel); % fe may be non zero in general
Ue=zeros(ndf*nen,nel);
for e=1:nel
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e) = g(i,ien(n,e));
        end
    end
    fe(:,e)=fe(:,e)-Ke(:,:,e)*Ue(:,e);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly operation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------
% build K and F
%----------------

K=zeros(neq,neq); 
F=zeros(neq,1);

    
% input the prescribed nodal forces in F

for N=1:nnp
    for i=1:ndf
        if (id(i,N) > 0)
            P=id(i,N);
            F(P)=f(i,N);
        end
    end
end

% compute global K and F

    for e=1:nel
        K = addstiff(K,id,Ke(:,:,e),ien(:,e),nen,ndf);
        F = addforce(F,id,fe(:,e),ien(:,e),nen,ndf);
    end


% Solve the system

    U=K\F;

            
            Ucomp=zeros(ndf,nnp);
            fe = zeros(ndf*nen,nel);
            
            for N=1:nnp
                for i=1:ndf
                    if (id(i,N) == 0)
                        Ucomp(i,N)=g(i,N);
                    else
                        P=id(i,N);
                        Ucomp(i,N)=U(P);
                    end
                end
            end
            
             ineq=0; % number of equations
            for i=1:ndf
                for N=1:nnp
                    if (idb(i,N) > 0) % assign an equation number to all prescribed nodes
                        ineq=ineq+1;
                        idb(i,N)=ineq;
                    end
                end
            end
            
            % Contribution of the displacement to the elemental force vector
            % fe=Ke*Ue
            %
            for e=1:nel
                Ue(:,e)=zeros(ndf*nen,1);
                for n=1:nen
                    for i=1:ndf
                        Ue(i+(n-1)*ndf,e)=Ucomp(i,ien(n,e));
                    end
                end
                fe(:,e)=Ke(:,:,e)*Ue(:,e);
            end
            
            % compute reactions R %
            
            R=zeros(ineq,1);
            for e=1:nel
                R = addforce(R,idb,fe(:,e),ien(:,e),nen,ndf);
            end
            
            % Collect reactions
            
            Rcomp=zeros(ndf,nnp);
            for N=1:nnp
                for i=1:ndf
                    if (idb(i,N) > 0)
                        Rcomp(i,N)=R(idb(i,N));
                    end
                end
            end
            
            Ue=zeros(ndf*nen,nel);
            str=zeros(1,nel);
            force=zeros(12,nel);
            
             for e=1:nel
                 Ue(:,e)=zeros(ndf*nen,1);
                 for n=1:nen
                     for i=1:ndf
                         Ue(i+(n-1)*ndf,e)=Ucomp(i,ien(n,e));
                     end
                 end               
                 
                 epsn(e)=B*[Qe(1,1:2,e),Qe(1,4:5,e);...
                     Qe(4,1:2,e),Qe(4,4:5,e);]*[Ue(1:2,e);Ue(7:8,e)];
                 
                 str(e)=epsn(e)*E(e);
                 
                 fe(:,e) = Ke(:,:,e)*Ue(:,e)-fnn(:,e); %Global force vector (**no relation with area**)
                 force(:,e)=Qe(:,:,e)*fe(:,e)-fnn(:,e); %Local force vector (**no relation with area**)
                 
             end
  
    
    ii=ii+1;
    

    fracture_flagR = zeros(2,nel);
    fracture_flag = zeros(2,nel);
    
    for e=1:nel
        if abs(str(e))>=0.1*E(e)%check local stresses if exceed yield stresses
            fracture_flag(:,e) = 1;
        end
        
    end
    fracture_flagR(:,:)=  fracture_flag(:,:)<1;
    ien( :, ~any(fracture_flagR,1)) = [];  %delete columns
    E( :, ~any(fracture_flagR,1)) = [];
    nel=size(ien,2);
   
    nelv(ii)=nel;
    
    E_x{ii}=E;
    ien_x{ii}=ien;
    force_x{ii}=fe;
    Qe_x{ii}=Qe;
    
    for i=1:nzd
        re(:,i)=sum(Rcomp(2,1+nnn*nxd+(i-1)*(nnn+1)*nxd:(nnn+1)*nxd+(i-1)*(nnn+1)*nxd));
    end
    
    Reaction(ii)=abs(sum(re)/(nxd*nzd));
    
    Ucomp_x(:,:,ii) = Ucomp;
%     xn_x(:,:,ii)=xn(:,:)+Ucomp(1:3,:);
    
    lmb(ii) = lambda;
   
end
 
      
lmb(1) = 0;
Reaction(1) = 0;
[val, idx] = max(Reaction);
dd=plot(lmb,Reaction,'r','LineWidth',2); 
xlabel('Elongation(δ)')
ylabel('Reaction(F)')
title('Reaction and Elongation of 3D Beam');
datatip(dd,'DataIndex',nnz(0.1:0.1:0.1*Length)+1);
datatip(dd,'DataIndex',idx);

Length    
EA
Area = trapz(lmb, Reaction)




%%%%%%%%%%%%%%%%%%%%
% plot the results %
%%%%%%%%%%%%%%%%%%%%
 plot_results_2021('beam',xn,f,idbx,Ucomp_x,Rcomp,ien_x,nelv,nen,nsd,ndf,nnp,force_x,fnn,Qe_x,0,g,E_x);

