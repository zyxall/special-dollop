
function [point,weight]=gauss(nglx,ngly,nglz,nsd)

   if nsd==1

      [point,weight]=gauss_1d(nglx);

   elseif nsd==2

      [point,weight]=gauss_2d(nglx,ngly);

   else

      [point,weight]=gauss_3d(nglx,ngly,nglz);

   end


function [point1,weight1]=gauss_1d(ngl)

%-------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients
%     of Gauss-Legendre quadrature for one-dimensional integration
%
%  Synopsis:
%     [point1,weight1]=gauss_1d(ngl)
%
%  Variable Description:
%     ngl - number of integration points
%     point1 - vector containing integration points
%     weight1 - vector containing weighting coefficients
%-------------------------------------------------------------------

%  initialization

   point1=zeros(ngl,1);
   weight1=zeros(ngl,1);

%  find corresponding integration points and weights

 if ngl==1           % 1-point quadrature rule
    point1(1)=0.0;
    weight1(1)=2.0;

 elseif ngl==2       % 2-point quadrature rule
    point1(1) = -0.577350269189626;
    point1(2) = -point1(1);
    weight1(1) = 1.0;
    weight1(2) = weight1(1);

 elseif ngl==3       % 3-point quadrature rule
    point1(1)=-0.774596669241483;
    point1(2)=0.0;
    point1(3)=-point1(1);
    weight1(1)=0.555555555555556;
    weight1(2)=0.888888888888889;
    weight1(3)=weight1(1);

 elseif ngl==4       % 4-point quadrature rule
    point1(1)=-0.861136311594053;
    point1(2)=-0.339981043584856;
    point1(3)=-point1(2);
    point1(4)=-point1(1);
    weight1(1)=0.347854845137454;
    weight1(2)=0.652145154862546;
    weight1(3)=weight1(2);
    weight1(4)=weight1(1);

else                 % 5-point quadrature rule
    point1(1)=-0.906179845938664;
    point1(2)=-0.538469310105683;
    point1(3)=0.0;
    point1(4)=-point1(2);
    point1(5)=-point1(1);
    weight1(1)=0.236926885056189;
    weight1(2)=0.478628670499366;
    weight1(3)=0.568888888888889;
    weight1(4)=weight1(2);
    weight1(5)=weight1(1);

end


function [point2,weight2]=gauss_2d(nglx,ngly)

%-------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients
%     of Gauss-Legendre quadrature for two-dimensional integration
%
%  Synopsis:
%     [point2,weight2]=gauss_2d(nglx,ngly)
%
%  Variable Description:
%     nglx - number of integration points in the x-axis
%     ngly - number of integration points in the y-axis
%     point2 - vector containing integration points
%     weight2 - vector containing weighting coefficients
%-------------------------------------------------------------------

%  determine the largest one between nglx and ngly

   if nglx > ngly
      ngl=nglx;
   else
      ngl=ngly;
   end

%  initialization

   point2=zeros(ngl,2);
   weight2=zeros(ngl,2);

%  find corresponding integration points and weights

 [pointx,weightx]=gauss_1d(nglx);     % quadrature rule for x-axis
 [pointy,weighty]=gauss_1d(ngly);     % quadrature rule for y-axis

%  quadrature for two-dimension

 for intx=1:nglx                     % quadrature in x-axis
   point2(intx,1)=pointx(intx);
   weight2(intx,1)=weightx(intx);
 end

 for inty=1:ngly                     % quadrature in y-axis
   point2(inty,2)=pointy(inty);
   weight2(inty,2)=weighty(inty);
 end


function [point3,weight3]=gauss_3d(nglx,ngly,nglz)

%-------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients
%     of Gauss-Legendre quadrature for three-dimensional integration
%
%  Synopsis:
%     [point3,weight3]=gauss_3d(nglx,ngly,nglz)
%
%  Variable Description:
%     nglx - number of integration points in the x-axis
%     ngly - number of integration points in the y-axis
%     nglz - number of integration points in the z-axis
%     point3 - vector containing integration points
%     weight3 - vector containing weighting coefficients
%-------------------------------------------------------------------

%  determine the largest one between nglx and ngly

   if nglx > ngly
     if nglx > nglz
       ngl=nglx;
     else
       ngl=nglz;
     end
   else
     if ngly > nglz
       ngl=ngly;
     else
       ngl=nglz;
     end
   end

%  initialization

   point3=zeros(ngl,3);
   weight3=zeros(ngl,3);

%  find corresponding integration points and weights

 [pointx,weightx]=gauss_1d(nglx);     % quadrature rule for x-axis
 [pointy,weighty]=gauss_1d(ngly);     % quadrature rule for y-axis
 [pointz,weightz]=gauss_1d(nglz);     % quadrature rule for z-axis

%  quadrature for two-dimension

 for intx=1:nglx                     % quadrature in x-axis
   point3(intx,1)=pointx(intx);
   weight3(intx,1)=weightx(intx);
 end

 for inty=1:ngly                     % quadrature in y-axis
   point3(inty,2)=pointy(inty);
   weight3(inty,2)=weighty(inty);
 end

 for intz=1:nglz                     % quadrature in z-axis
   point3(intz,3)=pointz(intz);
   weight3(intz,3)=weightz(intz);
 end

