function [B,B2] = B_beam(Xi,L)


B=zeros(1,4);

      B(1,1:4)=(1/L)*[6*Xi/L,3*Xi-1,-6*Xi/L,3*Xi+1];   
     
      B2(1,1:4)=(1/L)*[-6*Xi/L,3*Xi-1,6*Xi/L,3*Xi+1];  
end

