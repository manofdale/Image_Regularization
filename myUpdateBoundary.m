% A. Goze Polat, 1631092, CENG566, HW2
% Due: Sat 20.00

function I = myUpdateBoundary(J)

[r,c]=size(J);
I=J;
I(1,:)=J(2,:); %duplicate first row of image J in the new image

I(r,:)=J(r-1,:);

I(:,1)=I(:,2); 

I(:,c)=I(:,c-1); 

end