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

%% some old code in the class 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% for i=1:1:r,
%%	for j=1:1:y,
%%		I(i+r,j+1)=J(i,j);
%% end
%% end


%% for i=2:1:r, 
%%	I(i,1)=J(i,2);
%%	I(i,c+1)=J(i,c-1);
%% end

%% for i=2:1:c,
%%        I(1,i)=J(2,i);
%%        I(x+1,i)=J(r-1,i);
%% end
