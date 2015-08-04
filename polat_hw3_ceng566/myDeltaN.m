% A. Goze Polat 1631092
% find index of A(i,j)
% g=myDiffusivity(myGaussianBlur(u),) %m=rn,%n=cn
% h2=h^-2
% [m,n]=size(u)
% sgn =+1 if explicit, -1 if semiImplicit
function k=myDeltaN(deltaT,g,h2,i,j,m,n,sgn)
%disp('start myDeltaN')
	nDeltaT=4;	
	y=1+floor((j-1)/m);	
	x=1+mod((j-1),m);
	y0=1+floor((i-1)/m);
	x0=1+mod((i-1),m);
	k=0.0;
	g_i=g(x0,y0);
	g_j=g(x,y);
	if y==1 		%left	
		nDeltaT=nDeltaT-1;
	end
	if x==1 	   	%up
		nDeltaT=nDeltaT-1;
	end
	if x==m    		%down
		nDeltaT=nDeltaT-1;
	end
	if y==n		   	%right	
		nDeltaT=nDeltaT-1;
	end
	
	if i==j
	%-Sum_from_n=1_to_2 Sum_l_in_Nn(i) (g_i+g_j)/2h^2
		k=1.0-sgn*deltaT*nDeltaT*h2*(g_i+g_j)/2; %%
		%k=k*g(x,y);
	else if isNeighBour(m,n,i,j) == 1
	%(g_i+g_j)/2h^2
		k=sgn*deltaT*h2*(g_i+g_j)/2;
	end
	%% below was utterly wrong!!
	%end
	%if k==0
	%do nothing
	%else
	%k=k*g(x,y);%??
	%end
end
