%calculate A
function A=myA(u,g,h,deltaT,sgn)
	[m,n]=size(u);
	A=zeros(m*n);
	h2=h^2;	
	for i=1:m*n,
	     if mod(i,250)==0.0
	     i	     
	     end %display 
	    % y=1+floor((i-1)/m);	
	    % x=1+mod((i-1),m);
		for j=i:m*n,
			A(i,j)=myDeltaN(deltaT,g,h2,i,j,m,n,sgn);
			%%x1=1+mod((j-1),m);			
	  		%%y1=1+floor((i-1)/m);
	  		
			%%if k==0.0
			%%else
			% below was wrong:
			%5A(i,j)=k*(g(x,y)+g(x1,y1))/2.0;
			%%end
		end
	end
	%% A is symmetrical
	for i=1:m*n,
		for j=1:i,
			A(i,j)=A(j,i);
		end
	end
end
