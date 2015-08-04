%calculate A
function A=myA(m,n,deltaT)
	%A=zeros(m*n);
	A=eye(m*n)+myLaplace(m,n)*deltaT; %O(m*n)
	
	%{
	% below was too inefficient: O([m*n]^2)
	for i=1:m*n,
	     if mod(i,250)==0.0
	     i	     
	     end
		for j=i:m*n,
			A(i,j)=myDeltaN(deltaT,i,j,m,n);
		end
	end
	%% A is symmetrical
	for i=1:m*n,
		for j=1:i,
			A(i,j)=A(j,i);
		end
	end
	%}
end
