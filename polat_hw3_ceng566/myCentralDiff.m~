%gradient approximation to be used as the argument of myDiffusivity2
function y=myCentralDiff(u,h)
	[m,n]=size(u);
	y=blkdiag(0,u,0);
	h2=(2*h)^2;
	y=myUpdateBoundary(y);%%duplicate borders - (replace zeros)
	c=y;
	%find central differences
	for i=2:m+1,
		for j=2:n+1,
			c(i,j)=((y(i,j+1)-y(i,j-1))^2+(y(i+1,j)-y(i-1,j))^2)/h2;
		end
	end	
	y=c(2:m+1,2:n+1);
end
