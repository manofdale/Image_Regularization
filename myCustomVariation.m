function myCustomVariation(J)
	[m,n]=size(J);
	K=J;
	M=J;
	I=blkdiag(0,J,0);
	I=myUpdateBoundary(I);
	disp('normal variaton approx:')
	sum1=0;
	for i=2:m+1,
		for j=2:n+1,
			k=abs(I(i+1,j)-I(i,j))+abs(I(i,j+1)-I(i,j));
			sum1+=k;
			K(i-1,j-1)=k;
			%disp(strcat(num2str(i),',',num2str(j),': ',num2str(k)))
		end
	end
	K
	disp(num2str(sum1));
	sum1=0;
	disp('custom variaton approx:')
	for i=2:m+1,
		for j=2:n+1,
			k=abs(I(i+1,j)-I(i,j))+abs(I(i,j+1)-I(i,j))+abs(I(i,j-1)-I(i,j))+abs(I(i,j)-I(i-1,j));
			sum1+=k;
			M(i-1,j-1)=k;
			%disp(strcat(num2str(i),',',num2str(j),': ',num2str(k)))
		end
	end
	M
	disp(num2str(sum1));
end
