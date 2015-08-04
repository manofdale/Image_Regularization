function k=isNeighBour(m,n,i,j)
%disp('start isNeighBour')
	k=0;
	if abs(i-j) == m % same row
		k=1;
	else if abs(i-j) ==1 % difference is 1
		if floor((j-1)/m) == floor((i-1)/m)% same column
			k=1;
		end
	end
%disp('end isNeighBour')
end
