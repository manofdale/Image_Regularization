%prepare an m*nxm*n matrix to use as discrete Laplace operator
function L=myLaplace(m,n)
	L=sparse(m*n,m*n);
	%%% prepare diagonal %%%
	% edges -3
	d=-3.0*ones(m,n);
	% inside -4
	d(2:m-1,2:n-1)=-4.0*ones(m-2,n-2);
	% corners -2
	d(m,1)=-2.0;d(1,n)=-2.0;d(m,n)=-2.0;d(1,1)=-2.0;
	%%% prepare neighbors %%
	for i=1:n
		for j=1:m
			ix=(i-1)*m+j;
			% i==j %%%%%%			
			L(ix,ix)=d(j,i);
			% up %%%%%%%%
			if j>1
				jx=ix-1;
				L(ix,jx)=1.0;
			end
			% down %%%%%%
			if j<m
				jx=ix+1;
				L(ix,jx)=1.0;
			end	
			% left %%%%%%
			if i>1
				jx=ix-m;
				L(ix,jx)=1.0;
			end
			% right %%%%%
			if i<n
				jx=ix+m;
				L(ix,jx)=1.0;	
			end
		end
	end
end
