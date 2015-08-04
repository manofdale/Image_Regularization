function y = myDiffusivity(x,lambda,gType,plotMe)
%disp('start myDiffusivity')
switch lower(gType)
    case 'charbonnier'		%  1/ sqrt ( 1+| del u| ^2 / lambda ^2)
        y = 1./sqrt(1+x.^2/lambda^2); 
    case 'perona-malik1'   	% 1/(1+|del u|^2/lambda^2)
	y = 1./(1+x.^2/lambda^2);
    case 'perona-malik2 (exp)'	% exp(-|del u|^2/2sig^2)
    	y = exp(-1*x.^2/(2*lambda^2));    
    otherwise
        disp('Unknown method.')
        [r,c]=size(x);
        y=ones(r,c);
        gType='unkown method';
end
	if plotMe == 1
		figure(), plot(x,y),title(strcat("Diffusivity: ",gType));
	end
%disp('end myDiffusivity')
end

