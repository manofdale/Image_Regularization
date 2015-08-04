%take directly x2=x^2 
function y = myDiffusivity2(x2,lambda,gType,plotMe)
%disp('start myDiffusivity')
switch lower(gType)
    case 'charbonnier'		
        y = 1./sqrt(1+x2/lambda^2); 
    case 'perona-malik1'   	
	y = 1./(1+x2/lambda^2);
    case 'perona-malik2 (exp)'	
    	y = exp(-1*x2/(2*lambda^2));    
    otherwise
        disp('Unknown method.')
        [r,c]=size(x2);
        y=ones(r,c);
        gType='unkown method';
end
	if plotMe == 1
		figure(), plot(sqrt(x),y),title(strcat("Diffusivity: ",gType));
	end
%disp('end myDiffusivity')
end

