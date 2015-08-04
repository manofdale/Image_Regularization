% A. Goze Polat 1631092
function [u1,u2,u3]=myExplicitPM(u,h,niter,lambda,deltaT,blurSig)
%u=imread('LenaDark64.png');
sz=size(u,3);
if(sz == 3)  		% an rgb image
	u=rgb2gray(u);
	u=ind2gray(u,gray(1000));
else
	u=mat2gray(u);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1),imshow(u),title('original image');

% Create the gaussian filter
G = fspecial('gaussian',[blurSig*2+1 blurSig*2+1],blurSig);
%u=imfilter(u,fspecial('gaussian',[3,3],1),'same');
figure(2),imshow(imfilter(u,G,'same')),title('smoothed image');
u1=u;
[m,n]=size(u1);

u1=reshape(u1,m*n,1);
 if strcmp(class(u),'uint8')
    u1=double(u1);
end
means=zeros(niter);
%% uncomment below for PM 1 
%
disp('Perona malik type 1 explicit scheme');
%number of iterations n
for i=1:niter,
	u2=reshape(u1,m,n);
	
	figure(i+2),imshow(u2),title(strcat("perona-malik1-step: ",num2str(i-1)));
	% use central differences for gradient
	u3=myUpdateBoundary(blkdiag(0,u2,0)); % 
	imfilter(u3,G,'same');
	u3=u3(2:m+1,2:n+1);
	g=myDiffusivity2(myCentralDiff(u3,h),lambda,'perona-malik1',0);
	u1=myA(u2,g,h,deltaT,1)*u1;
	means(i)=mean(u1); % is this not supposed to be the same?
end
u1=reshape(u1,m,n);
figure(niter+3),imshow(u1),title('modified image perona-malik1');
figure(niter+4),plot(means),title('mean perona-malik1');

disp('Continue for the next scheme..'),pause()
disp('Perona malik type 2 explicit scheme');
u1=reshape(u,m*n,1);
 if strcmp(class(u),'uint8')
    u1=double(u1);
end
%
means=zeros(niter);
disp('Perona malik type 2 (exp) explicit scheme');
%number of iterations n
for i=1:niter,
	u2=reshape(u1,m,n);
	
	figure(i+2),imshow(u2),title(strcat("perona-malik2-step: ",num2str(i-1)));
	% use central differences for gradient
	u3=myUpdateBoundary(blkdiag(0,u2,0)); % 
	imfilter(u3,G,'same');
	u3=u3(2:m+1,2:n+1);
	g=myDiffusivity2(myCentralDiff(u3,h),lambda,'perona-malik2 (exp)',0);
	u1=myA(u2,g,h,deltaT,1)*u1;
	means(i)=mean(u1); % is this not supposed to be the same?
end
u1=reshape(u1,m,n);
figure(niter+3),imshow(u1),title('modified image perona-malik2');
figure(niter+4),plot(means),title('mean perona-malik2');;
%%%
disp('Continue for the next scheme..'),pause()
disp('Charbonnier explicit scheme');
u1=reshape(u,m*n,1);
 if strcmp(class(u),'uint8')
    u1=double(u1);
end
means=zeros(niter);
%number of iterations n
for i=1:niter,
	u2=reshape(u1,m,n);	
	figure(i+2),imshow(u2),title(strcat("charbonnier-step: ",num2str(i-1)));
	% use central differences for gradient
	u3=myUpdateBoundary(blkdiag(0,u2,0)); % 
	g=myDiffusivity2(myCentralDiff(u2,h),lambda,'charbonnier',0);
	u1=myA(u2,g,h,deltaT,1)*u1;
	means(i)=mean(u1); % is this not supposed to be the same?
end
u1=reshape(u1,m,n);
figure(niter+3),imshow(u1),title('final image charbonnier semi-implicit');
figure(niter+4),plot(means),title('mean charbonnier semi-implicit');

disp('Continue for the next scheme..'),pause()
disp('Perona malik type 1 semi-implicit scheme');
u1=reshape(u,m*n,1);
 if strcmp(class(u),'uint8')
    u1=double(u1);
end
means=zeros(niter);
%number of iterations n
for i=1:niter,
	u2=reshape(u1,m,n);
	
	figure(i+2),imshow(u2),title(strcat("perona-malik1-step: ",num2str(i-1)));
	% use central differences for gradient
	u3=myUpdateBoundary(blkdiag(0,u2,0)); % 
	imfilter(u3,G,'same');
	u3=u3(2:m+1,2:n+1);
	g=myDiffusivity2(myCentralDiff(u3,h),lambda,'perona-malik1',0);
	u1=inv(myA(u2,g,h,deltaT,-1))*u1;
	means(i)=mean(u1); % is this not supposed to be the same?
end
u1=reshape(u1,m,n);
figure(niter+3),imshow(u1),title('final image perona-malik1 semi-implicit');
figure(niter+4),plot(means),title('mean perona-malik1 semi-implicit');

disp('Continue for the next scheme..'),pause()
disp('Perona malik type 2 semi-implicit scheme');
u1=reshape(u,m*n,1);
 if strcmp(class(u),'uint8')
    u1=double(u1);
end
means=zeros(niter);
%number of iterations n
for i=1:niter,
	u2=reshape(u1,m,n);
	
	figure(i+2),imshow(u2),title(strcat("perona-malik2-step: ",num2str(i-1)));
	% use central differences for gradient
	u3=myUpdateBoundary(blkdiag(0,u2,0)); % 
	imfilter(u3,G,'same');
	u3=u3(2:m+1,2:n+1);
	g=myDiffusivity2(myCentralDiff(u3,h),lambda,'perona-malik2 (exp)',0);
	u1=inv(myA(u2,g,h,deltaT,-1))*u1;
	means(i)=mean(u1); % is this not supposed to be the same?
end
u1=reshape(u1,m,n);
figure(niter+3),imshow(u1),title('final image perona-malik2 semi-implicit');
figure(niter+4),plot(means),title('mean perona-malik2 semi-implicit');

disp('Continue for the next scheme..'),pause()
disp('Charbonnier semi-implicit scheme');
u1=reshape(u,m*n,1);
 if strcmp(class(u),'uint8')
    u1=double(u1);
end
means=zeros(niter);
%number of iterations n
for i=1:niter,
	u2=reshape(u1,m,n);
	
	figure(i+2),imshow(u2),title(strcat("charbonnier-step: ",num2str(i-1)));
	% use central differences for gradient
	u3=myUpdateBoundary(blkdiag(0,u2,0)); % 
	g=myDiffusivity2(myCentralDiff(u2,h),lambda,'charbonnier',0);
	u1=inv(myA(u2,g,h,deltaT,-1))*u1;
	means(i)=mean(u1); % is this not supposed to be the same?
end
u1=reshape(u1,m,n);
figure(niter+3),imshow(u1),title('final image charbonnier semi-implicit');
figure(niter+4),plot(means),title('mean charbonnier semi-implicit');
%fix 1: for charbonnier, smoothing is not necessary
%{
u2=reshape(u,m*n,1);
for i=1:n,	
	g=myDiffusivity(imfilter(u2,G,'same'),lambda,'perona-malik1',0);% x,lambda,gType,plotMe
	u2=myA(u2,g,deltaT)*u2;
end
u2=reshape(u2,m,n);
figure(),imshow(u2),title('modified image perona-malik1');
%}
end
