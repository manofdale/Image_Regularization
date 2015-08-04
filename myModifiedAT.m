%hw5, A. Goze Polat, 1631092
%Modified Ambrosio Tortelli approximation to Mumford Shah
function myModifiedAT(imName,noiseLvl,alphA,betA,deltaT,rhO,epsiloN,niter)
%%%%%%%%%%%%%%%% init u %%%%%%%%%%%%%%%%%%%%%%%%%
% u0=imnoise(u, 'salt & pepper',0.15); % add it outside
u=imread(imName);
if noiseLvl != 0
	u=imnoise(u, 'salt & pepper',noiseLvl);
end
u0=u;
Cstr=strsplit(imName,'.');
imName=Cstr{1};%discard .png
sz=size(u0,3);
if(sz == 3)  		% an rgb image
	u0=rgb2gray(u0);
	u0=ind2gray(u0,gray(1000));
else
	u0=mat2gray(u0);
end
if strcmp(class(u0),'uint8')
    u0=double(u0);
end
%%%%%%%%%%%%%%%%%% find us %%%%%%%%%%%%%%%%%%%%%%%%%
disp('initializing..');
[m,n]=size(u0);
us=reshape(myCentralDiff(u0,1),m*n,1);%h is 1

%%%%%%%%%%%%%%%%% initialize v %%%%%%%%%%%%%%%%%%%%%
% better than random init:
v=1.0-1.0./(1.0+us*(2.0*alphA*rhO));

%%% sigmA = #iterations for diffusion step of u %%%%
sigmA=ceil(sqrt(2*alphA/betA));
%%% rhOn = #iterations for diffusion step of v %%%%%
rhOn=ceil(rhO);

uk=reshape(u0,m*n,1);
uk0=(deltaT*betA/alphA)*uk;% == u0 or g
%constants to be used
c=1.0/(1.0+deltaT*betA/alphA);
b=1+deltaT/rhO^2;
d=2*deltaT*alphA/rhO;

%% if the image is large, can this save time????????
%%%%%%%%%%%%%%%%% initialize L %%%%%%%%%%%%%%%%%%%%%
titL=strcat('myL_',num2str(m),'x',num2str(n),'.mat');
if exist(strcat('myL_',num2str(m),'x',num2str(n),'.mat')) == 0
	disp('Calculating L..');
	L=myLaplace(m,n);%I+deltat*L
	%% save it to a file for later usage %%
	save(titL,'L');
else 
	load(titL,'L');
end
%%%%%%%%%%%%%%%%% initialize A %%%%%%%%%%%%%%%%%%%%%
titlA=strcat('myA_',num2str(m),'x',num2str(n),'_dt_',num2str(deltaT),'.mat');
if exist(strcat('myA_',num2str(m),'x',num2str(n),'_dt_',num2str(deltaT),'.mat')) == 0
	disp('Calculating A..');
	A=myA(m,n,deltaT);%I+deltat*L
	%% save it to a file for later usage %%
	save(titlA,'A');
else 
	%disp('Already calculated A, loading from file..');
	load(titlA,'A');
end
% iterate n times, if converged, break %%%%%%%%%%%%%
disp('iterating..');
deltaTover2=deltaT*0.5;
md=1+floor(niter/11);
subFigN=1;
for i=1:niter,
	%%% for sigmA times diffuse u %%%
	u_temp=uk;
	ttl1=num2str(i);%strcat('delta_t',num2str(deltaT),'_b',num2str(betA),'_i',num2str(niter));
	ttl2=ttl1; %%%%%%%%%%%!!!!!!!!
	ttl3=ttl1; %%%%%%%%%%%!!!!!!!!
	if mod(i,md)==0 || i == 1
		figure(1),subplot(3,4,subFigN),imshow(reshape(uk,m,n,1)),title(ttl1);
	end
	g=(1-v).^2;
	Lg=L*g;
	for j=1:sigmA,
		uk1=(uk0+uk+(g.*(L*uk)-uk.*Lg+L*(g.*uk))*deltaTover2)*c;
		uk=uk1;
	end
	us=reshape(myCentralDiff(reshape(uk,m,n),1),m*n,1);
	%B=rhO^2./(b-us*d);
	dB=1.0./(us*d+1+deltaT/rhO^2);%diagonal vec of inv(diag(B))
	%%% for rhOn times diffuse v %%%
	for j=1:rhOn,
		vk1=dB.*(A*v+us*d);%faster than inverse %%%
		v=vk1;
	end
	if mod(i,md)==0 || i == 1
		figure(2),subplot(3,4,subFigN),imshow(reshape(v,m,n,1)),title(ttl2);
		uDiff=abs(u_temp-uk);
		figure(3),subplot(3,4,subFigN),imshow(reshape(uDiff,m,n)),title(ttl3);
		subFigN=subFigN+1;
	end
	%break if converged	
	if abs(u_temp-uk) < abs(u_temp)*epsiloN 
		disp('converged');		
		break;
	end
end
%fill the remaining subplots
for k=subFigN:12
			figure(1),subplot(3,4,k),imshow(reshape(uk,m,n,1)),title(ttl1);
			figure(2),subplot(3,4,k),imshow(reshape(v,m,n,1)),title(ttl2);
			figure(3),subplot(3,4,k),imshow(reshape(uDiff,m,n)),title(ttl3);
end
	figttl=strcat(imName,'a',num2str(alphA),'b',num2str(betA),'d',num2str(deltaT),'r',num2str(rhO),'n',num2str(niter));
	figU=figure(1);	  % to observe u
	saveas(figU,strcat(figttl,'U'),'png');
	figV=figure(2);   % to observe v
	saveas(figV,strcat(figttl,'V'),'png');
	figDiff=figure(3);% to observe convergence & what is lost
	saveas(figDiff,strcat(figttl,'Diff'),'png');
end % of function myModifiedAT
	
%uk1=(deltaT,alphA,betA,uk,v)
%us2=

%v=(uk,deltaT,rho,alpha,vn)
%vn=v;
%u0
%end
