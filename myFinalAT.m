%Final Project, Ceng 566, A. Goze Polat
% This function implements the second modification
% with the commented out values
function [U,V,H]=myFinalAT(imName,noiseLvl,epsiloN,niter,deltaT,alphA,betA,rhO,omegA,initH)
	%gammA=1;
%% ?? ( +%gamma(u-v)^2 {abs {nabla h}}^2)
	%%% sigmA = #iterations for diffusion step of u %%%%
	sigmA=1; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%sigmA=ceil(sqrt(2*alphA*(1-omegA)/betA));
	%%% rhOn = #iterations for diffusion step of v %%%%%
	rhOn=1;  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%rhOn=ceil(rhO);
	gammaN=1;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	%gammaN=sigmA;
	
	u=imread(imName);
	if noiseLvl != 0
		u=imnoise(u, 'salt & pepper',noiseLvl);
	end
	Cstr=strsplit(imName,'.');
	imName=Cstr{1};	%discard .png
	imSz=size(u,3);
	if(imSz == 3)  	% an rgb image
		u0=rgb2gray(u);
		u0=ind2gray(u,gray(1000));
	else
		u0=mat2gray(u);
	end
	if strcmp(class(u0),'uint8')
	    u0=double(u0);
	end
	% finding us = smoothed u
	disp('initializing..');
	[m,n]=size(u0);
	us=reshape(myCentralDiff(u0,1),m*n,1);%h is 1
	figure(10),imshow(reshape(us,m,n)),title('initial cDiff image')
	%%%%%%%%%%%%%%%%% initialize L %%%%%%%%%%
	titL=strcat('myL_',num2str(m),'x',num2str(n),'.mat');
	if exist(strcat('myL_',num2str(m),'x',num2str(n),'.mat')) == 0
		disp('Calculating L..');
		L=myLaplace(m,n);%I+deltat*L
		%% save it to a file for later usage %%
		save(titL,'L');
	else
		load(titL,'L');
	end
	%%%%%%%%%%%%%%%%% initialize A %%%%%%%%%%%
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
	%%%%%%%%%%%%%%%%% initialize v %%%%%%%%%%%%%%%%%%%%%
	% better than random init:
	v=zeros(m*n,1);%1.0./(1.0+us*(alphA*rhO));
	%k=1000;
	%{	
	rhO2=0.5*k*rhO;
	alphA2=0.1*alphA/k;
	betA2=0.7*betA/k;
	vr=1.0./(1.0+us*(alphA2*rhO2));%1.0./(1.0+us*(alphA/1000*rhO2));
	%}
	figure(21),imshow(reshape(v,m,n)),title('initial v image')

	
	%rhO2n=ceil(rhO2);
	%constants to be used
	db_a=deltaT*betA/alphA;
	%db_a2=deltaT*betA2/alphA2;
	%c2=1.0/(1.0+db_a2);
	b=1+deltaT/rhO^2;
	%b2=1+deltaT/rhO2^2;	
	%d2=2*deltaT*alphA2/rhO2;
	%% uk0 & uk
	uk=reshape(u0,m*n,1);
	u_original=reshape(u0,m*n,1);
	h=uk;
	if initH == -1
		h=1-uk;
	end
	vs=v;%reshape(myCentralDiff(reshape(v,m,n),1),m*n,1);
	
	uk0=uk*(db_a);% == u0 or g times constant
	%uk2=uk;
	%%% start iterating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('iterating..');
	deltaTover2=deltaT*0.5;
	rd_2=rhO*deltaTover2;
	dao_r=deltaT*alphA*omegA/rhO;
	dao_r2=dao_r/rhO;
	d_r2=deltaT/rhO^2;
	md=1+floor(niter/11);
	subFigN=1;
	sq=sqrt(0.5);
	figure(1),subplot(3,4,1),imshow(reshape(uk,m,n,1)),title('u');
	for i=1:niter,
		hs=reshape(myCentralDiff(reshape(h,m,n),1),m*n,1);
		h2=h.^2;
		c=1.0./(1.0.+hs*db_a);
		%%% sigmA times diffuse u %%%
		u_temp=uk;
	%	u_temp2=uk2;
		ttl1=num2str(i);%strcat('delta_t',num2str(deltaT),'_b',num2str(betA),'_i',num2str(niter));
		ttl2=ttl1; %%%%%%%%%%%!!!!!!!!
		ttl3=ttl1; %%%%%%%%%%%!!!!!!!!		
		g=(v.^2)*(1-omegA)+(h.^2)*omegA;
	%	g2=vr.^2;
		Lg=L*g;
	%	Lg2=L*g2;
		for j=1:sigmA,
			uk1=(uk0.*hs+uk+(g.*(L*uk)-uk.*Lg+L*(g.*uk))*deltaTover2).*c;
			uk=uk1;
		end
		us=reshape(myCentralDiff(reshape(uk,m,n),1),m*n,1);
	%{
		for j=1:sigmA2,
			uk3=(uk0+uk2+(g2.*(L*uk2)-uk2.*Lg2+L*(g2.*uk2))*deltaTover2)*c;
			uk2=uk3;
		end
	%}		
	%	us2=reshape(myCentralDiff(reshape(uk2,m,n),1),m*n,1);
		%B=rhO^2./(b-us*d);
	%	dB2=1.0./(us*d2+1+deltaT/rhO2^2);%diagonal vec of inv(diag(B))
	%us*d*alphA+hs*omegA*d
	
		dB=1.0./(us*((2*deltaT*alphA*(1-omegA))/rhO)+(1+d_r2));%diagonal vec of inv(diag(B))
	%{
		for j=1:rhO2n,
			vr1=dB2.*(A*vr+deltaT/rhO2^2);
			vr=vr1;
		end
	%}
		%h=1.0./(vr+1.0);%1.0./(1.0+reshape(myCentralDiff(reshape(v.^2,m,n),1),m*n,1)*(alphA*rhO*rhO));%ones(m*n,1);%
		%h=(h-0.5).*2;
		%h=v;
		%reshape(myCentralDiff(reshape(h,m,n),1),m*n,1); # (cv)^2??;
		Lh2=L*h2;
		%%% for rhOn times diffuse v %%%
		for j=1:rhOn,
			vk1=dB.*(v+(h2.*(L*v)-v.*Lh2+L*(h2.*v))*deltaTover2+d_r2);%faster than inverse %%%
			v=vk1;
		end
		
		vs=reshape(myCentralDiff(reshape(v,m,n),1),m*n,1);
		%dH=1.0./(vs*(deltaT*rhO/(2*alphA*omegA))+((uk-u_original).^2)*(deltaT*betA/(alphA*omegA))+1.0);
		dH=1.0./((vs*rhO+us*(2*alphA*omegA))*(deltaT/(2*betA))+1);
		%%% for rhOn times diffuse v %%%
		ug2=(uk-u_original).^2;%v.^2;%.*us;
		Lug2=L*ug2;
		%aov2=v.^2*ao_ao;
		%Laov2=L*aov2;
		for j=1:gammaN,
			h1=dH.*(h+(ug2.*(L*h)-h.*Lug2+L*(ug2.*h))*deltaTover2);%faster than inverse %%%
			h=h1;
		end
		%mean(h)
		
		if mod(i,md)==0 || i == 1
			if subFigN == 1
				%figure(1),subplot(3,4,subFigN),imshow(reshape(uk,m,n,1)),title('u');
				figure(2),subplot(3,4,subFigN),imshow(reshape(v,m,n,1)),title('v');
				figure(7),subplot(3,4,subFigN),imshow(reshape(h,m,n,1)),title('h');				
				figure(8),subplot(3,4,subFigN),imshow(reshape(vs,m,n)),title('vs');
				figure(4),subplot(3,4,subFigN),imshow(histeq(reshape(v,m,n,1))),title('histeqEdge');
				uDiff=abs(u_temp-uk);
				%figure(3),subplot(3,4,subFigN),imshow(reshape(uDiff,m,n)),title('diff');
			else
				%omegA*=1.1;
				disp(subFigN)
				figure(1),subplot(3,4,subFigN),imshow(reshape(uk,m,n,1)),title(ttl1);
				figure(2),subplot(3,4,subFigN),imshow(reshape(v,m,n,1)),title(ttl2);
				figure(7),subplot(3,4,subFigN),imshow(reshape(h,m,n,1)),title(ttl2);
				figure(8),subplot(3,4,subFigN),imshow(reshape(vs,m,n)),title(ttl2);
				figure(4),subplot(3,4,subFigN),imshow(histeq(reshape(v,m,n,1))),title(ttl2);
				uDiff=abs(u_temp-uk);
				%figure(3),subplot(3,4,subFigN),imshow(reshape(uDiff,m,n)),title(ttl3);
			end
			subFigN=subFigN+1;
		end
		%break if converged
		if abs(u_temp-uk) < abs(u_temp)*epsiloN
			disp('converged');
			break;		
		else 
			if i==niter
				epsi=max(max(u_temp))/max(max(u_temp-uk));
				disp(strcat('epsiloN:',num2str(epsi)));
			end
		
		end
	end
	%fill the remaining subplots
	for k=subFigN:12
		figure(1),subplot(3,4,k),imshow(reshape(uk,m,n,1)),title(ttl1);
		figure(2),subplot(3,4,k),imshow(reshape(v,m,n,1)),title(ttl2);
		%figure(3),subplot(3,4,k),imshow(reshape(uDiff,m,n)),title(ttl3);
		figure(4),subplot(3,4,k),imshow(histeq(reshape(v,m,n,1))),title(ttl1);
		figure(7),subplot(3,4,k),imshow(reshape(h,m,n)),title(ttl1);
		figure(8),subplot(3,4,k),imshow(reshape(vs,m,n)),title(ttl2);
	end
	figttl=strcat(imName,'a',num2str(alphA),'b',num2str(betA),'d',num2str(deltaT),'r',num2str(rhO),'n',num2str(niter),'o',num2str(omegA));
	figU=figure(1);	  % to observe u
	saveas(figU,strcat('modification2_',figttl,'U'),'png');
	figV=figure(2);   % to observe v
	saveas(figV,strcat('modification2_',figttl,'V'),'png');
	figH=figure(7);% to observe convergence & what is lost
	saveas(figH,strcat('modification2_',figttl,'H'),'png');
	figure(5),imshow(reshape(uk,m,n))
	%figure(6),imshow(reshape(1.0./(1.0+((1-v).^2)*(alphA*rhO)),m,n));
	figure(100),imshow(histeq(reshape(uk,m,n)));
	figUhisteq=figure(100);
	saveas(figV,strcat('modification2_',figttl,'U_histeq'),'png');
	figure(101),imshow(histeq(reshape(h,m,n)));
	figHhisteq=figure(101);
	saveas(figV,strcat('modification2_',figttl,'H_histeq'),'png');
	figure(102),imshow(histeq(reshape(v,m,n)));
	figHhisteq=figure(102);
	saveas(figV,strcat('modification2_',figttl,'V_histeq'),'png');
end
