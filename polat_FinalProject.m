% CENG 566 Final Project implemented by A. Goze Polat, 1631092
% The code was tested with octave v3.6.4, it may or may not work on matlab
% image package is necessary for this code to run
pkg load all

%disp('for perona malik experiment, see polat_hw3.m')
% or enter polat_hw3 folder and run:
% >> polat_hw3
% or
% >> I=imread('NoiseVelopMinix2.png');
% >> myCLMC(I+1,1,500,0.05,0.02,1); 
% note that this folder has some older versions of various modules such as myA

%%%%% DEFAULTS %%%%%%%%%
imName='noisevelop.png';
%imName='tigerling.png';
noiseLvl=0;
epsiloN=0.0005;
niter=500;
deltaT=0.25;
%%%%%%%%%%%%%%%%%%%%%%%%%
disp('*****************************************************************');
disp('Final Project by A. Goze Polat, for Ceng 566 Image Processing')
disp('Original AT approximation running..');
disp('*****************************************************************');
alphA=5000;
betA=200;
rhO=0.05;
myAT('noisevelop.png',noiseLvl,epsiloN,niter,deltaT,alphA,betA,rhO);
disp('*****************************************************************');
disp('First modification to AT: add diffusion speed (temperature h) to edge set. The speed of heat diffusion for h depends on the error between u and g. The edge strengthening example is running.. (Observe the figures for h,v and u.)');
disp('*****************************************************************');
alphA=4000;
betA=1;
rhO=0.08;
myHeatModifiedAT(imName,noiseLvl,epsiloN,2000,deltaT,alphA,betA,rhO,1,1);
disp('*****************************************************************');
disp('First modification to AT: add diffusion speed (temperature h) to edge set. The speed of heat diffusion for h depends on the error between u and g. The texture removal example is running..');
disp('*****************************************************************');
alphA=0.005;
betA=1;
rhO=0.5;
disp('*****************************************************************');
myHeatModifiedAT('tigerling.png',noiseLvl,epsiloN,2000,deltaT,alphA,betA,rhO,1,1);
alphA=10000.0;
betA=5000.0;
rhO=0.2;
omegA=0.2;
disp('Second modification to AT: Heat edge feedback to u, to get produce a more "curvy" image. h initialized as g. The model is running..');
disp('*****************************************************************');
% to make edge strengthening on h more visible, increase alpha,e.g. myFinalAT('noisevelop.png',0,0.0005,1500,0.25,50000,10000,0.05,0.1,1)
myFinalAT('noisevelop.png',noiseLvl,0.00001,500,deltaT,alphA,betA,rhO,omegA,1);
disp('Second modification to AT: Heat edge feedback to u, to get produce a more "curvy" image. h initialized as 1-g. The model is running..');
disp('*****************************************************************');
myFinalAT('noisevelop.png',noiseLvl,0.00001,500,deltaT,alphA,betA,rhO,omegA,-1);
disp('*****************************************************************');
disp('End of examples..');
disp('*****************************************************************');
disp('Thanks for the great course!');
disp('And special thanks for always sharing your scientific wisdom with us in the class :)');
