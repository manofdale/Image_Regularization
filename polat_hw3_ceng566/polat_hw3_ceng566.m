%A. Goze Polat 1631092
%Modified slightly myCLMC and myCentralDiff
pkg load all
I=imread('NoiseVelopMinix2.png');
myCLMC(I+1,1,500,0.05,0.02,1);
