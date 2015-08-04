CENG 566 Final Project implemented by A. Goze Polat

*********************************************************
PROJECT REPORT AND POSTER

	polat_ceng566_FinalProject.pdf
	polat_ceng566_FinalPoster.pdf
	
*********************************************************
IMAGE REGULARIZATION
*********************************************************
A common inverse problem of image processing is regularization of images. As any regularization problem, the 
solution requires one to introduce apriori information. For example by assuming that high 
variation in an image is not natural, denoising and texture removal can be achieved. To do that, one needs to use 
a formal representation that can incorporate such assumptions. For this purpose, energy functionals are used.
High variation and other undesired properties increase the energy. Thus, the whole regularization task 
can be considered as a minimization problem. 

Although images can have high variation, this is often undesired since the computational resources to process 
and classify an image and the amount of information to represent it increase. Therefore,  reducing the amount of 
variation within an image and reducing the edge count in a meaningful way (i.e. acquiring a cartoon image = a 
simpler and more compact representation of the image) is a prominent interest area for image processing.

In this project, new sets of partial differential equations are introduced that modify and improve Ambrosio
Tortelli approximation of Mumford Shah formula to achieve better image regularization.
*********************************************************
IMPLEMENTATION 

Implementation of the modification 1 idea is available in:
	myHeatModifiedAT.m
	
Implementation of the modification 2(&1) idea is available in:
	myFinalAT.m
	
For the examples and other details regarding the implementation:
	polat_FinalProject.m
	
*********************************************************
RUNNING THE EXAMPLES

After installing octave and necessary packages, 

Enter into the folder where polat_FinalProject.m exists from the command prompt and then start octave by running

octave

(I always use "octave --traditional" when I start octave, to make it more compatible with matlab, however this may not be helpful.)

Or start octave from the Applications and when the octave command prompt opens, enter into the right folder by

 cd some_path/to/Image_Regularization

You can run the examples by executing the script
 polat_FinalProject

Note that there will be a lot of figures when all the examples are finished. To close all the figures at once, simply exit octave by typing
>> quit

*********************************************************
INSTALLATION USING PACKAGE MANAGERS

Install XCode via the Mac App Store. Once installed, install the Command Line Tools from XCode's Apple Menu > Preferences > Downloads.

INSTALLATION of FINK PACKAGE MANAGER

Follow Fink's installation instructions (http://www.finkproject.org/download/srcdist.php). For those who prefer it, there is a GUI available for Fink, Fink Commander.

To install the latest octave, type 
 sudo fink install octave-atlas 
at the Terminal's command line. (for 64 bit version of Fink)
If you use 32 bit version, type
 sudo fink install octave

INSTALL PACKAGES USING FINK:
    Fink should also be used to install packages for Octave. For example, the control systems package for Lion may be installed by typing 
  fink install control-atlas-oct362 
at terminal. 

If this does not work, try to install packages from octave command prompt by running:

 pkg install -forge general
 pkg install -forge control
 pkg install -forge specfun
 pkg install -forge signal
 pkg install -forge image
 
********************************************************************
Other alternative package managers such as Homebrew and Macports can also be used for installing octave. (See http://wiki.octave.org/Octave_for_MacOS_X#Binary_installer_for_OSX_10.9.1)

********************************************************************
MANUAL INSTALLATION FROM BINARY FILES

A. INSTALLATION of OCTAVE

1. To install Octave on your Mac, begin by downloading the latest binary release from http://sourceforge.net/projects/octave/files/Octave%20MacOSX%20Binary/2013-12-30%20binary%20installer%20of%20Octave%203.8.0%20for%20OSX%2010.9.1%20%28beta%29/

	(A more stable but older version is availabe at http://sourceforge.net/projects/octave/files/Octave%20MacOSX%20Binary/2011-04-21%20binary%20of%20Octave%203.4.0/)

2. In the Finder, double-click on the downloaded GNU_Octave_3.8.0-6.dmg file in your "Downloads" folder, and wait for the disk image to mount.

3. Once the disk image has mounted, drag the Octave application, which is located inside, to your "Applications" folder. 

B. INSTALLATION of GNUPLOT

1. The Gnuplot application is available in the "Extras" folder on the disk image that you downloaded to install the Octave application.

2. Double-click on the gnuplot.dmg file located in the "Extras" folder of the disk image.

3. When the disk image has mounted drag the Gnuplot application to your "Applications" folder. 

C. INSTALLATION of  AQUATERM

1. Download aquaterm from http://sourceforge.net/projects/aquaterm/files/

2. Double-click to mount the downloaded disk image.

3. Inside the mounted disk image, double click on the AquaTerm.pkg file to install AquaTerm. Agree to the license and follow the prompts to install AquaTerm in your "Applications" folder. 

Restart your computer after the installation process is complete. This will help prevent some of the most common problems that you might encounter after installing Octave. 

D: INSTALLATION of PACKAGES

1. Double-click the Octave application located in your "Applications" folder.

After a few seconds the Mac Terminal application should open and display an Octave command prompt.

3. run the commands:
	pkg install -forge general
	pkg install -forge control
	pkg install -forge specfun
	pkg install -forge signal
	pkg install -forge image

***************************************************************
If there is any problem with the installation, alternative ways to install octave may exist. See http://wiki.octave.org/Octave_for_MacOS_X#Binary_installer_for_OSX_10.9.1 
for further detail.
***************************************************************
