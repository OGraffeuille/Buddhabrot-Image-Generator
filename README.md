# Buddhabrot
This repository contains a c++ program which generates images using the "Buddhabrot" algorithm. This algorithm is a modified version of the well known Mandelbrot fractal. More information on this algorithm can be found at https://en.wikipedia.org/wiki/Buddhabrot.

The Buddhabrot is created by sampling random points on the complex plane and iterating them using the Mandelbrot sequence, and plotting the path of points which diverge after a specified range of iterations. Specifying different number of iterations to plot will change the appearance of the Buddhabrot. It's worth noting that since this is an inherently random process, there is no 'true' buddhabrot, just generated instances of it. 

This code is a new version of a Matlab project I completed years ago. Re-writing this in C++ allowed me to make the algorithm more efficient, and allowed me to generate images of much higher resolution.

Print quality images are generally ~7500x7500 bmp files. A few lower resolution (2500x2500) images are shown here:

![buddhabrot example image 1](https://drive.google.com/uc?id=1GZLf3f8yhi0Rbg2gmroRac-YHl6Y6JIz)
![buddhabrot example image 2](https://drive.google.com/uc?id=1MPGPFIHTT1p2shx1ah1gV5pk_tfKaN8F)
![buddhabrot example image 2](https://drive.google.com/uc?id=1Rg3zEAskgumz0M7coD6GvwIkpRz_O85D)

Note:
'Lopsided' Buddhabrot images are made by changing the range of possible starting values for the first iteration of the algorithm, such that both real and imaginary components are positive.
