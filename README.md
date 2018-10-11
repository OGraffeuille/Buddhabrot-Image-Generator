# Buddhabrot
A c++ program to generate images using the "Buddhabrot" algorithm, a modified version of the Mandelbrot pattern. This code is a new version of a Matlab project I worked on a number of years ago, in a (eventually successful) attempt to make the algorithm more efficient to be able to make print quality images.

Print quality images are generally ~7500x7500 bmp files. A few lower resolution (2500x2500) images are shown here:

![buddhabrot example image 1](https://drive.google.com/uc?id=1GZLf3f8yhi0Rbg2gmroRac-YHl6Y6JIz)
![buddhabrot example image 2](https://drive.google.com/uc?id=1MPGPFIHTT1p2shx1ah1gV5pk_tfKaN8F)
![buddhabrot example image 2](https://drive.google.com/uc?id=1Rg3zEAskgumz0M7coD6GvwIkpRz_O85D)


Note:
'Lopsided' Buddhabrot images are made by changing the range of possible starting values for the first iteration of the algorithm, such that both real and imaginary components are positive.
