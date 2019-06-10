# MATLAB Code for the Research in Graduate School

## What is Synthetic-aperture radar (SAR)?

Synthetic-aperture radar (SAR) is a radar that create two-dimensional images of landscapes under any weather circumstance.

And, what is polarimetric SAR (PolSAR)?
Note that objects reflect radar waves with different intensities with respect to the polarization of the incident waves.
So, we can obtain the information of the objects by analyzing the characteristic of reflect waves.

SAR Wiki: https://en.wikipedia.org/wiki/Synthetic-aperture_radar

PolSAR Wiki: https://en.wikipedia.org/wiki/Polarimetry

## Folder - VelocitySAR

The first research topic in graduate school. The research is aimed to detect the velocity of moving target in the radar image.
It is known that the moving target makes images blur in synthetic aperture radar (SAR) image. To fix blurring in image, we need the information of velocity of the target.

## Folder - PolSAR

##

The source code of radar image toolbox for polarimetric synthetic aperture radar (PolSAR) was including in this folder.
They are inside the folders 
[PolSAR/src/@PolSAR_AnalyzeTool](https://github.com/eesaber/Matlab/tree/master/PolSAR/src/%40PolSAR_AnalyzeTool)
and [PolSAR/src/@SeaIce](https://github.com/eesaber/Matlab/tree/master/PolSAR/src/%40SeaIce), respectively.
The toolbox was implemented to generate the experiment results for the master thesis.

### PolSAR_AnalyzeTool

`PolSAR_AnalyzeTool`is the class for analyzing the characteristic of PolSAR images. 
It wraps the radar images, classification results as data member, and implements the member function for analyzing PolSAR images.

1. The conventional radar image analysis algorithms were implemented as class member function, including `eigenDecomposition()`, `fourComponentDecomposition()` and `pauliDecomposition()`. 
The conventional radar image analysis algorithms were used to classified the objects in the image.

2. We can assess polarimetric characteristic of the object in the image by using the member functions `para<name>()`, 
including `paraCoPolCorrelation()`, `paraMoisture()` and `paraRoughness1()`,  etc. 
For example, we can find the level of roughness of the ground by calling `paraRoughness1()`.

3. The machine learning algorithms, including support vector machine (SVM), fuzzy *c*-means (FCM) and *k*-means were wrapped in the member functions, `mySVM()`, `myFCM()` and `myKmeans()`, respectively. 
By calling these function, the supervised model, SVM, will trained on image and then show the classified results; 
the unsupervised model, FCM and *k*-means will show the classified results.

The deep learning model, convolutional neural network (CNN), is implemented by python, thus, in [other repository](https://github.com/eesaber/PolSAR_ML).

#### How to Use
---------------
* Initialization a `PolSAR_AnalyzeTool` class object as 
```Matlab
img_a = PolSAR_AnalyzeTool()
```

* To apply "four-component decomposition", a kind of radar image analysis, on the radar image stored in `img_a`, 
call the member function `fourComponentDecomposition()` as
```Matlab
img_a.fourComponentDecomposition()
```

* To find the usage of a member function, type 
```console
> help img_a.fourComponentDecomposition
```
in the console. Then, the usage of `fourComponentDecomposition` is shown as 
```console
> help img_a.fourComponentDecomposition
> --- help for PolSAR_AnalyzeTool/fourComponentDecomposition ---

  fourComponentDecomposition implements the original 4-component decomposition.
 
  syntax:
  fourComponentDecomposition(varargin)
  Specifies function properties using one or more Name,Value pair arguments. 
  Name-value pair settings apply to all the lines plotted. 
 
  Inputs:
  Name-Value Pair Arguments:
     ('vv', s) - The pixel above SNR threshold is marked. If s is ture, 
                 function will show the figure.
     ('hh', s) - The pixel above SNR threshold is marked. If s is ture, 
                 function will show the figure.
     ('hv', s) - The pixel above SNR threshold is marked. If s is ture, 
                 function will show the figure.
 
  Example: 
     1.
 
  Other m-files required: plot_para
  Subfunctions: none
  MAT-files required: none
 
 ------------- <<<<< >>>>>--------------
  Author: K.S. Yang
  email: fbookzone@gmail.com
 ------------- <<<<< >>>>>--------------
```

### SeaIce
`SeaIce` was the derived class from `PolSAR_AnalyzeTool`. This class was designed for sea-ice classification. The usage is the same as
`PolSAR_AnalyzeTool`, and the member functions, which were specific for image classification, were implemented.

#### How To Use
-------------
* Initialization a class object as
```Matlab
img_s = SeaIce()
```

## Folder - Homework
The MATLAB code used in course homework.
