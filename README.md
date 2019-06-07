# MATLAB Code for the Research in Graduate School

## Folder - VelocitySAR

The first research topic in graduate school. The research is aimed to detect the velocity of moving target in the radar image.
It is known that the moving target makes images blur in synthetic aperture radar (SAR) image. To fix blurring in image, we need the information of velocity of the target.

## Folder - PolSAR

The source code of radar image toolbox for polarimetric synthetic aperture radar (PolSAR) is including in this folder.
The toolbox was implemented to generate the experiment results for the master thesis.

### PolSAR_AnalyzeTool

`PolSAR_AnalyzeTool`is the class for PolSAR analysis. Algorithms for radar image analysis were impelmented as memeber function, including `eigenDecomposition()`, `fourComponentDecomposition()` and `pauliDecomposition()`. 
And, memeber functions `para<name>()`, including `paraCoPolCorrelation()`, `paraMoisture()` and `paraRoughness1()`,  etc. 
were used to extract the characteristic of each pixel in radar image.

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
     ('vv', s) - The pixel above SNR thershold is marked. If s is ture, 
                 function will show the figure.
     ('hh', s) - The pixel above SNR thershold is marked. If s is ture, 
                 function will show the figure.
     ('hv', s) - The pixel above SNR thershold is marked. If s is ture, 
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

How To Use
-------------
* Initialization a class object as
```Matlab
img_s = SeaIce()
```

## Folder - Homework
The MATLAB code used in course homework.
