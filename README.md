3D_Retrieval_scan2search
========================

Sometimes stack overflow problem occurs. Maybe try not to use array larger than 100 elements.  
溢出问题解决  
改为堆变量?或者实验BitSet
kd-tree can be used to fast retrive in the database

##Structure
1. Pre-processing:
Normlize  
Rasterize  
Denoise  

2. Descriptor:
Spherical harmonics descriptor  

3. Retrieve:
quick-sort?/kd-tree?  
Show candidate(can be more) (clustering???)  

4. Verification:
Denoise-test: Bilateral filter  
Rotate-test: SH rotate invariant  



##Spherical harmonics descriptor 

Pseudo code for spherical harmonics:

    sort vertex according to radius
    for each frequency l (idx_l)
        for each rasterized vertex in one radius region
            calculate F_lr = F(idx_l,idx_r);
        get sum of F_lr in one radius ragion = a_ml
        spherical harmonics descriptor SH(idx_l,idx_r) += abs(a_ml)^2;
    SH = sqrt(SH);
        
*where F_lr is the following equation, l is m in the equation
![SH](https://github.com/mincongzhang/3D_Retrieval_scan2search/raw/master/spherical harmonics.jpg)    


Function: double gsl_sf_legendre_sphPlm (int l, int m, double x)  
Function: int gsl_sf_legendre_sphPlm_e (int l, int m, double x, gsl_sf_result * result)  
These routines compute the normalized associated Legendre polynomial \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x) suitable for use in spherical harmonics. The parameters must satisfy m >= 0, l >= m, |x| <= 1. Theses routines avoid the overflows that occur for the standard normalization of P_l^m(x).

The Legendre polynomial P(n,x) can be defined by:

    P(0,x) = 1
    P(1,x) = x
    P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)


##notes
computing SH very slow: precompute?   
(1)divide and concour  http://www.ams.org/journals/mcom/2002-71-238/S0025-5718-01-01386-2/  
(2)FFT to fourier then SH:http://connection.ebscohost.com/c/articles/67655125/3d-objects-retrieval-using-spherical-harmonics-feature-vector-method  
(3)face based? http://liris.cnrs.fr/Documents/Liris-2276.pdf  

##Reference:
1. Spherical harmonics and Legendre polynomials ,involving solution when m is negative:  
http://blog.sciencenet.cn/blog-548663-715825.html  
2. Rodrigues' formula: equation dn/dxn = (d/dx)n:  
http://wenku.baidu.com/view/72c151ef102de2bd960588f2  
3. Rotation Invariant Spherical Harmonic of 3D Shape:  
http://www.chenkuantong.com/?p=1210  
4. An explaination of a_mn component(arbitrary constants?):  
https://www.ngs.noaa.gov/PUBS_LIB/Geodesy4Layman/TR80003F.HTM
5. Spherical harmonics library:  
http://www.cs.dartmouth.edu/~geelong/sphere/  
6. Spherical Harmonics Visual Representation:  
http://afj-phd.blogspot.co.uk/2008/11/spherical-harmonics-visual.html  
7. GSL library to calculate Spherical harmonics:  
https://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html  
8. Shape Descriptors from John Hopkins  
http://www.cs.jhu.edu/~misha/Code/ShapeSPH/ShapeDescriptor/
