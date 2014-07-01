3D_Retrieval_scan2search
========================

Spherical harmonics descriptor


##Functions
Pseudo code for spherical harmonics:

    dff
    fff

Function: double gsl_sf_legendre_sphPlm (int l, int m, double x)  
Function: int gsl_sf_legendre_sphPlm_e (int l, int m, double x, gsl_sf_result * result)  
These routines compute the normalized associated Legendre polynomial \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x) suitable for use in spherical harmonics. The parameters must satisfy m >= 0, l >= m, |x| <= 1. Theses routines avoid the overflows that occur for the standard normalization of P_l^m(x).

The Legendre polynomial P(n,x) can be defined by:

    P(0,x) = 1
    P(1,x) = x
    P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)


##Reference:
1. Spherical harmonics and Legendre polynomials ,involving solution when m is negative:  
http://blog.sciencenet.cn/blog-548663-715825.html  
2. Rodrigues' formula: equation dn/dxn = (d/dx)n:  
http://wenku.baidu.com/view/72c151ef102de2bd960588f2  
3. Rotation Invariant Spherical Harmonic of 3D Shape:  
http://www.chenkuantong.com/?p=1210  
4. Spherical harmonics library:  
http://www.cs.dartmouth.edu/~geelong/sphere/  
5. GSL library to calculate Spherical harmonics:  
https://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html
