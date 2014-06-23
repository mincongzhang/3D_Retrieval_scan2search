3D_Retrieval_scan2search
========================

Spherical harmonics descriptor

Function: double gsl_sf_legendre_Plm (int l, int m, double x)
Function: int gsl_sf_legendre_Plm_e (int l, int m, double x, gsl_sf_result * result)

The Legendre polynomial P(n,x) can be defined by:

    P(0,x) = 1
    P(1,x) = x
    P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)

reference:  
1. Spherical harmonics and Legendre polynomials ,involving solution when m is negative:  
http://blog.sciencenet.cn/blog-548663-715825.html  
2. Rodrigues' formula: equation dn/dxn = (d/dx)n:  
http://wenku.baidu.com/view/72c151ef102de2bd960588f2  
3. Rotation Invariant Spherical Harmonic of 3D Shape:  
http://www.chenkuantong.com/?p=1210  
4. Spherical harmonics library:  
http://www.cs.dartmouth.edu/~geelong/sphere/  
