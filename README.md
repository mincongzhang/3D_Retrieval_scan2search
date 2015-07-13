3D Shape Retrieval "Scan to search"
========================

##Demo

[Video Demo](https://www.youtube.com/watch?v=IC5RN_ZZ9Zg)

Input

<a style="float:left;">
    <img src="https://github.com/mincongzhang/3D_Retrieval_scan2search/raw/master/Thesis/input_engine_scanned_scantosearch_test.jpg" alt="Scan" title="Scan" height="200"/>
</a>
    <a style="float:left;">
    <img src="https://github.com/mincongzhang/3D_Retrieval_scan2search/raw/master/Thesis/input_engine_scantosearch_test.jpg" alt="Input" title="Input" height="200"/>
</a>

Output

<img src="https://github.com/mincongzhang/3D_Retrieval_scan2search/raw/master/Thesis/output_engine_scantosearch_test.jpg" alt="Output" title="Output" width="600" align="middle"/>

##Main algorithm - Spherical Harmonics
Spherical Harmonics Basic Idea:  
(1)Fouries transform can describe the distribution of a signal, and how the signal changes along time/space domain  
(2)Spherical Harmonics can describe the distribution of a "signal" in space(sphere with certain radius), and how the "signal" changes along space domain  


##Structure
1. Pre-processing:  
(1)Normalize  
-make the center of mass of the model be at the point (R,R,R)  
-scale so that the average distance from vertices to the center of mass is R/2  
-(scale invariance and against outliers)  
(2)Denoise  
-use 3D bilateral filter to denoise  
-(scanned noise invariance)  
(3)Rasterize  
-rasterize in to a 2R*2R*2R voxel grid (normally choose R to be ~32)  
-(provide adequate granularity for discriminating shapes while filtering out high-freq noise)

2. Descriptor:  
Spherical harmonics descriptors(SH)  
Distance histogram descriptors(DH)  

3. Retrieval: 
SH and DH database  
Retrieval with same weight (with 0.5 on SH and DH respectively)

##Rasterization
O(n) solution:

    /* fill voxel grid and get rasterized points in O(n) */
    bitmap = new bitmap(2R*2R*2R)               //temporary bitmap
    for each voxel to be filled
        if(bitmap(voxel.x,voxel.y,voxel.z))==0  //has not been filled or registered
            bitmap(voxel.x,voxel.y,voxel.z)=1
            grid_point.push_back(voxel)
        end
    end
    delete bitmap;

##Spherical harmonics descriptor 

Pseudo code for spherical harmonics:

    //SH is spherical harmonics descriptor 
    sort vertex according to radius
    for each frequency l
        for each rasterized vertex in one radius range r 
            calculate F_lr = F(theta,phi)
        end
        a_ml = sum of F_lr in one radius ragion
        SH(l,r) += abs(a_ml)^2    
    end
    SH = sqrt(SH)
        
*where F_lr is the following equation, l is m in the equation
![SH](https://github.com/mincongzhang/3D_Retrieval_scan2search/raw/master/spherical harmonics.jpg)    


Function: double gsl_sf_legendre_sphPlm (int l, int m, double x)  
Function: int gsl_sf_legendre_sphPlm_e (int l, int m, double x, gsl_sf_result * result)  
These routines compute the normalized associated Legendre polynomial \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x) suitable for use in spherical harmonics. The parameters must satisfy m >= 0, l >= m, |x| <= 1. Theses routines avoid the overflows that occur for the standard normalization of P_l^m(x).

The Legendre polynomial P(n,x) can be defined by:

    P(0,x) = 1
    P(1,x) = x
    P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)

<!--
##Interesting descovery
1. Bilateral filter  
When implementing the bilateral filter I find that it is really slow when I try to get one-ring neighbours for every vertex. Then I realize that the vertex is stored as list, and so the complexity for one-ring neighbours is O(n^2). Thus I build a kd-tree to get some nearest neighbours, and filter out some outliers to approximately get the one-ring neighbours for further denoising. With the kd-tree the complexity can be reduced to O(n)+O(nlog(n)).
-->

##Further improvement
1. kd-tree can be used to fast retrive in the database
 
2. database clustering
 
3. rasterization algorithm updates  

<!--
(idea: 目前的就对每一个面填充3个点,但是有可能有些面有4个点, 所以会出现空出一个三角形的情况)  
(solution: 选出第一个点, 将其按顺序和之后的点俩俩组合成三角形再填充,能保证整个面都填满)


<!--
distance histogram normalization
devide by the total voxels number
:已修改, 但是对10号data有个异常值很奇怪

发现spherical harmonics不normal的话, 
scanned的SH会明显大很多, 但总体趋势差不多,
是不是normalize后会更好一些?
值得一试(如果3DI公司给我offer的话才干)
-->

4. Spherical Harmonics transform speed up  
(1)divide and concour  
http://www.ams.org/journals/mcom/2002-71-238/S0025-5718-01-01386-2/  
(2)FFT to fourier then SH:  
http://connection.ebscohost.com/c/articles/67655125/3d-objects-retrieval-using-spherical-harmonics-feature-vector-method  
(3)O(n) solution for spherical harmonics:  
http://liris.cnrs.fr/Documents/Liris-2276.pdf   

<!--
按照paper"A 3D search engine"里, 它的频率上限就设到了R/2=16, 所以可以速度快一些
不改变原有框架的情况下把频率上限调低一些能快很多, 因为高频循环量非常大
-->

5. New descriptors

<!--
新想法
3D histogram: 3D histogram因为受rotation影响, 不妨以最大的bin为基准对齐?校准?
搜索Shape Histograms的paper
--> 

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
