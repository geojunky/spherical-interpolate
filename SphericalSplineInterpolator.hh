/*
 * =====================================================================================
 *
 *       Filename:  SphericalSplineInterpolator.hh
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/07/14 11:22:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rakib Hassan (rakib.hassan@sydney.edu.au), 
 *        Company:  The University of Sydney
 *
 * =====================================================================================
 */

#ifndef SPHERICAL_SPLINE_INTERPOLATOR
#define SPHERICAL_SPLINE_INTERPOLATOR

#ifdef __cplusplus
extern "C" {
#endif

#include <f2c.h>
int intrc0_(integer *n, real *plat, real *plon, real *x, 
    real *y, real *z__, real *w, integer *list, integer *lptr, integer *
    lend, integer *ist, real *pw, integer *ier);

int intrc1_(integer *n, real *plat, real *plon, real *x, 
    real *y, real *z__, real *f, integer *list, integer *lptr, integer *
    lend, integer *iflgs, real *sigma, integer *iflgg, real *grad, 
    integer *ist, real *fp, integer *ier);

int trmesh_(integer *n, real *x, real *y, real *z__, integer 
    *list, integer *lptr, integer *lend, integer *lnew, integer *near__, 
    integer *next, real *dist, integer *ier);

int trans_(integer *n, real *rlat, real *rlon, real *x, real 
    *y, real *z__);

int gradg_(integer *n, real *x, real *y, real *z__, real *f, 
    integer *list, integer *lptr, integer *lend, integer *iflgs, real *
    sigma, integer *nit, real *dgmax, real *grad, integer *ier);

int getsig_(integer *n, real *x, real *y, real *z__, real *
    h__, integer *list, integer *lptr, integer *lend, real *grad, real *
    tol, real *sigma, real *dsmax, integer *ier);
#ifdef __cplusplus
}
#endif

#include <vector>
#include <set>

using namespace std;

class SphericalSplineInterpolator
{
    public:
    
    SphericalSplineInterpolator(float *lons, float *lats, float *vals, int n);
    
    void InterpolateC0(int n, float *queryLons, float *queryLats, float *results);
    
    ~SphericalSplineInterpolator();

    private:
    int m_numNodes;
    vector<float> m_lons;
    vector<float> m_lats;
    vector<float> m_vals;
    vector<float> m_xs;
    vector<float> m_ys;
    vector<float> m_zs;

    /*-----------------------------------------------------------------------------
     * Triangulation-related arrays 
     *-----------------------------------------------------------------------------*/
    vector<int> m_list;
    vector<int> m_lptr;
    vector<int> m_lend;
    vector<int> m_near;
    vector<int> m_next;
    vector<float> m_dist;
};

#endif

