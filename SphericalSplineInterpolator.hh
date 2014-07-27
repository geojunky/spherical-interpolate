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
int intrc0_(integer *n, doublereal *plat, doublereal *plon, doublereal *x, 
    doublereal *y, doublereal *z__, doublereal *w, integer *list, integer *lptr, integer *
    lend, integer *ist, doublereal *pw, integer *ier);

int intrc1_(integer *n, doublereal *plat, doublereal *plon, doublereal *x, 
    doublereal *y, doublereal *z__, doublereal *f, integer *list, integer *lptr, integer *
    lend, integer *iflgs, doublereal *sigma, integer *iflgg, doublereal *grad, 
    integer *ist, doublereal *fp, integer *ier);

int trmesh_(integer *n, doublereal *x, doublereal *y, doublereal *z__, integer 
    *list, integer *lptr, integer *lend, integer *lnew, integer *near__, 
    integer *next, doublereal *dist, integer *ier);

int trans_(integer *n, doublereal *rlat, doublereal *rlon, doublereal *x, doublereal 
    *y, doublereal *z__);

int gradg_(integer *n, doublereal *x, doublereal *y, doublereal *z__, doublereal *f, 
    integer *list, integer *lptr, integer *lend, integer *iflgs, doublereal *
    sigma, integer *nit, doublereal *dgmax, doublereal *grad, integer *ier);

int getsig_(integer *n, doublereal *x, doublereal *y, doublereal *z__, doublereal *
    h__, integer *list, integer *lptr, integer *lend, doublereal *grad, doublereal *
    tol, doublereal *sigma, doublereal *dsmax, integer *ier);
#ifdef __cplusplus
}
#endif

#include <vector>
#include <set>

using namespace std;

class SphericalSplineInterpolator
{
    public:
    
    SphericalSplineInterpolator(double *lons, double *lats, double *vals, int n);
    
    void InterpolateC0(int n, double *queryLons, double *queryLats, double *results);
    void InterpolateC1(int n, double *queryLons, double *queryLats, double *results);
    
    ~SphericalSplineInterpolator();

    private:
    int m_numNodes;
    vector<double> m_lons;
    vector<double> m_lats;
    vector<double> m_vals;
    vector<double> m_xs;
    vector<double> m_ys;
    vector<double> m_zs;

    /*-----------------------------------------------------------------------------
     * Triangulation-related arrays 
     *-----------------------------------------------------------------------------*/
    vector<int> m_list;
    vector<int> m_lptr;
    vector<int> m_lend;
    vector<int> m_near;
    vector<int> m_next;
    vector<double> m_dist;
};

#endif

