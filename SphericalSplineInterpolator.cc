/*
 * =====================================================================================
 *
 *       Filename:  SphericalSplineInterpolator.cc
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

#include "SphericalSplineInterpolator.hh"
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

int myrandom (int i) 
{ 
    return rand()%i;
}

SphericalSplineInterpolator::SphericalSplineInterpolator(float *lons, float *lats, float *vals, int n)
:m_numNodes(n)
{
    /*-----------------------------------------------------------------------------
     * Resize arrays
     *-----------------------------------------------------------------------------*/
    m_lons.resize(m_numNodes);
    m_lats.resize(m_numNodes);
    m_vals.resize(m_numNodes);
    m_xs.resize(m_numNodes);
    m_ys.resize(m_numNodes);
    m_zs.resize(m_numNodes);
    
    vector<int> order(m_numNodes);
    for(int i=0; i<m_numNodes; i++)
    {
        m_lons[i] = lons[i];
        m_lats[i] = lats[i];
        m_vals[i] = vals[i];

        order[i]=i;
    }
    
    /*-----------------------------------------------------------------------------
     * Shuffle data to ensure the first 3 points are not colinear
     *-----------------------------------------------------------------------------*/
    random_shuffle(order.begin(), order.end(), myrandom);

    for(int i=0; i<m_numNodes; i++) 
    {
        float lat = m_lats[i];
        m_lats[i] = m_lats[order[i]];
        m_lats[order[i]] = lat;
        
        float lon = m_lons[i];
        m_lons[i] = m_lons[order[i]];
        m_lons[order[i]] = lon;
        
        float val = m_vals[i];
        m_vals[i] = m_vals[order[i]];
        m_vals[order[i]] = val;
    }
   
    trans_( &m_numNodes, 
            m_lats.data(), 
            m_lons.data(), 
            m_xs.data(), 
            m_ys.data(), 
            m_zs.data() );

    /*-----------------------------------------------------------------------------
     * Filter out nodes close to the poles
     *-----------------------------------------------------------------------------*/
    {
        vector<float> fx;
        vector<float> fy;
        vector<float> fz;
        vector<float> fv;
        
        int counter = 0;
        for(int i=0; i<m_numNodes; i++)
        {
            if(fabs(1-fabs(m_zs[i])) < 1e-5){}
            else
            {
                fx.push_back(m_xs[i]);
                fy.push_back(m_ys[i]);
                fz.push_back(m_zs[i]);
                fv.push_back(m_vals[i]);

                counter++;
            }
        }

        m_numNodes = counter;
        m_xs.resize(m_numNodes);
        m_ys.resize(m_numNodes);
        m_zs.resize(m_numNodes);
        m_vals.resize(m_numNodes);

        copy(fx.begin(), fx.end(), m_xs.begin());
        copy(fy.begin(), fy.end(), m_ys.begin());
        copy(fz.begin(), fz.end(), m_zs.begin());
        copy(fv.begin(), fv.end(), m_vals.begin());
    }

    /*-----------------------------------------------------------------------------
     * Create triangulation 
     *-----------------------------------------------------------------------------*/
    m_list.resize(m_numNodes*6);
    m_lptr.resize(m_numNodes*6);
    m_lend.resize(m_numNodes);
    m_near.resize(m_numNodes);
    m_next.resize(m_numNodes);
    m_dist.resize(m_numNodes);
     
    int lnew  = 0;
    int error = 0;
    trmesh_( &m_numNodes, 
             m_xs.data(), 
             m_ys.data(), 
             m_zs.data(), 
             m_list.data(), 
             m_lptr.data(), 
             m_lend.data(), 
             &lnew, 
             m_near.data(), 
             m_next.data(), 
             m_dist.data(), 
             &error );

    cout << "Triangulation: " << error << endl;
}

SphericalSplineInterpolator::~SphericalSplineInterpolator()
{

}

void SphericalSplineInterpolator::InterpolateC0(int n, float *queryLons, float *queryLats, float *results)
{
    int ist = 1; /* start from first triangle vertex */
    int ier = 0; /* error reporter */
    
    int ierSum = 0;
    for(int i=0; i<n; i++)
    {
        float plon = queryLons[i];
        float plat = queryLats[i];
        
        intrc0_( &m_numNodes, 
                 &plat, 
                 &plon, 
                 m_xs.data(), 
                 m_ys.data(), 
                 m_zs.data(), 
                 m_vals.data(), 
                 m_list.data(), 
                 m_lptr.data(), 
                 m_lend.data(), 
                 &ist, 
                 &results[i], 
                 &ier );
        ierSum += ier;
    }

    cout << "Interpolation: " << ierSum/n << endl;
}

