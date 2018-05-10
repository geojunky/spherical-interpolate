/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/07/14 10:59:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Rakib Hassan (rakib.hassan@sydney.edu.au), 
 *        Company:  The University of Sydney
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SphericalSplineInterpolator.hh"
#include <iostream>

using namespace std;
const double PI = 3.14159265359;

void ReadGrid(string fileName, int *nNodes, double **lons, double **lats, double **vals)
{
    char buffer[1024] = {};
    int numRows = 0;

    /*-----------------------------------------------------------------------------
     * Find out number of entries in the file
     *-----------------------------------------------------------------------------*/
    {
        FILE *f = fopen(fileName.c_str(), "r");
        while(fgets(buffer, sizeof(buffer), f)!=NULL) numRows++;
        fclose(f);
    }

    
    /*-----------------------------------------------------------------------------
     * Allocate memory
     *-----------------------------------------------------------------------------*/
    *nNodes = numRows;
    *lons = new double[numRows];
    *lats = new double[numRows];
    *vals = new double[numRows];

    /*-----------------------------------------------------------------------------
     * Read values 
     *-----------------------------------------------------------------------------*/
    {
        double eps=1e-5;
        FILE *f = fopen(fileName.c_str(), "r");
        int counter = 0;
        double sf = atan(1.)/45.;

        while(fgets(buffer, sizeof(buffer), f)!=NULL) 
        {
            double lon,lat,val;

            sscanf(buffer, "%lf %lf %lf", &lon, &lat, &val);
            
            (*lons)[counter] = lon * sf; /* convert to radians */
            (*lats)[counter] = lat * sf; /* convert to radians */
            (*vals)[counter] = val;
            
            if( ((*lons)[counter] < 0)      ||
                ((*lons)[counter] > (2*PI+eps))   ||
                ((*lats)[counter] < (-PI/2-eps)) || 
                ((*lats)[counter] > (PI/2+eps)) )
            {
                cout << (*lons)[counter] << endl;
                cout << (*lats)[counter] << endl;

                cout << "Error in input: lat not within [-PI/2, PI/2] or lon not within [0, 2*PI]" << endl;
                exit(0);
            }

            counter++;
        }
        
        fclose(f);
    }     
}

int main(int argc, char **argv)
{
    if(argc != 3)
    {
        cout << "\nUsage: ./sphericalinterpolate <triplet-data-file> <output-file> \n\nNote: data-file should have the format <lon lat val>. Latitudes should be within [-PI/2, PI/2]. Longitudes should be within [0, 2*PI]" << endl;
        return 0;
    }

    int nNodes = 0;
    double *lons, *lats, *vals;
    

    /*-----------------------------------------------------------------------------
     * Read data
     *-----------------------------------------------------------------------------*/
    ReadGrid(argv[1], &nNodes, &lons, &lats, &vals);

    
    /*-----------------------------------------------------------------------------
     * Create Spherical Interpolation Object
     *-----------------------------------------------------------------------------*/
    SphericalSplineInterpolator SPI(lons, lats, vals, nNodes);
    

    /*-----------------------------------------------------------------------------
     * Create user-defined evaluation points and test interpolation functionality
     *-----------------------------------------------------------------------------*/
    const int NLON = 361;
    const int NLAT = 181;
    
    double *userLons             = new double[NLON*NLAT];
    double *userLats             = new double[NLON*NLAT];
    double *interpolatedValues   = new double[NLON*NLAT];
    double sf                    = atan(1.)/45.;
    
    for(int i=0; i<NLON; i++)
    {
        for(int j=0; j<NLAT; j++)
        {
            userLons[i*NLAT + j] = i * sf;
            userLats[i*NLAT + j] = (j-90) * sf;
        }
    }
    
    /* Interpolate values at user-defined locations */
    SPI.InterpolateC0(NLON*NLAT, userLons, userLats, interpolatedValues);

    /* Output user-defined locations and interpolated values */
    {
        FILE *fo = fopen(argv[2], "w");
        
        for(int i=0; i<NLON; i++)
        {
            for(int j=0; j<NLAT; j++)
            {
                fprintf(fo, "%e %e %e\n", userLons[i*NLAT+j]/sf, userLats[i*NLAT+j]/sf, interpolatedValues[i*NLAT+j]);
            }
        }

        fclose(fo);
    }

    /*-----------------------------------------------------------------------------
     * Cleanup
     *-----------------------------------------------------------------------------*/
    delete [] lons;
    delete [] lats;
    delete [] vals;
    
    delete [] userLons;
    delete [] userLats;
    delete [] interpolatedValues;

    return 1;
}

