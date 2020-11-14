#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"


    //  Kuzmin disk structure
    struct Kuzmin_in {
      Real M; // Total Mass
      Real a; // Scale Length
    };
    struct Kuzmin_out {
      Real force_r;     // grad Phi (r component)
      Real force_theta; // grad Phi (theta component)
    };


    //*********************************************************************************
    // Get parameters (total mass and scale length) for Kuzmin disk potential
    // Input: 
    //    Two dimentional array stellarMass[row][col] containing
    //    stellar masses (without ghost cells) between upper and lower   
    //    cells in theta direction (row) and all cells in radial direction (col)
    // Output:
    //    Kuzmin_in structure containing total mass (M) and length scale (a)
    //*********************************************************************************
    template <int rows, int cols>
    Kuzmin_in getKuzminParams(MeshBlock *pmb, Real (&stellarMass)[rows][cols])
    {
      Real sumMass[cols]={0};
      Real halfMass, ihalfRadius;  
      Kuzmin_in s;

      // Sum over upper and lower 2 cells in theta direction
      for (int j=0; j<cols; j++)
        for (int i=0; i<rows; i++){
          sumMass[j]+=stellarMass[i][j];
        }

      // Calculate accumulated mass array
      for (int idx=1; idx<cols; idx++){
        sumMass[idx] += sumMass[idx-1];
      }

      halfMass = sumMass[cols-1] / 2.;

      // Calculate half mass radius index
      for (int idx=0; idx<cols; idx++)
        if (sumMass[idx] <= halfMass)
            ihalfRadius = idx;

      s.M = sumMass[cols-1];
      s.a = pmb->pcoord->x1v(ihalfRadius + pmb->is) / sqrt(3); 
      return s;
    }
   



    //*********************************************************************************
    // Calculate gradient of the Kuzmin potential  
    // Input:
    //        i = index to radial component of the cells where gradient phi is computed
    //        j = index to theta component of the cells where gradient phi is computed
    //        Kuzim_in = structure contains the Kuzim disk parameters
    // Output:
    //        Kuzmin_out structure contains radial and azimuthal components of force per
    //        unit mass (negative of analytical gradient of Kuzmin disk potential) at 
    //        location indicated by radial and azimutal index i & j
    //********************************************************************************* 
    Kuzmin_out getForce(MeshBlock *pmb, int i, int j, Kuzmin_in s_in, Real G)
    {
        Real r, theta, rsin, rcos, denom;
        Kuzmin_out s_out;

        r = pmb->pcoord->x1v(i);
        theta = pmb->pcoord->x2v(j);
        rsin = r * sin(theta);
        rcos = r * cos(theta);
        denom = pow(SQR(rsin) + SQR(s_in.a + abs(rcos)), 1.5);
        s_out.force_r =-1 * G * s_in.M * (r + s_in.a * abs(cos(theta))) / denom;
        s_out.force_theta = ((rcos >=0)?1.:-1.) * (G * s_in.M * s_in.a * rsin)/denom/r;
        return s_out;
    };



    //*********************************************************************************
    // Calculate value of the Kuzmin potential  
    // Input:
    //        i = index to radial component of the cells where potential is computed
    //        j = index to theta component of the cells where potential is computed
    //        Kuzim_in = structure contains the Kuzim disk parameters
    //*********************************************************************************
    Real getKumzminPotential(MeshBlock *pmb, int i, int j, Kuzmin_in s_in, Real G)
    {
        Real r, theta, rsin, rcos, denom;

        r = pmb->pcoord->x1v(i);
        theta = pmb->pcoord->x2v(j);
        rsin = r * sin(theta);
        rcos = r * cos(theta);
        denom = sqrt(SQR(rsin) + SQR(abs(rcos)+s_in.a));
        return -1 * G * s_in.M / denom;
    }