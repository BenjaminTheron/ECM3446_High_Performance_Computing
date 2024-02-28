/*******************************************************************************
2D advection example program which advects a Gaussian u(x,y) at a fixed velocity



Outputs: initial.dat - inital values of u(x,y) 
         final.dat   - final values of u(x,y)

         The output files have three columns: x, y, u

         Compile with: gcc -o advection2D -std=c99 advection2D.c -lm

Notes: The time step is calculated using the CFL condition

********************************************************************************/

/*********************************************************************
                     Include header files 
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/*********************************************************************
                      Function propotypes
**********************************************************************/
float horizontal_velocity(float z, float z0, float u_star, float k);

/*********************************************************************
                      Main function
**********************************************************************/
int main(){

  /* Grid properties */
  const int NX=1000;    // Number of x points
  const int NY=1000;    // Number of y points
  const float xmin=0.0; // Minimum x value
  const float xmax=30.0; // Maximum x value
  const float ymin=0.0; // Minimum y value
  const float ymax=30.0; // Maximum y value
  
  /* Parameters for the Gaussian initial conditions */
  const float x0=3.0;                    // Centre(x)
  const float y0=15.0;                    // Centre(y)
  const float sigmax=1.00;               // Width(x)
  const float sigmay=5.00;               // Width(y)
  const float sigmax2 = sigmax * sigmax; // Width(x) squared
  const float sigmay2 = sigmay * sigmay; // Width(y) squared

  /* Boundary conditions */
  const float bval_left=0.0;    // Left boudnary value
  const float bval_right=0.0;   // Right boundary value
  const float bval_lower=0.0;   // Lower boundary
  const float bval_upper=0.0;   // Upper bounary
  
  /* Time stepping parameters */
  const float CFL=0.9;   // CFL number 
  /* DOES NOT COME FROM CHANGING THE CFL NUMBER */

  const int nsteps=800; // Number of time steps

  /* Velocity */
  // As per the requirements of Q2.3 horizontal velocity is no longer
  // a constant value but instead varies according to the height
  //const float velx=1.00; // Velocity in x direction
  const float vely=0.00; // Velocity in y direction
  
  /* Arrays to store variables. These have NX+2 elements
     to allow boundary values to be stored at both ends */
  float x[NX+2];          // x-axis values
  float y[NX+2];          // y-axis values
  float u[NX+2][NY+2];    // Array of u values
  float dudt[NX+2][NY+2]; // Rate of change of u

  float x2;   // x squared (used to calculate iniital conditions)
  float y2;   // y squared (used to calculate iniital conditions)
  
  /* Calculate distance between points */
  float dx = (xmax-xmin) / ( (float) NX);
  float dy = (ymax-ymin) / ( (float) NY);
 
  /* Parameters for the logarithmic profile required by question 2.3 */
  const float u_star = 0.2;  // friction velocity                                         
  const float z0 = 1.0;      // roughness length                                           
  const float k = 0.41;      // Von Karman's constant
 
  /* Calculate time step using the CFL condition */
  /* The fabs function gives the absolute value in case the velocity is -ve */
  float dt = CFL / ( (fabs(horizontal_velocity(-1, z0, u_star, k)) / dx) + (fabs(vely) / dy) );

  /*** Report information about the calculation ***/
  printf("Grid spacing dx     = %g\n", dx);
  printf("Grid spacing dy     = %g\n", dy);
  printf("CFL number          = %g\n", CFL);
  printf("Time step           = %g\n", dt);
  printf("No. of time steps   = %d\n", nsteps);
  printf("End time            = %g\n", dt*(float) nsteps);
  printf("Distance advected x = %g\n", horizontal_velocity(-1, z0, u_star, k)*dt*(float) nsteps);
  printf("Distance advected y = %g\n", vely*dt*(float) nsteps);

  /*** Place x points in the middle of the cell ***/
  /* LOOP 1 */
  /* CAN BE PARALLELISED */
  #pragma omp parallel for default(none) shared(x, dx)
  for (int i=0; i<NX+2; i++) {
    x[i] = ( (float) i - 0.5) * dx;
  }

  /*** Place y points in the middle of the cell ***/
  /* LOOP 2 */
  /* CAN BE PARALLELISED */
  #pragma omp parallel for default(none) shared(y, dy)
  for (int j=0; j<NY+2; j++) {
    y[j] = ( (float) j - 0.5) * dy;
  }

  /*** Set up Gaussian initial conditions ***/
  /* LOOP 3 */
  /* CAN BE PARALLELISED */
  #pragma omp parallel for collapse(2) default(none) shared(x, y, u) private(x2, y2)
  for (int i=0; i<NX+2; i++) {
    for (int j=0; j<NY+2; j++) {
      x2      = (x[i]-x0) * (x[i]-x0);
      y2      = (y[j]-y0) * (y[j]-y0);
      u[i][j] = exp( -1.0 * ( (x2/(2.0*sigmax2)) + (y2/(2.0*sigmay2)) ) );
    }
  }

  /*** Write array of initial u values out to file ***/
  FILE *initialfile;
  initialfile = fopen("initial.dat", "w");
  /* LOOP 4 */
  /* CANNOT BE PARALLELISED - This creates an output dependency (write after write)
  /* as executing this loop in parallel leads to the results being written into */
  /* the file in the incorrect order. (Here the results must be printed in the exact */
  /* same order they stored in memory - next to one another) */
  for (int i=0; i<NX+2; i++) {
    for (int j=0; j<NY+2; j++) {
      fprintf(initialfile, "%g %g %g\n", x[i], y[j], u[i][j]);
    }
  }
  fclose(initialfile);

  /*** Update solution by looping over time steps ***/
  /* LOOP 5 */
  /* CANNOT BE PARALLELISED - The timesteps need to be calculated sequentially */
  /* in order to produce the correct results. There is a flow data dependency as */
  /* time step i requires time step i - 1 (the u value) to have been calculated.*/
  for (int m=0; m<nsteps; m++) {
   
    /*** Apply boundary conditions at u[0][:] and u[NX+1][:] ***/
    /* LOOP 6 */
    /* CAN BE PARALLELISED */
    #pragma omp parallel for default(none) shared(u)
    for (int j=0; j<NY+2; j++) {
      u[0][j]    = bval_left;
      u[NX+1][j] = bval_right;
    }

    /*** Apply boundary conditions at u[:][0] and u[:][NY+1] ***/
    /* LOOP 7 */
    /* CAN BE PARALLELISED */
    #pragma omp parallel for default(none) shared(u)
    for (int i=0; i<NX+2; i++) {
      u[i][0]    = bval_lower;
      u[i][NY+1] = bval_upper;
    }
    
    /*** Calculate rate of change of u using leftward difference ***/
    /* Loop over points in the domain but not boundary values */
    /* LOOP 8 */
    /* CAN BE PARALLELISED */
    /* The u array is fully populated in a prior loop, so there is no risk of a */
    /* flow data dependency as u[i] does not need u[i-1] before it is updated */
    #pragma omp parallel for default(none) shared(dx, dy, u, dudt, y)
    for (int i=1; i<NX+1; i++) {
      for (int j=1; j<NY+1; j++) {
        dudt[i][j] = -horizontal_velocity(y[j], z0, u_star, k) * (u[i][j] - u[i-1][j]) / dx
		     - vely * (u[i][j] - u[i][j-1]) / dy;
      }
    }
    /* Decrease the size of the time step so that that graph is output correctly by*/ 
    /* using the correct horizontal velocity value. This occurs when applying the maximum vertical*/ 
    /* velocity to the logarithmic profile used for Q2.3, as horizontal velocity now varies with height. */
    dt = CFL / ( (fabs(horizontal_velocity(y[NX+1], z0, u_star, k)) / dx) + (fabs(vely) / dy) );
    
    /*** Update u from t to t+dt ***/
    /* Loop over points in the domain but not boundary values */
    /* LOOP 9 */
    /* CAN BE PARALLELISED  */
    /* The dudt array is fully populated and as the index being looked at is */
    /* the same for the u array as it is for the dudt arrary, there is no output*/
    /* or flow data dependency */
    #pragma omp parallel for collapse(2) default(none) shared(u, dudt, dt)
    for(int i=1; i<NX+1; i++) {
      for (int j=1; j<NY+1; j++) {
        u[i][j] = u[i][j] + dudt[i][j] * dt;
      }
    }
    
  } // time loop
  
  /*** Write array of final u values out to file ***/
  FILE *finalfile;
  finalfile = fopen("final.dat", "w");
  /* LOOP 10 */
  /* CANNOT BE PARALLELISED - This creates an output dependency (write after write) */
  /* as executing this loop in parallel leads to the results being written into */
  /* the file in the incorrect order. (Here the results must be printed in the exact */
  /* same order they stored in memory - next to one another) */
  for (int i=0; i<NX+2; i++) {
    for (int j=0; j<NY+2; j++) {
      fprintf(finalfile, "%g %g %g\n", x[i], y[j], u[i][j]);
    }
  }
  fclose(finalfile);

  /* CALCULATE THE VERTICALLY AVERAGED DISTRIBUTION OF u(x,y)
  /* Calculate the vertical average of the values across the whole vertical domain */
  float vertical_average[NX+2];
  float total = 0.0;
  for (int i = 0; i < NX+2; i++) {
    for (int j = 1; j < NY+1; j++) {
      total += u[i][j];
    }
    vertical_average[i] = total/NY;
    total = 0.0;
  }

  /* Write the vertically averaged distribution to a file */
  FILE *vertical_distribution_file;
  vertical_distribution_file = fopen("vert_avg_dist_final.dat", "w");

  /* If no file exists with the provided name, the file pointer will point to null */
  if (vertical_distribution_file == NULL) {
    printf("ERROR: FILE COULD NOT BE OPENED %s\n", "vert_avg_dist_final.dat");
    exit(8);
  }

  /* Write the values stored in the vertical_average array to the file */
  for (int i = 0; i < NX+2; i++) {
    for (int j = 0; j < NY+2; j++) {
      fprintf(vertical_distribution_file, "%f %f\n", x[i], vertical_average[i]);
    }
  }
  fclose(vertical_distribution_file);

  return 0;
}


float horizontal_velocity(float z, float z0, float u_star, float k) {
  if (z > z0) {
    return (u_star/k) * log(z/z0);
  }
  /* Used to calculate the intial value of dt (when velx = 1) */
  else if (z == -1) {
    return 1.0;
  }
  else {
    return 0;
  }
}


/* End of file ******************************************************/
