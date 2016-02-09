//
// Created by awbrenn on 1/20/16.
//
#include <cmath>
#include "cfd.h"
#include "cfdUtility.h"


// dealing with negative results
int mod(int a, int b)
{
  return (a%b+b)%b;
}


cfd::cfd(const int nx, const int ny, const float dx)
{
  Nx = nx;
  Ny = ny;
  Dx = dx;
  nloops = 10;
  gravityX = 0.0;
  gravityY = (float)(-9.8);
  density1 = new float[Nx*Ny]();
  density2 = new float[Nx*Ny]();
  velocity1 = new float[Nx*Ny*2]();
  velocity2 = new float[Nx*Ny*2]();
  color1 = new float[Nx*Ny*3]();
  color2 = new float[Nx*Ny*3]();
  divergence = new float[Nx*Ny]();
  pressure = new float[Nx*Ny]();
  densitySourceField = 0;
  colorSourceField = 0;
}


cfd::~cfd()
{
  delete density1;
  delete density2;
  delete velocity1;
  delete velocity2;
  delete color1;
  delete color2;
}


void cfd::bilinearlyInterpolate(const int ii, const int jj, const float x, const float y)
{
  // get index of sample
  const int i = mod((int) (x/Dx), Nx);
  const int j = mod((int) (y/Dx), Ny);

  // get weights of samples
  const float ax = std::abs(x/Dx - ((int)(x/Dx)));
  const float ay = std::abs(y/Dx - ((int)(y/Dx)));
  const float w1 = (1-ax) * (1-ay);
  const float w2 = ax * (1-ay);
  const float w3 = (1-ax) * ay;
  const float w4 = ax * ay;

  density2[dIndex(ii,jj)]    = density1[dIndex(i,j)]                           * w1 +
                               density1[dIndex((i + 1) % Nx, j)]               * w2 +
                               density1[dIndex(i, (j + 1) % Ny)]               * w3 +
                               density1[dIndex((i + 1) % Nx, (j + 1) % Ny)]    * w4;
  velocity2[vIndex(ii,jj,0)] = velocity1[vIndex(i,j,0)]                        * w1 +
                               velocity1[vIndex((i + 1) % Nx, j,0)]            * w2 +
                               velocity1[vIndex(i, (j + 1) % Ny,0)]            * w3 +
                               velocity1[vIndex((i + 1) % Nx, (j + 1) % Ny,0)] * w4;
  velocity2[vIndex(ii,jj,1)] = velocity1[vIndex(i,j,1)]                        * w1 +
                               velocity1[vIndex((i + 1) % Nx, j,1)]            * w2 +
                               velocity1[vIndex(i, (j + 1) % Ny,1)]            * w3 +
                               velocity1[vIndex((i + 1) % Nx, (j + 1) % Ny,1)] * w4;
  color2[cIndex(ii,jj,0)]    = color1[cIndex(i,j,0)]                           * w1 +
                               color1[cIndex((i + 1) % Nx, j,0)]               * w2 +
                               color1[cIndex(i, (j + 1) % Ny,0)]               * w3 +
                               color1[cIndex((i + 1) % Nx, (j + 1) % Ny,0)]    * w4;
  color2[cIndex(ii,jj,1)]    = color1[cIndex(i,j,1)]                           * w1 +
                               color1[cIndex((i + 1) % Nx, j,1)]               * w2 +
                               color1[cIndex(i, (j + 1) % Ny,1)]               * w3 +
                               color1[cIndex((i + 1) % Nx, (j + 1) % Ny,1)]    * w4;
  color2[cIndex(ii,jj,2)]    = color1[cIndex(i,j,2)]                           * w1 +
                               color1[cIndex((i + 1) % Nx, j,2)]               * w2 +
                               color1[cIndex(i, (j + 1) % Ny,2)]               * w3 +
                               color1[cIndex((i + 1) % Nx, (j + 1) % Ny,2)]    * w4;
}


void cfd::advect(const float dt)
{
  float x, y;

  // advect each grid point
  for (int j=0; j<Ny; ++j)
  {
    for (int i=0; i<Nx; ++i)
    {
      x = i*Dx - velocity1[vIndex(i,j,0)]*dt;
      y = j*Dx - velocity1[vIndex(i,j,1)]*dt;
      bilinearlyInterpolate(i, j, x,y);
    }
  }

  swapFloatPointers(&density1, &density2);
  swapFloatPointers(&velocity1, &velocity2);
  swapFloatPointers(&color1, &color2);
}


void cfd::addSourceColor()
{
  if (colorSourceField != 0)
  {
    for (int j=0; j<Ny; ++j)
    {
      for (int i=0; i<Nx; ++i)
      {
        color1[cIndex(i,j,0)] += colorSourceField[cIndex(i,j,0)];
        color1[cIndex(i,j,1)] += colorSourceField[cIndex(i,j,1)];
        color1[cIndex(i,j,2)] += colorSourceField[cIndex(i,j,2)];
      }
    }
    // re-initialize colorSourceField
    Initialize(colorSourceField, Nx*Ny*3, 0.0);
    colorSourceField = 0;
  }
}


void cfd::addSourceDensity()
{
  int index;

  if (densitySourceField != 0)
  {
    for (int j=0; j<Ny; ++j)
    {
      for (int i=0; i<Nx; ++i)
      {
        index = dIndex(i,j);
        density1[index] += densitySourceField[index];
        density1[index] += densitySourceField[index];
        density1[index] += densitySourceField[index];
      }
    }
    // re-initialize colorSourceField
    Initialize(densitySourceField, Nx*Ny, 0.0);
    densitySourceField = 0;
  }
}


void cfd::computeVelocity(float force_x, float force_y)
{
  for (int j=0; j<Ny; ++j)
  {
    for (int i=0; i<Nx; ++i)
    {
      velocity1[vIndex(i,j,0)] = force_x * density1[dIndex(i,j)];
      velocity1[vIndex(i,j,1)] = force_y * density1[dIndex(i,j)];
    }
  }
}


void cfd::computeDivergence()
{
  for (int j = 0; j < Ny; ++j)
  {
    for (int i = 0; i < Nx; ++i)
    {
      divergence[dIndex(i,j)] = ((velocity1[vIndex(i+1,j,0)] - velocity1[vIndex(i-1,j,0)]) / (2*Dx)) +
                                ((velocity1[vIndex(i,j+1,1)] - velocity1[vIndex(i,j-1,0)]) / (2*Dx));
    }
  }
}


void cfd::computePressure()
{
  Initialize(pressure, Nx*Ny, 0.0);

  for(int k = 0; k < nloops; ++k)
  {
    for (int j = 0; j < Ny; ++j)
    {
      for (int i = 0; i < Nx; ++i)
      {
        pressure[pIndex(i,j)] = (float)((pressure[pIndex(i+1,j)] + pressure[pIndex(i-1,j)]  +
                                         pressure[pIndex(i,j+1)] + pressure[pIndex(i,j-1)]) *
                                         (0.25) - (Dx*Dx/4.0)*divergence[dIndex(i,j)]);
      }
    }
  }
}


void cfd::computePressureForces(int i, int j, float* force_x, float* force_y)
{
  *force_x = (pressure[pIndex(i+1,j)] - pressure[pIndex(i-1,j)])/(2*Dx);
  *force_y = (pressure[pIndex(i,j+1)] - pressure[pIndex(i,j-1)])/(2*Dx);
}


void cfd::computeVelocityBasedOnPressureForces()
{
  float force_x, force_y;

  for (int j = 0; j < Ny; ++j)
  {
    for (int i = 0; i < Nx; ++i)
    {
      computePressureForces(i, j, &force_x, &force_y);
      velocity1[vIndex(i,j,0)] -= force_x;
      velocity1[vIndex(i,j,1)] -= force_y;
    }
  }
}


void cfd::sources()
{
  // add sources
  addSourceColor();
  addSourceDensity();

  // compute sources
  computeVelocity(gravityX, gravityY);
  computeDivergence();
  computePressure();
  computeVelocityBasedOnPressureForces();
}