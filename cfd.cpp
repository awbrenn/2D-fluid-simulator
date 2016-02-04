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
  gravityX = 0.0;
  gravityY = (float)(-9.8);
  density1 = new float[Nx*Ny]();
  density2 = new float[Nx*Ny]();
  velocity1 = new float[Nx*Ny*2]();
  velocity2 = new float[Nx*Ny*2]();
  color1 = new float[Nx*Ny*3]();
  color2 = new float[Nx*Ny*3]();
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

void cfd::sources()
{
  // add sources
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

  // compute sources
  for (int j=0; j<Ny; ++j)
  {
    for (int i=0; i<Nx; ++i)
    {
      velocity1[vIndex(i,j,0)] = gravityX * density1[dIndex(i,j)];
      velocity1[vIndex(i,j,1)] = gravityY * density1[dIndex(i,j)];
    }
  }
}