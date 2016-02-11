//
// Created by awbrenn on 1/20/16.
//
#include <cmath>
#include "cfd.h"
#include "cfdUtility.h"


cfd::cfd(const int nx, const int ny, const float dx, const float Dt, int Nloops)
{
  Nx = nx;
  Ny = ny;
  Dx = dx;
  dt = Dt;
  nloops = Nloops;
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


float cfd::getDensity(int i, int j)
{
  if (i < Nx && i >=0 && j < Ny && j >= 0)
    return density1[dIndex(i,j)];
  else
    return 0.0;
}


float cfd::getVelocity(int i, int j, int c)
{
  if (i < Nx && i >=0 && j < Ny && j >= 0)
    return velocity1[vIndex(i,j,c)];
  else
    return 0.0;
}


float cfd::getColor(int i, int j, int c)
{
  if (i < Nx && i >=0 && j < Ny && j >= 0)
    return color1[cIndex(i,j,c)];
  else
    return 0.0;
}


float cfd::getPressure(int i, int j)
{
  if (i < Nx && i >=0 && j < Ny && j >= 0)
    return pressure[pIndex(i,j)];
  else
    return 0.0;
}


float cfd::getDivergence(int i, int j)
{
  if (i < Nx && i >=0 && j < Ny && j >= 0)
    return divergence[dIndex(i,j)];
  else
    return 0.0;
}


float cfd::InterpolateDensity(int i, int j, float w1, float w2, float w3, float w4)
{
  return getDensity(i    , j)     * w1 +
         getDensity(i + 1, j)     * w2 +
         getDensity(i    , j + 1) * w3 +
         getDensity(i + 1, j + 1) * w4;
}

float cfd::InterpolateVelocity(int i, int j, int c, float w1, float w2, float w3, float w4)
{
  return getVelocity(i    , j,     c) * w1 +
         getVelocity(i + 1, j,     c) * w2 +
         getVelocity(i,     j + 1, c) * w3 +
         getVelocity(i + 1, j + 1, c) * w4;
}

float cfd::InterpolateColor(int i, int j, int c, float w1, float w2, float w3, float w4)
{
  return getColor(i    , j,     c) * w1 +
         getColor(i + 1, j,     c) * w2 +
         getColor(i,     j + 1, c) * w3 +
         getColor(i + 1, j + 1, c) * w4;
}


void cfd::bilinearlyInterpolate(const int ii, const int jj, const float x, const float y)
{
  // get index of sample
  const int i = (int) (x/Dx);
  const int j = (int) (y/Dx);

  // get weights of samples
  const float ax = std::abs(x/Dx - i);
  const float ay = std::abs(y/Dx - j);
  const float w1 = (1-ax) * (1-ay);
  const float w2 = ax * (1-ay);
  const float w3 = (1-ax) * ay;
  const float w4 = ax * ay;

  if (i < Nx && i >=0 && j < Ny && j >= 0) {

    density2[dIndex(ii, jj)] = InterpolateDensity(i, j, w1, w2, w3, w4);

    velocity2[vIndex(ii, jj, 0)] = InterpolateVelocity(i, j, 0, w1, w2, w3, w4);
    velocity2[vIndex(ii, jj, 1)] = InterpolateVelocity(i, j, 1, w1, w2, w3, w4);

    color2[cIndex(ii, jj, 0)] = InterpolateColor(i, j, 0, w1, w2, w3, w4);
    color2[cIndex(ii, jj, 1)] = InterpolateColor(i, j, 1, w1, w2, w3, w4);
    color2[cIndex(ii, jj, 2)] = InterpolateColor(i, j, 2, w1, w2, w3, w4);
  }
  else // you are out of bounds
  {
    density2[dIndex(ii, jj)]   = 0.0f;
    velocity2[vIndex(ii,jj,0)] = 0.0f;
    velocity2[vIndex(ii,jj,1)] = 0.0f;
    color2[cIndex(ii,jj,0)]    = 0.0f;
    color2[cIndex(ii,jj,1)]    = 0.0f;
    color2[cIndex(ii,jj,2)]    = 0.0f;
  }

}


void cfd::advect()
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
      velocity1[vIndex(i,j,0)] += force_x * density1[dIndex(i,j)]*dt;
      velocity1[vIndex(i,j,1)] += force_y * density1[dIndex(i,j)]*dt;
    }
  }
}


void cfd::computeDivergence()
{
  for (int j = 0; j < Ny; ++j)
  {
    for (int i = 0; i < Nx; ++i)
    {
      divergence[dIndex(i,j)] = (getVelocity(i+1, j,   0) -
                                 getVelocity(i-1, j,   0)) / (2*Dx) +
                                (getVelocity(i,   j+1, 1) -
                                 getVelocity(i,   j-1, 1)) / (2*Dx);
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
        pressure[pIndex(i,j)] = ((getPressure(i+1, j)     +
                                   getPressure(i-1, j)    +
                                   getPressure(i,   j+1)  +
                                   getPressure(i,   j-1)) *
                                   0.25f) - ((Dx*Dx/4.0f) * getDivergence(i,j));
      }
    }
  }
}


void cfd::computePressureForces(int i, int j, float* force_x, float* force_y)
{
  *force_x = (getPressure(i+1, j) - getPressure(i-1, j)) / (2*Dx);
  *force_y = (getPressure(i, j+1) - getPressure(i, j-1)) / (2*Dx);
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