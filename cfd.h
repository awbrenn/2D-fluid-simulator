//
// Created by awbrenn on 1/20/16.
//

#ifndef ADVECTION_CFD_H
#define ADVECTION_CFD_H


class cfd
{
  public:
    // constructors/destructors
    cfd(const int nx, const int ny, const float dx, const float dt);
    ~cfd();

    // public methods
    void advect();
    void sources();
    void addSourceColor();
    void addSourceDensity();

    // getters
    int getNx()           const { return Nx; }
    int getNy()           const { return  Ny; }
    float getDx()         const { return Dx; }
    float getGravityX()   const { return gravityX;}
    float getGravityY()   const { return gravityY;}
    float* getDensity()  const { return density1; }
    float* getVelocity1() const { return velocity1; }
    float* getVelocity2() const { return velocity2; }
    float* getColor()    const { return color1; }
    float* getDensitySourceField() { return densitySourceField; }

    // setters
    void setDensitySourceField(float* dsrc) { densitySourceField = dsrc; }
    void setColorSourceField(float* csrc)   { colorSourceField = csrc; }

    // indexing
    int dIndex(int i, int j)        const { return i+Nx*j; }
    int pIndex(int i, int j)        const { return i+Nx*j; }
    int vIndex(int i, int j, int c) const { return (i+Nx*j)*2+c; }
    int cIndex(int i, int j, int c) const { return (i+Nx*j)*3+c; }
    int clampUpperBound(int index, int upper_bound) { return index<upper_bound  ? index : upper_bound; }
    int clampLowerBound(int index, int lower_bound) { return index>=lower_bound ? index : lower_bound; }

  private:
    int     Nx, Ny;
    int     nloops;
    float   Dx;
    float   dt;
    float   gravityX, gravityY;
    float   *density1, *density2;
    float   *velocity1, *velocity2;
    float   *color1, *color2;
    float   *divergence;
    float   *pressure;
    float   *densitySourceField;
    float   *colorSourceField;

    // private methods
    void bilinearlyInterpolate(const int ii, const int jj, const float x, const float y);
    void computeVelocity(float force_x, float force_y);
    void computeDivergence();
    void computePressure();
    void computePressureForces(int i, int j, float* force_x, float* force_y);
    void computeVelocityBasedOnPressureForces();
    float InterpolateColor(int i, int j, int c, float w1, float w2, float w3, float w4);
    float InterpolateVelocity(int i, int j, int c, float w1, float w2, float w3, float w4);
    float InterpolateDensity(int i, int j, float w1, float w2, float w3, float w4);
};

#endif //ADVECTION_CFD_H
