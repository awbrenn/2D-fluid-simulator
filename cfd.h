//
// Created by awbrenn on 1/20/16.
//

#ifndef ADVECTION_CFD_H
#define ADVECTION_CFD_H


class cfd
{

  public:
    // constructors/destructors
    cfd(int nx, int ny, float dx);
    ~cfd();

    // public methods
    void advect(const float dt = (float)(1.0/24.0));
    void sources();

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

  private:
    int     Nx, Ny;
    int     nloops;
    float   Dx;
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
    void addSourceColor();
    void addSourceDensity();
    void computeVelocity(float force_x, float force_y);
    void computeDivergence();
    void computePressure();
    void computePressureForces(int i, int j, float* force_x, float* force_y);
    void computeVelocityBasedOnPressureForces();
};

#endif //ADVECTION_CFD_H
