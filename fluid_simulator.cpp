//------------------------------------------------
//
//  cfd_paint
//
//
//-------------------------------------------------

//-------------------------------------------------
//
//  usage:
//
//  cfd_paint is an interactive paint program
//  in which the user paints density, color, and
//  or divergence sources that flow using
//  computational fluid dynamics and react with 
//  obstructions in the space.
//
//  There are two paint modes.  Typing 'o' puts the
//  program in obstruction painting mode. When you
//  hold down the left mouse button and paint, you
//  will see a black obstruction painted.  This 
//  obstruction may be any shape.
//
//  Typing 's' puts the program in source painting 
//  mode.  Now painting with the left mouse button
//  down injects density into the simulation.
//  The flow it produces evolves as you
//  continue to paint.  The flow bounces off any
//  obstructions that have been painted or are
//  subsequently painted.
//
//  Typing 'b' clears all obstructions, flow, density,
//  and color.
//
//  Typing '=' and '-' brightens and darkens the display.
//
//  Pressing the spacebar starts and stops the flow 
//  evolution. While the evolution is stopped, you
//  can continue painting obstructions.
//
//
//
//-------------------------------------------------




#include <cmath>
#include "CmdLineFind.h"

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.
#include "cfd.h"

#include <iostream>
#include <OpenImageIO/imageio.h>
#include <omp.h>


using namespace std;
using namespace lux;
OIIO_NAMESPACE_USING

int iwidth, iheight;
float* display_map;
float* density_source;
float* color_source;
float* obstruction_source;
cfd *fluid;

int paint_mode;
enum{ PAINT_OBSTRUCTION, PAINT_SOURCE, PAINT_DIVERGENCE, PAINT_COLOR };

bool toggle_animation_on_off;

float scaling_factor;

#define BRUSH_SIZE 11
float obstruction_brush[BRUSH_SIZE][BRUSH_SIZE];
float source_brush[BRUSH_SIZE][BRUSH_SIZE];

int xmouse_prev, ymouse_prev;

////////  OpenImageIO reader

int readOIIOImage( const char* fname)
{
  int channels;
  ImageInput *in = ImageInput::create (fname);
  if (! in) { return -1; }
  ImageSpec spec;
  in->open (fname, spec);
  iwidth = spec.width; // note iwidth and iheight are set to to image size
  iheight = spec.height;
  channels = spec.nchannels;
  float* pixels = new float[iwidth*iheight*channels];
  color_source = new float[iwidth*iheight*channels]; // allocate appropriate space for image

  in->read_image (TypeDesc::FLOAT, pixels);
  long index = 0;
  for( int j=0;j<iheight;j++)
  {
    for( int i=0;i<iwidth;i++ )
    {
      for( int c=0;c<channels;c++ )
      {
        color_source[ (i + iwidth*(iheight - j - 1))*channels + c ] = pixels[index++];
      }
    }
  }
  delete pixels;

  in->close ();
  delete in;

  return 0;
}


//--------------------------------------------------------
//
//  Initialization routines
//
//  
// Initialize all of the fields to zero

void InitializeBrushes()
{
  int brush_width = (BRUSH_SIZE-1)/2;
  for( int j=-brush_width;j<=brush_width;j++ )
  {
    int jj = j + brush_width;
    float jfactor =  (float(brush_width) - (float)fabs(j) )/float(brush_width);
    for( int i=-brush_width;i<=brush_width;i++ )
    {
      int ii = i + brush_width;
      float ifactor =  (float(brush_width) - (float)fabs(i) )/float(brush_width);
      float radius = (float) ((jfactor * jfactor + ifactor * ifactor) / 2.0);
      source_brush[ii][jj] = pow(radius,0.5);
      obstruction_brush[ii][jj] = (float)(1.0 - pow(radius, 1.0/4.0));
    }
  }
}


void setNbCores( int nb )
{
  omp_set_num_threads( nb );
}


//----------------------------------------------------

void ConvertToDisplay()
{
  float *color = fluid->getColorPointer();
  for( int j=0;j<iheight;j++ )
  {
#pragma omp parallel for
    for(int i=0;i<iwidth;i++ )
    {
      int index = i + iwidth*j;
      float r,g,b;
      r = color[index*3];
      g = color[index*3+1];
      b = color[index*3+2];
      display_map[3*index  ] = r * scaling_factor;
      display_map[3*index+1] = g * scaling_factor;
      display_map[3*index+2] = b * scaling_factor;
    }
  }
}


//------------------------------------------
//
//  Painting and display code
//

void resetScaleFactor( float amount )
{
   scaling_factor *= amount;
}



void DabSomePaint( int x, int y )
{
  int brush_width = (BRUSH_SIZE-1)/2;
  int xstart = x - brush_width;
  int ystart = y - brush_width;
  if( xstart < 0 ){ xstart = 0; }
  if( ystart < 0 ){ ystart = 0; }

  int xend = x + brush_width;
  int yend = y + brush_width;
  if( xend >= iwidth ){ xend = iwidth-1; }
  if( yend >= iheight ){ yend = iheight-1; }

  if( paint_mode == PAINT_OBSTRUCTION )
  {
    for(int ix=xstart;ix <= xend; ix++)
    {
      for( int iy=ystart;iy<=yend; iy++)
      {
        int index = ix + iwidth*(iheight-iy-1);
//        color_source[3 * index] *= obstruction_brush[ix - xstart][iy - ystart];
//        color_source[3 * index + 1] *= obstruction_brush[ix - xstart][iy - ystart];
//        color_source[3 * index + 2] *= obstruction_brush[ix - xstart][iy - ystart];
        obstruction_source[index] *= obstruction_brush[ix - xstart][iy - ystart];
//        if (obstruction_source[index] != 1.0)
//          cout << index << endl;
      }
    }
//  fluid->setColorSourceField(color_source);
    fluid->setObstructionSourceField(obstruction_source);
  }
  else if( paint_mode == PAINT_SOURCE )
  {
    for(int ix=xstart;ix <= xend; ix++)
    {
      for( int iy=ystart;iy<=yend; iy++)
      {
        int index = ix + iwidth*(iheight-iy-1);
        color_source[3 * index] += source_brush[ix - xstart][iy - ystart];
        color_source[3 * index + 1] += source_brush[ix - xstart][iy - ystart];
        color_source[3 * index + 2] += source_brush[ix - xstart][iy - ystart];
        density_source[index] += source_brush[ix-xstart][iy-ystart];
      }
    }
    fluid->setColorSourceField(color_source);
    fluid->setDensitySourceField(density_source);
  }

  return;
}


//----------------------------------------------------
//
//  GL and GLUT callbacks
//
//----------------------------------------------------


void cbDisplay( void )
{
  glClear(GL_COLOR_BUFFER_BIT );
  glDrawPixels( iwidth, iheight, GL_RGB, GL_FLOAT, display_map );
  glutSwapBuffers();
}


void update()
{
  fluid->advect();
  fluid->sources();
}

// animate and display new result
void cbIdle()
{
  update();
  ConvertToDisplay();
  glutPostRedisplay(); 
}


void cbOnKeyboard( unsigned char key, int x, int y )
{
  switch (key) 
  {
    case '-': case '_':
    resetScaleFactor( 0.9 );
    break;

    case '+': case '=':
    resetScaleFactor( (float)(1.0/0.9) );
    break;

    case 'r':
    scaling_factor = 1.0;
    break;

    case ' ':
    toggle_animation_on_off = !toggle_animation_on_off;

    case 'o':
    paint_mode = PAINT_OBSTRUCTION;
    break;

    case 's':
    paint_mode = PAINT_SOURCE;
    break;

    default:
    break;
  }
}


void cbMouseDown( int button, int state, int x, int y )
{
  if( button != GLUT_LEFT_BUTTON ) { return; }
  if( state != GLUT_DOWN ) { return; }
  xmouse_prev = x;
  ymouse_prev = y;
  DabSomePaint( x, y );
}


void cbMouseMove( int x, int y )
{
  xmouse_prev = x;
  ymouse_prev = y;
  DabSomePaint( x, y ); 
}


void PrintUsage()
{
  cout << "cfd_paint keyboard choices\n";
  cout << "s       turns on painting source strength\n";
  cout << "o       turns on painting obstructions\n";
  cout << "+/-     increase/decrease brightness of display\n";
  cout << "r       resets brightness to default\n";
}


//---------------------------------------------------

int main(int argc, char** argv)
{
  CmdLineFind clf(argc, argv);

  iwidth = clf.find("-NX", 512, "Horizontal grid points");
  iheight = clf.find("-NY", iwidth, "Vertical grid points");

  int nloops = clf.find("-nloops", 3, "Number of loops over pressure.");

  setNbCores(4);

  string imagename = clf.find("-image", "none", "Image to drive color");

  clf.usage("-h");
  clf.printFinds();
  PrintUsage();

  // initialize a few variables
  scaling_factor = 1.0;
  toggle_animation_on_off = true;

  // if reading the image fails we need to allocate space for color_source
  if (readOIIOImage(imagename.c_str()) != 0)
    color_source = new float[iwidth*iheight*3]();

  density_source = new float[iwidth*iheight]();
  obstruction_source = new float[iwidth*iheight];
  //Initialize(obstruction_source, iwidth*iheight, 1.0);
  for(int i=0;i<iwidth*iheight;i++ ) { obstruction_source[i] = 1.0; }

  // initialize fluid
  fluid = new cfd(iwidth, iheight, 1.0, (float)(1.0/24.0), nloops);
  fluid->setColorSourceField(color_source);

  display_map = new float[iwidth*iheight*3];

  InitializeBrushes();

  paint_mode = PAINT_SOURCE;

  // GLUT routines
  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(iwidth, iheight);

  // Open a window 
  char title[] = "Fluid Simulator";
  glutCreateWindow(title);

  glClearColor(1, 1, 1, 1);

  glutDisplayFunc(&cbDisplay);
  glutIdleFunc(&cbIdle);
  glutKeyboardFunc(&cbOnKeyboard);
  glutMouseFunc(&cbMouseDown);
  glutMotionFunc(&cbMouseMove);

  glutMainLoop();
  return 1;
}
