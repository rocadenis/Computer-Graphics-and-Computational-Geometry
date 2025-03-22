/*
  This program plots different 2D functions.
*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <string>
#include <complex>
using namespace std::complex_literals;

//#include "glut.h" //MSVC local library install
#include <GL/glut.h> //system-wide install (or compiler default path)

double circle = atan(1) * 8; 
double halfCircle = atan(1) * 4;
double tau = circle; // 2 * PI = TAU
double pi = halfCircle; // TAU / 2 = PI

int g_w = 800, g_h = 800;

unsigned char g_prevKey;

int g_recursionMax = 8, g_recursionCurrent = 2;
double g_jfa = -0.82, g_jfb = -0.17; //Julia-Fatou a and b values.

//----------------Utility functions----------------------

void bitmapString(void* font, const char* str) {
  //Draw a string, character-by-character.
  char cp;
  for(const char* c = str; *c != 0; ++c) {
    cp = *c; //to respect const
    glutBitmapCharacter(font, cp);
  }
}

void drawBitmapString(const char* str, float x = -2, float y = -2) {
  //Draw a string, optionally setting raster position.
  /*
    We define the convetion that both values -2 mean 'do not change
    raster position'.
  */
  if((-2 != x) || (-2 != y)) {
    glRasterPos2f(x, y);
  }
  //freeglut, not old glut: glutBitmapString(GLUT_BITMAP_8_BY_13, str);
  bitmapString(GLUT_BITMAP_8_BY_13, str);
}

template <typename Numtype>
void drawBitmapNumber(Numtype number, float x = -2, float y = -2) {
  //Convert a number to a string, then draw it.
  //We need the template so we don't display '2' as '2.000000'.
  if((-2 != x) || (-2 != y)) {
    glRasterPos2f(x, y);
  }
  bitmapString(GLUT_BITMAP_8_BY_13, std::to_string(number).c_str());
}

void drawRecursionLevel() {
  //Simple utility function.
  drawBitmapString("Recursion Level: ", -0.98, -0.98);
  /*if we don't set explicit raster positions,
    drawing characters increments the paster position appropriately
  */
  drawBitmapNumber(g_recursionCurrent); 
}

void drawJfConstants() {
  drawBitmapString("Julia-Fatou constants: a = ", -0.98, -0.98);
  drawBitmapNumber(g_jfa);
  drawBitmapString(", b = ");
  drawBitmapNumber(g_jfb);
}
//^^^^^^^^^^^^^^^^^Utility functions^^^^^^^^^^^^^^^^^^


class Turtle {
/*
  Turtle Graphics:
  draw using points, directions and distances.
  (Radial coordinates.)
*/
protected:
  double m_x, m_y;
  double m_angle;

public:
  Turtle(double x = 0, double y = 0):
    m_x(x),
    m_y(y),
    m_angle(0) { }

  void rotate(double angle) {
    m_angle += angle;
  }

  void move(double distance) {
    //Move the Turtle without drawing.
    /*
      We convert from Radial coordinates
      to Cartesian coordinates.
     */
    m_x += distance * cos(m_angle);
    m_y += distance * sin(m_angle);
  }

  void draw(double distance) {
    //Move the Turtle and draw its path.
    glBegin(GL_LINES); {
      glVertex2d(m_x, m_y);
      move(distance);
      glVertex2d(m_x, m_y);
    }
    glEnd();
  }

  void resetPos() {m_x = 0; m_y = 0;}
  void resetRotation() {m_angle = 0;}
};

void drawCircle(double cx, double cy, double radius, int segments) {
  //How to draw a circle with Turtle graphics.
  Turtle t;
  //Arrive at the Cartesian coordinates of the centre.
  t.move(cx);
  t.rotate(pi/2);
  t.move(cy);
  //Reset rotation.
  t.rotate(-pi/2);

  //Arrive on the circle, at angle 0.
  t.move(radius);
  /*
    Up, the tangent on the circle
    (in the trigonometric direction).
  */
  t.rotate(pi/2);
  //2 * pi / segments
  double angle = tau / double(segments);
  //2 * pi * radius
  double segmentLength = (tau * radius) / segments;
  /*
    In order to properly fit segments one in the
    continuation of the other, we draw them with
    the angle of the middle of the circle surface
    they replace (not the start).
    This is how we get the least approximation error.
    Try setting this to zero, see what happens.
  */
  double midAngle = tau * double(0.5) / double(segments);
  t.rotate(midAngle);
  //<= so we make a loop, by overlapping the first and last segments.
  for(int ii = 0; ii <= segments; ++ii) {
    t.draw(segmentLength);
    t.rotate(angle);
  }
}

void drawSquare(Turtle t, float distance) {
  /*
    We assume the lower-left point of the square as the starting point,
    and the distance as the side length.
    (so: draw by moving forward and turning left)
  */
  t.draw(distance);

  t.rotate(pi/2);
  t.draw(distance);

  t.rotate(pi/2);
  t.draw(distance);

  t.rotate(pi/2);
  t.draw(distance);
}

void fractalKochCurve(Turtle t, float distance, int recursionsLeft = 1) {
  if(recursionsLeft > 0) {
    --recursionsLeft;
    distance /=3;
    
    //Draw straight forwards: '_'
    fractalKochCurve(t, distance, recursionsLeft); 
    t.move(distance);

    //Turn left: '_/'
    t.rotate(pi/3);
    fractalKochCurve(t, distance, recursionsLeft);
    t.move(distance);

    //Turn right: '_/\'
    t.rotate(- 2 * pi/3);
    fractalKochCurve(t, distance, recursionsLeft);
    t.move(distance);

    //Turn left: '_/\_'
    t.rotate(pi/3);
    fractalKochCurve(t, distance, recursionsLeft);
    //t.move(distance);
    /*
      ^ No need to move the equivalent distance,
      since no more segments are left.
    */
  } else {
    t.draw(distance);
  }
}

void Display1() {
  glColor3f(1, 0, 0);
  drawRecursionLevel();
  //Size of the fractal - radius of the circle circumscribing the starting triangle.
  double radius = 0.95;
  //Draw the circle containing the whole fractal.
  drawCircle(0, 0, radius, 36);
  //Start from the centre.
  Turtle t0(0, 0);
  Turtle t1 = t0;
  Turtle t2 = t0;
  //Rotate towards the 3 vertices of an equilateral triangle. (assuming we start from the centre).
  t0.rotate(0 * pi / 3);
  t1.rotate(2 * pi / 3);
  t2.rotate(4 * pi / 3);

  //Move onto the verices.
  t0.move(radius);
  t1.move(radius);
  t2.move(radius);
  
  //Rotate in the correct direction to draw edges from the vertices.
  t0.rotate(-pi/3 - pi/2);
  t1.rotate(-pi/3 - pi/2);
  t2.rotate(-pi/3 - pi/2);

  //Why sqrt(3)?
  fractalKochCurve(t0, sqrt(3) * radius, g_recursionCurrent);
  fractalKochCurve(t1, sqrt(3) * radius, g_recursionCurrent);
  fractalKochCurve(t2, sqrt(3) * radius, g_recursionCurrent);
}

void fractalBinaryTree(Turtle t, float distance, int recursionsLeft = 1) {
    if(recursionsLeft > 0) {
      --recursionsLeft;
      t.draw(distance);
      Turtle tLeft = t, tRight = t;
      tRight.rotate(-pi/4);
      tLeft.rotate(pi/4);
      fractalBinaryTree(tRight, distance/2, recursionsLeft);
      fractalBinaryTree(tLeft, distance/2, recursionsLeft);
    } else {
      t.draw(distance);
    }
}

void Display2() {
  glColor3f(1, 0, 0);
  drawRecursionLevel();
  Turtle t(0, -0.95);
  t.rotate(pi/2); //up
  fractalBinaryTree(t, 0.95, g_recursionCurrent);
}

void fractalSquare(Turtle t, float size, int recursionsLeft) {
    if(recursionsLeft > 0) {
      --recursionsLeft;
      //patratul initial
      t.rotate(pi/2);
        t.move(size);
      t.rotate(-pi/2);
        t.move(size);
      drawSquare(t, size);

      //ma mut in jurul patratului 
      t.move(4/3*size);//dreapta mijloc
      fractalSquare(t, size/3, recursionsLeft);
      t.rotate(-pi/2);//orientez in jos
      t.move(1/3*size);//mut originea patratului
      fractalSquare(t, size/3, recursionsLeft);
      t.rotate(-pi/2);//orientez in stanga
      t.move(1/3*size);
      fractalSquare(t, size/3, recursionsLeft);
      t.move(4/3*size);
      fractalSquare(t, size/3, recursionsLeft);
      t.rotate(-pi/2);
      t.move(1/3*size);
      fractalSquare(t, size/3, recursionsLeft);
      t.move(4/3*size);
      fractalSquare(t, size/3, recursionsLeft);
      t.rotate(-pi/2);
      t.move(1/3*size);
      fractalSquare(t, size/3, recursionsLeft);
      t.move(4/3*size);
      fractalSquare(t, size/3, recursionsLeft);
    }
  
}

void Display3() {
  //Draw the recursive-square fractal here.
  Turtle t(-0.95, -0.95);
  double x=1.90;
  drawSquare(t, x);

  fractalSquare(t, x/3, g_recursionCurrent);

  glColor3f(1, 0, 0);
  drawRecursionLevel();
  
}

//am substituit  conform unei regului  pentru fiecare valoare de A si B 
void fractalLSystemRecursive(Turtle &t, char symbol, float distance, int recursionsLeft) {
  if(recursionsLeft == 0) {
      if(symbol == 'A' ) {
          t.draw(distance);
        } else if(symbol == 'B') {
          t.draw(distance);
      } else if(symbol == '+') {
          t.rotate(pi / 3);
      } else if(symbol == '-') {
          t.rotate(-pi / 3); 
      }
  } else {
      if(symbol == 'A') {
          fractalLSystemRecursive(t, 'B', distance, recursionsLeft - 1);
          fractalLSystemRecursive(t, '-', distance, recursionsLeft - 1);
          fractalLSystemRecursive(t, 'A', distance, recursionsLeft - 1);
          fractalLSystemRecursive(t, '-', distance, recursionsLeft - 1);
          fractalLSystemRecursive(t, 'B', distance, recursionsLeft - 1);
      } else if(symbol == 'B') {
          fractalLSystemRecursive(t, 'A', distance, recursionsLeft - 1);
          fractalLSystemRecursive(t, '+', distance, recursionsLeft - 1);
          fractalLSystemRecursive(t, 'B', distance, recursionsLeft - 1);
          fractalLSystemRecursive(t, '+', distance, recursionsLeft - 1);
          fractalLSystemRecursive(t, 'A', distance, recursionsLeft - 1);
      } else if(symbol == '+' || symbol == '-') {
          fractalLSystemRecursive(t, symbol, distance, 0);
      }
  }
}

void Display4() {
  glColor3f(1, 0, 0);
  drawRecursionLevel();
  
  Turtle t(-0.95, -0.95);
  if (g_recursionCurrent % 2 != 0) {
    t.rotate(pi/3);
    fractalLSystemRecursive(t, 'A', 1.90 / pow(2, g_recursionCurrent), g_recursionCurrent);
  } else {
  fractalLSystemRecursive(t, 'A', 1.90 / pow(2, g_recursionCurrent), g_recursionCurrent);
}
}

template <typename FloatType>
class JF {
protected:
  //The x and y mathematical bounds of the fractal slice we're displaying.
  FloatType m_xmin, m_xmax, m_ymin, m_ymax;
  //The constant we're biasing the JF fractal with.
  std::complex<FloatType> m_c;
  //The radius around the origin we're using to detect divergence.
  FloatType m_maxRadius;
  //How many iterations we'll do to allow the number sequence to
  //exceed the limit.
  int m_maxIteration;

  virtual inline int test(std::complex<FloatType> z, std::complex<FloatType> c, double maxRadius = 2, int maxIteration = 50) {
    /*
      Compute the Julia-Fatou set in a point in 4D (x, y, a, b). Return the iterations *left*
      upon radius breach. So, a return value of 0 means estimated-divergence, other values
      mean speed of estimated convergence.
    */
    //We create a number sequence, and estimate its limit.
    for(int ii = maxIteration; ii > 0; --ii) {
      z = z * z + c;
      if(abs(z) > maxRadius)
	return(ii);
    }
    return 0;
  }
  
public:
  JF(FloatType xmin, FloatType xmax, FloatType ymin, FloatType ymax, FloatType a = 0, FloatType b = 0, FloatType maxRadius = 20, int maxIteration = 150):
    m_xmin(xmin),
    m_xmax(xmax),
    m_ymin(ymin),
    m_ymax(ymax),
    m_c(a, b),
    m_maxRadius(maxRadius),
    m_maxIteration(maxIteration) {
  }

  void draw(FloatType l, FloatType r, FloatType b, FloatType t, int samplePointsHorizontal, int samplePointsVertical) {
    /*
      Draw the current slice of the JF set onto the screen.
      Left, right, bottom, top, and the steps for each axis.
    */
    glPointSize(1);
    FloatType stepx = (m_xmax - m_xmin) / FloatType(samplePointsHorizontal);
    FloatType stepy = (m_ymax - m_ymin) / FloatType(samplePointsVertical);
    FloatType steph = (r      - l)      / FloatType(samplePointsHorizontal);
    FloatType stepv = (t      - b)      / FloatType(samplePointsVertical);
    int iterations;
    std::complex<FloatType> z;
    glBegin(GL_POINTS);
    /*
      We need to move both on screen pixels and in the mathematical plane -
      at the same time.
    */
    for(FloatType jj = 0, y = m_ymin, v = b; jj < samplePointsVertical; jj += 1, y += stepy, v += stepv) {
      z.imag(y);
      for(FloatType ii = 0, x = m_xmin, h = l; ii < samplePointsHorizontal; ii += 1, x += stepx, h += steph) {
	z.real(x);
	iterations = test(z, m_c, m_maxRadius, m_maxIteration);
	if(0 == iterations) {
	  glColor3f(1, 0, 0);
	  glVertex2d(h, v);	  
	}
      }
    }
    glEnd();
  }
};

void Display5() {
  glColor3f(1, 0, 0);
  drawJfConstants();
  float drawSize = 0.95;
  JF<double> jf(-2, 2, -2, 2, g_jfa, g_jfb);
  jf.draw(-drawSize, drawSize, -drawSize, drawSize, g_w, g_h);
}


template <typename FloatType>
  class MB: public JF<FloatType> {
  public:
    MB(FloatType xmin, FloatType xmax, FloatType ymin, FloatType ymax, FloatType a = 0, FloatType b = 0, FloatType maxRadius = 20, int maxIteration = 150):
      JF<FloatType>(xmin, xmax, ymin, ymax, a, b, maxRadius, maxIteration) {}
  
    void compute_mandelbrot(double left, double right, double top, double bottom) {
      const int MAX_ITERATIONS = 500;
      const int width = glutGet(GLUT_WINDOW_WIDTH);
      const int height = glutGet(GLUT_WINDOW_HEIGHT);
  
      glBegin(GL_POINTS);
      for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
          std::complex<double> c(left + (x * (right - left) / width), bottom + ((height - y) * (top - bottom) / height));
  
          std::complex<double> z(0.0, 0.0);
          int iterations = 0;
          while (abs(z) < 2.0 && iterations < MAX_ITERATIONS) {
            z = z * z + c;
            ++iterations;
          }
  
          float r, g, b;
  
          if (iterations == MAX_ITERATIONS) {
            glColor3f(0.0f, 0.0f, 0.0f);
          } else {
            float colorIndex = (iterations % 17) / 16.0f;
          
            if (colorIndex < 0.1) {
              glColor3f(0.0706f, 0.1294f, 0.9333f);
            } else if (colorIndex < 0.15) {
              glColor3f(0.165f, 0.263f, 0.835f);
            } else if (colorIndex < 0.2) {
              glColor3f(0.2f, 0.4f, 0.8f);
            } else if (colorIndex < 0.25) {
              glColor3f(0.267f, 0.533f, 0.729f);
            } else if (colorIndex < 0.3) {
              glColor3f(0.337f, 0.667f, 0.667f);
            } else if (colorIndex < 0.35) {
              glColor3f(0.4f, 0.796f, 0.6f);
            } else if (colorIndex < 0.4) {
              glColor3f(0.4f, 0.796f, 0.6f);
            } else if (colorIndex < 0.45) {
              glColor3f(0.471f, 0.933f, 0.529f);
            } else if (colorIndex < 0.5) {
              glColor3f(0.604f, 0.796f, 0.396f);
            } else if (colorIndex < 0.55) {
              glColor3f(0.667f, 0.667f, 0.333f);
            } else if (colorIndex < 0.6) {
              glColor3f(0.83f, 1.0f, 0.16f);
            } else if (colorIndex < 0.65) {
              glColor3f(0.733f, 0.533f, 0.263f);
            } else if (colorIndex < 0.7) {
              glColor3f(0.8f, 0.396f, 0.2f);
            } else {
              glColor3f(0.945f, 0.110f, 0.055f);
            }
          }
          
          glVertex2i(x, y);
        }
      }
      glEnd();
    }
  };
  


MB<double>* mandelbrot;
void Display6()
{
    glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    const int width = glutGet( GLUT_WINDOW_WIDTH );
    const int height = glutGet( GLUT_WINDOW_HEIGHT );
    glOrtho( 0, width, 0, height, -1, 1 );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    mandelbrot->compute_mandelbrot( -2.0, 2.0, -2.0, 2.0 );
    glutSwapBuffers();
}



void Display7() {
}

void Display8() {
}

void Display9() {
}

void Display10() {
}

void init(void) {
  glColor3f(1, 0, 0); //Just a starting default drawing colour.
  glClearColor(1.0,1.0,1.0,1.0);
  glLineWidth(1);
  glPointSize(1);
  //glPolygonMode(GL_FRONT, GL_LINE);
  //As we want pixel-perfect display for JF fractals, don't enable point smoothing.
  glEnable(GL_SMOOTH);
  //glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glHint(GL_NICEST, GL_POINT_SMOOTH_HINT);
  glHint(GL_NICEST, GL_LINE_SMOOTH_HINT);
  glHint(GL_NICEST, GL_POLYGON_SMOOTH_HINT);
  glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

  //Alpha-blending
  glEnable(GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

}

void Display(void) {
  // Clear the buffer. See init();
  glClear(GL_COLOR_BUFFER_BIT);

  switch(g_prevKey) {
  case '1':
    Display1();
    break;
  case '2':
    Display2();
    break;
  case '3':
    Display3();
    break;
  case '4':
    Display4();
    break;
  case '5':
    Display5();
    break;
  case '6':
    Display6();
    break;
  case '7':
    Display7();
    break;
  case '8':
    Display8();
    break;
  case '9':
    Display9();
    break;
  case '0':
    Display10();
    break;
  default:
    break;
  }
  glFlush();
}

void Reshape(int w, int h) {
  g_w = w;
  g_h = h;
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}
void KeyboardFunc(unsigned char key, int x, int y) {
  switch(key) {
  case 27: // escape
    exit(0);
    break;
  case '+':
    ++g_recursionCurrent;
    if(g_recursionCurrent > g_recursionMax)
      g_recursionCurrent = g_recursionMax;
    break;
  case '-':
    --g_recursionCurrent;
    if(g_recursionCurrent < 0)
      g_recursionCurrent = 0;
    break;
  case 'j':
    g_jfa -= 0.01;
    if(g_jfa < -2)
      g_jfa = -2;
    break;
  case 'l':
    g_jfa += 0.01;
    if(g_jfa > 2)
      g_jfa = 2;
    break;
  case 'k':
    g_jfb -= 0.01;
    if(g_jfb < -2)
      g_jfb = -2;
    break;
  case 'i':
    g_jfb += 0.01;
    if(g_jfb > 2)
      g_jfb = 2;
    break;
  default:
    //Only change the image if a 'special' key wasn't pressed.
    g_prevKey = key;
  }

  //The proper way to ask glut to redraw the window.
  glutPostRedisplay();
}

/*
  Callback upon mouse press or release.
  The button can be:
  GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, GLUT_RIGHT_BUTTON
  (and further for mousewheel and other mouse buttons)
  The state can be either GLUT_DOWN or  GLUT_UP, for
  a pressed or released button.
  (x, y) are the coordinates of the mouse.
*/
void MouseFunc(int button, int state, int x, int y) {
  std::cout<< "Mouse button ";
  std::cout<<( (button == GLUT_LEFT_BUTTON) ? "left" : ((button == GLUT_RIGHT_BUTTON) ? "right": "middle") ) << " ";
  std::cout<< ( (state == GLUT_DOWN) ? "pressed" : "released" );
  std::cout<< " at coordinates: " << x <<" x " << y << std::endl;
}

int main(int argc, char** argv) {
  glutInit(&argc, argv);
  glutInitWindowSize(g_w, g_h);
  glutInitWindowPosition(-1, -1);
  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGBA); 
  glutCreateWindow (argv[0]);
  init();
  glutReshapeFunc(Reshape);
  glutKeyboardFunc(KeyboardFunc);
  glutMouseFunc(MouseFunc);
  glutDisplayFunc(Display);
  //glutIdleFunc(Display);
  glutMainLoop();

  return 0;
}
