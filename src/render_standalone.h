
#ifndef __render_standalone_h__
#define __render_standalone_h__

#define USE_SOFTWARE_OPENGL 1

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "OpenGL/glu.h"
#else
#ifdef USE_SOFTWARE_OPENGL
#include "GL/osmesa.h"
#else
//#include "vtk_glew.h"
#endif
#include "GL/glu.h"
#endif

typedef double FieldData;

#endif // __render_standalone_h__

