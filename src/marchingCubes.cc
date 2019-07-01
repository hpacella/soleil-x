//
// Marching Cubes Example Program
// by Cory Bloyd (corysama@yahoo.com)
//
// A simple, portable and complete implementation of the Marching Cubes
// and Marching Tetrahedrons algorithms in a single source file.
// There are many ways that this code could be made faster, but the
// intent is for the code to be easy to understand.
//
// For a description of the algorithm go to
// http://astronomy.swin.edu.au/pbourke/modelling/polygonise/
//
// This code is public domain.
//

#include <stdio.h>
#include <math.h>
#include <iostream>

#include "render_standalone.h"
#include "renderImage.h"

//This program requires the OpenGL and GLUT libraries
// You can obtain them for free from http://www.opengl.org

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "OpenGL/glu.h"
#else
#ifdef USE_SOFTWARE_OPENGL
#include "GL/osmesa.h"
#else
#include "vtk_glew.h"

#endif
#include "GL/glu.h"
#endif

#define DEBUG 1


int gNumFluidX, gNumFluidY, gNumFluidZ;
const FieldData *gRho, *gPressure, *gTemperature;
const FieldData3 *gVelocity, *gCenterCoordinates;
FieldData *gDomainMin, *gDomainMax;
VisualizationField gVisualizationField;
FieldData* gIsosurfaceScale;
int gDrawnTriangles;


struct GLvector
{
  GLfloat fX;
  GLfloat fY;
  GLfloat fZ;
};

//These tables are used so that everything can be done in little loops that you can look at all at once
// rather than in pages and pages of unrolled code.

//a2fVertexOffset lists the positions, relative to vertex0, of each of the 8 vertices of a cube
static const GLfloat a2fVertexOffset[8][3] =
{
  {0.0, 0.0, 0.0},{1.0, 0.0, 0.0},{1.0, 1.0, 0.0},{0.0, 1.0, 0.0},
  {0.0, 0.0, 1.0},{1.0, 0.0, 1.0},{1.0, 1.0, 1.0},{0.0, 1.0, 1.0}
};

//a2iEdgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
static const GLint a2iEdgeConnection[12][2] =
{
  {0,1}, {1,2}, {2,3}, {3,0},
  {4,5}, {5,6}, {6,7}, {7,4},
  {0,4}, {1,5}, {2,6}, {3,7}
};

//a2fEdgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
static const GLfloat a2fEdgeDirection[12][3] =
{
  {1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
  {1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
  {0.0, 0.0, 1.0},{0.0, 0.0, 1.0},{ 0.0, 0.0, 1.0},{0.0,  0.0, 1.0}
};

//a2iTetrahedronEdgeConnection lists the index of the endpoint vertices for each of the 6 edges of the tetrahedron
static const GLint a2iTetrahedronEdgeConnection[6][2] =
{
  {0,1},  {1,2},  {2,0},  {0,3},  {1,3},  {2,3}
};

//a2iTetrahedronEdgeConnection lists the index of verticies from a cube
// that made up each of the six tetrahedrons within the cube
static const GLint a2iTetrahedronsInACube[6][4] =
{
  {0,5,1,6},
  {0,1,2,6},
  {0,2,3,6},
  {0,3,7,6},
  {0,7,4,6},
  {0,4,5,6},
};

static const GLfloat afAmbientWhite [] = {0.25, 0.25, 0.25, 1.00};
static const GLfloat afAmbientRed   [] = {0.25, 0.00, 0.00, 1.00};
static const GLfloat afAmbientGreen [] = {0.00, 0.25, 0.00, 1.00};
static const GLfloat afAmbientBlue  [] = {0.00, 0.00, 0.25, 1.00};
static const GLfloat afDiffuseWhite [] = {0.75, 0.75, 0.75, 1.00};
static const GLfloat afDiffuseRed   [] = {0.75, 0.00, 0.00, 1.00};
static const GLfloat afDiffuseGreen [] = {0.00, 0.75, 0.00, 1.00};
static const GLfloat afDiffuseBlue  [] = {0.00, 0.00, 0.75, 1.00};
static const GLfloat afSpecularWhite[] = {1.00, 1.00, 1.00, 1.00};
static const GLfloat afSpecularRed  [] = {1.00, 0.25, 0.25, 1.00};
static const GLfloat afSpecularGreen[] = {0.25, 1.00, 0.25, 1.00};
static const GLfloat afSpecularBlue [] = {0.25, 0.25, 1.00, 1.00};


GLenum    ePolygonMode = GL_FILL;
GLfloat   fTargetValue;
GLfloat   fTime = 0.0;
GLvector  sSourcePoint[3];
GLfloat v0[3], v1[3], v2[3], v3[3], v4[3], v5[3], v6[3], v7[3];
GLfloat gScaleX, gScaleY, gScaleZ;
GLfloat* gCameraLookAt;

GLfloat fSample(GLfloat fX, GLfloat fY, GLfloat fZ);

GLvoid vMarchingCubes();
GLvoid vMarchCube(GLfloat fX, GLfloat fY, GLfloat fZ, GLfloat fScale);


void initializeMarchingCubes(GLfloat lightPosition[4])
{
  GLfloat afPropertiesAmbient [] = {0.50, 0.50, 0.50, 1.00};
  GLfloat afPropertiesDiffuse [] = {0.75, 0.75, 0.75, 1.00};
  GLfloat afPropertiesSpecular[] = {1.00, 1.00, 1.00, 1.00};
  
#if 1
  glClearColor( 0, 0, 0, 1 );
#else
  glClearColor( 1, 1, 1, 1 );
#endif
  
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glPolygonMode(GL_FRONT_AND_BACK, ePolygonMode);
  
  glLightfv( GL_LIGHT0, GL_AMBIENT,  afPropertiesAmbient);
  glLightfv( GL_LIGHT0, GL_DIFFUSE,  afPropertiesDiffuse);
  glLightfv( GL_LIGHT0, GL_SPECULAR, afPropertiesSpecular);
  glLightfv( GL_LIGHT0, GL_POSITION, lightPosition);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
  printf("light position %f %f %f\n", lightPosition[0], lightPosition[1], lightPosition[2]);
  
  glEnable( GL_LIGHT0 );
  
  glMaterialfv(GL_BACK,  GL_AMBIENT,   afAmbientGreen);
  glMaterialfv(GL_BACK,  GL_DIFFUSE,   afDiffuseGreen);
  glMaterialfv(GL_FRONT, GL_AMBIENT,   afAmbientBlue);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   afDiffuseBlue);
  glMaterialfv(GL_FRONT, GL_SPECULAR,  afSpecularWhite);
  glMaterialf( GL_FRONT, GL_SHININESS, 25.0);
  
  
}




void vDrawScene(int numFluidX,
                int numFluidY,
                int numFluidZ,
                const FieldData* rho,
                const FieldData* pressure,
                const FieldData3* velocity,
                const FieldData3* centerCoordinates,
                const FieldData* temperature,
                FieldData domainMin[3],
                FieldData domainMax[3],
                VisualizationField visualizationField,
                FieldData targetValue,
                FieldData isosurfaceScale[2],
                GLfloat cameraLookAt[3])
{
  gNumFluidX = numFluidX;
  gNumFluidY = numFluidY;
  gNumFluidZ = numFluidZ;
  gRho = rho;
  gPressure = pressure;
  gVelocity = velocity;
  gCenterCoordinates = centerCoordinates;
  gTemperature = temperature;
  gDomainMin = domainMin;
  gDomainMax = domainMax;
  gVisualizationField = visualizationField;
  gIsosurfaceScale = isosurfaceScale;
  fTargetValue = targetValue;
  gCameraLookAt = cameraLookAt;
  
  gScaleX = (domainMax[0] - domainMin[0]) / numFluidX;
  gScaleY = (domainMax[1] - domainMin[1]) / numFluidY;
  gScaleZ = (domainMax[2] - domainMin[2]) / numFluidZ;
  
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  
  glPushMatrix();
  vMarchingCubes();
  glPopMatrix();
  
  glFinish();
  
}


//fGetOffset finds the approximate point of intersection of the surface
// between two points with the values fValue1 and fValue2
GLfloat fGetOffset(GLfloat fValue1, GLfloat fValue2, GLfloat fValueDesired)
{
  GLdouble fDelta = fValue2 - fValue1;
  
  if(fDelta == 0.0)
  {
    return 0.5;
  }
  return (fValueDesired - fValue1)/fDelta;
}

#if 0 // OLD

//vGetColor generates a color from a given position and normal of a point
GLvoid vGetColor(GLvector &rfColor, GLvector &rfPosition, GLvector &rfNormal)
{
  GLfloat fX = rfNormal.fX;
  GLfloat fY = rfNormal.fY;
  GLfloat fZ = rfNormal.fZ;
  rfColor.fX = (fX > 0.0 ? fX : 0.0) + (fY < 0.0 ? -0.5*fY : 0.0) + (fZ < 0.0 ? -0.5*fZ : 0.0);
  rfColor.fY = (fY > 0.0 ? fY : 0.0) + (fZ < 0.0 ? -0.5*fZ : 0.0) + (fX < 0.0 ? -0.5*fX : 0.0);
  rfColor.fZ = (fZ > 0.0 ? fZ : 0.0) + (fX < 0.0 ? -0.5*fX : 0.0) + (fY < 0.0 ? -0.5*fY : 0.0);
}

GLvoid vNormalizeVector(GLvector &rfVectorResult, GLvector &rfVectorSource)
{
  GLfloat fOldLength;
  GLfloat fScale;
  
  fOldLength = sqrtf( (rfVectorSource.fX * rfVectorSource.fX) +
                     (rfVectorSource.fY * rfVectorSource.fY) +
                     (rfVectorSource.fZ * rfVectorSource.fZ) );
  
  if(fOldLength == 0.0)
  {
    rfVectorResult.fX = rfVectorSource.fX;
    rfVectorResult.fY = rfVectorSource.fY;
    rfVectorResult.fZ = rfVectorSource.fZ;
  }
  else
  {
    fScale = 1.0/fOldLength;
    rfVectorResult.fX = rfVectorSource.fX*fScale;
    rfVectorResult.fY = rfVectorSource.fY*fScale;
    rfVectorResult.fZ = rfVectorSource.fZ*fScale;
  }
}



int findNearestCoordinate(GLfloat fx, GLfloat fy, GLfloat fz)
{
  int x_ = 0;
  int y_ = 0;
  int z_ = 0;
  FieldData best;
  int bestX = -99;
  int bestY = -99;
  int bestZ = -99;
  
  for(int i = 0; i < gNumFluidX; i++) {
    int index = x_ + gNumFluidX * y_ + gNumFluidX * gNumFluidY * z_;
    const FieldData3* coordinate = gCenterCoordinates + index;
    float distance = fabs(fx - coordinate->x[0]) + fabs(fy - coordinate->x[1]) + fabs(fz - coordinate->x[2]);
    if(bestX < 0 || best > distance) {
      bestX = x_;
      best = distance;
    }
    x_++;
  }
  
  for(int i = 0; i < gNumFluidY; i++) {
    int index = x_ + gNumFluidX * y_ + gNumFluidX * gNumFluidY * z_;
    const FieldData3* coordinate = gCenterCoordinates + index;
    float distance = fabs(fx - coordinate->x[0]) + fabs(fy - coordinate->x[1]) + fabs(fz - coordinate->x[2]);
    if(bestY < 0 || best > distance) {
      bestY = y_;
      best = distance;
    }
    y_++;
  }
  
  for(int i = 0; i < gNumFluidZ; i++) {
    int index = x_ + gNumFluidX * y_ + gNumFluidX * gNumFluidY * z_;
    const FieldData3* coordinate = gCenterCoordinates + index;
    float distance = fabs(fx - coordinate->x[0]) + fabs(fy - coordinate->x[1]) + fabs(fz - coordinate->x[2]);
    if(bestZ < 0 || best > distance) {
      bestZ = z_;
      best = distance;
    }
    z_++;
  }
  int bestIndex = bestX + gNumFluidX * bestY + gNumFluidX * gNumFluidY * bestZ;
  return bestIndex;
}



GLfloat fSample(GLfloat fX, GLfloat fY, GLfloat fZ)
{
  int index = findNearestCoordinate(fX, fY, fZ);
  const FieldData* data;
  switch(gVisualizationField) {
    case rhoField:
      data = gRho;
      break;
    case pressureField:
      data = gPressure;
      break;
    case temperatureField:
      data = gTemperature;
      break;
    default:
      std::cerr << "invalid visualization field " << gVisualizationField << std::endl;
      std::cerr << "fields " << rhoField << " " << pressureField << " " << temperatureField << std::endl;
      data = gRho;
      break;
  }
  return data[index];
}



//vGetNormal() finds the gradient of the scalar field at a point
//This gradient can be used as a very accurate vertx normal for lighting calculations
GLvoid vGetNormal(GLvector &rfNormal, GLfloat fX, GLfloat fY, GLfloat fZ)
{
  rfNormal.fX = fSample(fX-0.01, fY, fZ) - fSample(fX+0.01, fY, fZ);
  rfNormal.fY = fSample(fX, fY-0.01, fZ) - fSample(fX, fY+0.01, fZ);
  rfNormal.fZ = fSample(fX, fY, fZ-0.01) - fSample(fX, fY, fZ+0.01);
  vNormalizeVector(rfNormal, rfNormal);
}

#endif // OLD


// draw a cube, for when all vertices match the isosurface value
GLvoid drawCube(int index)
{
  const FieldData3* coordinates = gCenterCoordinates + index;
  GLfloat fX = coordinates->x[0];
  GLfloat fY = coordinates->x[1];
  GLfloat fZ = coordinates->x[2];
  
#if 1
  std::cout << "draw cube at " << fX << " " << fY << " " << fZ << std::endl;
#endif

  GLfloat sx = gScaleX * 0.5;
  GLfloat sy = gScaleY * 0.5;
  GLfloat sz = gScaleZ * 0.5;

  GLfloat v0[] = {  sx,  sy,  sz };
  GLfloat v1[] = { -sx,  sy,  sz };
  GLfloat v2[] = { -sx, -sy,  sz };
  GLfloat v3[] = {  sx, -sy,  sz };
  GLfloat v4[] = {  sx, -sy, -sz };
  GLfloat v5[] = {  sx,  sy, -sz };
  GLfloat v6[] = { -sx,  sy, -sz };
  GLfloat v7[] = { -sx, -sy, -sz };

  glPushMatrix();
  glTranslatef(fX, fY, fZ);
  glBegin(GL_TRIANGLES);

  // cube vertices follow example at http://www.songho.ca/opengl/gl_vertexarray.html

  // front face
  GLfloat n0[] = { 0, 0, 1 };
  glNormal3fv(n0);

    glVertex3fv(v0);    // v0-v1-v2
    glVertex3fv(v3);
    glVertex3fv(v2);

    glVertex3fv(v2);    // v2-v3-v0
    glVertex3fv(v1);
    glVertex3fv(v0);

    // right face =================
  GLfloat n1[] = { 1, 0, 0 };
  glNormal3fv(n1);

    glVertex3fv(v0);
    glVertex3fv(v5);
    glVertex3fv(v4);

    glVertex3fv(v4);
    glVertex3fv(v3);
    glVertex3fv(v0);

    // top face ===================
  GLfloat n2[] = { 0, 1, 0 };
  glNormal3fv(n2);

    glVertex3fv(v5);    
    glVertex3fv(v0);
    glVertex3fv(v1);

    glVertex3fv(v1);   
    glVertex3fv(v6);
    glVertex3fv(v5);

  // rear face
  GLfloat n3[] = { 0, 0, -1 };
  glNormal3fv(n3);

    glVertex3fv(v4);
    glVertex3fv(v5);
    glVertex3fv(v6);

    glVertex3fv(v6);
    glVertex3fv(v7);
    glVertex3fv(v4);

  // left face
  GLfloat n4[] = { -1, 0, 0 };
  glNormal3fv(n4);

    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v7);

    glVertex3fv(v7);
    glVertex3fv(v6);
    glVertex3fv(v1);

  // bottom face
  GLfloat n5[] = { 0, -1, 0 };
  glNormal3fv(n5);

    glVertex3fv(v3);
    glVertex3fv(v4);
    glVertex3fv(v7);

    glVertex3fv(v7);
    glVertex3fv(v2);
    glVertex3fv(v3);

  glEnd();
  glPopMatrix();

  gDrawnTriangles += 12;

}

#if 0

//vMarchCube1 performs the Marching Cubes algorithm on a single cube
GLvoid vMarchCubeOLD(GLfloat fX, GLfloat fY, GLfloat fZ, GLfloat scaleX, GLfloat scaleY, GLfloat scaleZ)
{

  extern GLint aiCubeEdgeFlags[256];
  extern GLint a2iTriangleConnectionTable[256][16];

#if DEBUG
  std::cout << "vMarchCube " << fX << " " << fY << " " << fZ << std::endl;
#endif
  
  GLint iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
  GLfloat fOffset;
  GLvector sColor;
  GLfloat afCubeValue[8];
  GLvector asEdgeVertex[12];
  GLvector asEdgeNorm[12];
  
  //Make a local copy of the values at the cube's corners
  for(iVertex = 0; iVertex < 8; iVertex++)
  {
    afCubeValue[iVertex] = fSample(fX + a2fVertexOffset[iVertex][0]*scaleX,
                                   fY + a2fVertexOffset[iVertex][1]*scaleY,
                                   fZ + a2fVertexOffset[iVertex][2]*scaleZ);
#if DEBUG
    std::cout << "vertex " << iVertex << " at " << (fX + a2fVertexOffset[iVertex][0]*scaleX) << " " << (fY + a2fVertexOffset[iVertex][1]*scaleY) << " " << (fZ + a2fVertexOffset[iVertex][2]*scaleZ) << std::endl;
    std::cout << "sample " << afCubeValue[iVertex] << std::endl;
#endif
  }
  
  //Find which vertices are inside of the surface and which are outside
  iFlagIndex = 0;
  unsigned equivalences = 0;
  for(iVertexTest = 0; iVertexTest < 8; iVertexTest++)
  {
#if DEBUG
    std::cout << "vertex value " << afCubeValue[iVertexTest] << " <=? " << fTargetValue << " ";
#endif
    if(afCubeValue[iVertexTest] == fTargetValue) equivalences++;
    if(afCubeValue[iVertexTest] <= fTargetValue) {
      iFlagIndex |= 1<<iVertexTest;
#if DEBUG
      std::cout << "YES iFlagIndex " << iFlagIndex;
#endif
    }
#if DEBUG
    std::cout << std::endl;
#endif
  }
  
  //Find which edges are intersected by the surface
  iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

#if DEBUG
  std::cout << "iEdgeFlags " << iEdgeFlags << std::endl;
#endif
  
  //If the cube is entirely inside or outside of the surface, then there will be no intersections
  if(iEdgeFlags == 0)
  {
    if(equivalences != 8) {
#if DEBUG
      std::cout << "cube is empty" << std::endl;
#endif
      return;
    }
#if DEBUG
    std::cout << "all vertices were equal" << std::endl;
#endif
    drawCube(fX, fY, fZ, scaleX, scaleY, scaleZ);
    return;
  }
  
  //Find the point of intersection of the surface with each edge
  //Then find the normal to the surface at those points
  for(iEdge = 0; iEdge < 12; iEdge++)
  {
    //if there is an intersection on this edge
    if(iEdgeFlags & (1<<iEdge))
    {
      fOffset = fGetOffset(afCubeValue[ a2iEdgeConnection[iEdge][0] ],
                           afCubeValue[ a2iEdgeConnection[iEdge][1] ], fTargetValue);

#if DEBUG
      std::cout << "iEdge " << iEdge << " fOffset " << fOffset << std::endl;
#endif
      
      asEdgeVertex[iEdge].fX = fX + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][0]  +  fOffset * a2fEdgeDirection[iEdge][0]) * scaleX;
      asEdgeVertex[iEdge].fY = fY + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][1]  +  fOffset * a2fEdgeDirection[iEdge][1]) * scaleY;
      asEdgeVertex[iEdge].fZ = fZ + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][2]  +  fOffset * a2fEdgeDirection[iEdge][2]) * scaleZ;
      
      vGetNormal(asEdgeNorm[iEdge], asEdgeVertex[iEdge].fX, asEdgeVertex[iEdge].fY, asEdgeVertex[iEdge].fZ);
#if DEBUG
      std::cout << "asEdgeVertex " << asEdgeVertex[iEdge].fX << " " << asEdgeVertex[iEdge].fY << " " << asEdgeVertex[iEdge].fZ << std::endl;
#endif
    }
  }
  
  
  //Draw the triangles that were found.  There can be up to five per cube
  glBegin(GL_TRIANGLES);
  for(iTriangle = 0; iTriangle < 5; iTriangle++)
  {
    if(a2iTriangleConnectionTable[iFlagIndex][3*iTriangle] < 0) {
      break;
    }
    
    for(iCorner = 0; iCorner < 3; iCorner++)
    {
      iVertex = a2iTriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];
      
      vGetColor(sColor, asEdgeVertex[iVertex], asEdgeNorm[iVertex]);
      glMaterialfv(GL_FRONT, GL_AMBIENT,   color);
      glMaterialfv(GL_FRONT, GL_DIFFUSE,   color);
      glMaterialfv(GL_FRONT, GL_SPECULAR,  afSpecularWhite);
      glColor3f(sColor.fX, sColor.fY, sColor.fZ);
      glNormal3f(asEdgeNorm[iVertex].fX, asEdgeNorm[iVertex].fY,   asEdgeNorm[iVertex].fZ);
      glVertex3f(asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);

#if DEBUG
      printf("tri %g %g %g\n", asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);
#endif

      gDrawnTriangles++;
    }
  }  
  glEnd();
}

#endif


// handles wraparound boundary conditions
GLint getOffsetIndex(int index, int offsetX, int offsetY, int offsetZ) {
  int ix = index / (gNumFluidZ * gNumFluidY);
  int tmp = index - ix * gNumFluidZ * gNumFluidY;
  int iy = tmp / gNumFluidZ;
  tmp -= iy * gNumFluidZ;
  int iz = tmp;
  int nextX = ix + offsetX;
  if(nextX < 0) nextX = gNumFluidX - 1;
  if(nextX >= gNumFluidX) nextX = 0;
  int nextY = iy + offsetY;
  if(nextY < 0) nextY = gNumFluidY - 1;
  if(nextY >= gNumFluidY) nextY = 0;
  int nextZ = iz + offsetZ;
  if(nextZ < 0) nextZ = gNumFluidZ - 1;
  if(nextZ >= gNumFluidZ) nextZ = 0;
  int nextIndex = nextZ + nextY * gNumFluidZ + nextX * gNumFluidZ * gNumFluidY;
#if DEBUG
  std::cout << "cell " << index << " (" << ix << "," << iy << "," << iz << ")\tneighbor " << nextIndex << " (" << nextX << "," << nextY << "," << nextZ << ")" << std::endl;
#endif
  return nextIndex;
}


// gather the indices of points surrounding a cell corner
GLvoid getNeighborIndices(int index, int offsetX, int offsetY, int offsetZ, int indices[])
{
  unsigned indexCount = 0;
  for(unsigned ix = 0; ix < 2; ++ix) {
    int offX = ix ? offsetX : 0;
    for(unsigned iy = 0; iy < 2; ++iy) {
      int offY = iy ? offsetY : 0;
      for(unsigned iz = 0; iz < 2; ++iz) {
        int offZ = iz ? offsetZ : 0;
        indices[indexCount++] = getOffsetIndex(index, offX, offY, offZ);
      }
    }
  }
}


GLfloat getFluidSample(int index) {
  const FieldData* data;
  switch(gVisualizationField) {
    case rhoField:
      data = gRho;
      break;
    case pressureField:
      data = gPressure;
      break;
    case temperatureField:
      data = gTemperature;
      break;
    default:
      std::cerr << "invalid visualization field " << gVisualizationField << std::endl;
      std::cerr << "fields " << rhoField << " " << pressureField << " " << temperatureField << std::endl;
      data = gRho;
      break;
  }
#if DEBUG
  std::cout << "fluid sample at " << index << " = " << data[index] << std::endl;
#endif
  return data[index];
}



//vMarchCube performs the Marching Cubes algorithm on a single cube
GLvoid vMarchCube(int index)
{
  extern GLint aiCubeEdgeFlags[256];
  extern GLint a2iTriangleConnectionTable[256][16];
  
#if DEBUG
  std::cout << "------------------------------" << std::endl;
  std::cout << "vMarchCube " << index << std::endl;
#endif

  // Make a local copy of the values at the cube's corners
  const unsigned numCorners = 8;
  GLfloat afCubeValue[numCorners];
  for(unsigned corner = 0; corner < numCorners; ++corner) {
    int offsetX = (corner & 1) ? 1 : -1;
    int offsetY = (corner & 2) ? 1 : -1;
    int offsetZ = (corner & 4) ? 1 : -1;
    int neighborIndices[numCorners];
    getNeighborIndices(index, offsetX, offsetY, offsetZ, neighborIndices);
    GLfloat cornerSum = 0;
    for(unsigned i = 0; i < numCorners; ++i) {
      cornerSum += getFluidSample(neighborIndices[i]);
    }
    afCubeValue[corner] = cornerSum / numCorners;
#if DEBUG
  std::cout << "cube corner " << corner << " has value " << afCubeValue[corner] << std::endl;
#endif
  }
  
  //Find which vertices are inside of the surface and which are outside
  int iFlagIndex = 0;
  unsigned equivalences = 0;
  for(unsigned iVertexTest = 0; iVertexTest < numCorners; iVertexTest++)
  {
    if(afCubeValue[iVertexTest] == fTargetValue) equivalences++;
    if(afCubeValue[iVertexTest] <= fTargetValue) {
      iFlagIndex |= 1 << iVertexTest;
    }
  }

  //Find which edges are intersected by the surface
  int iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];
  
#if DEBUG
  std::cout << "iEdgeFlags " << iEdgeFlags << std::endl;
#endif

  //If the cube is entirely inside or outside of the surface, then there will be no intersections
  if(iEdgeFlags == 0)
  {
    if(equivalences == numCorners) {
      drawCube(index);
    }
    return;
  }

  //Find the point of intersection of the surface with each edge
  const unsigned numEdges = 12;
  GLvector asEdgeVertex[numEdges];
  GLvector asEdgeNorm[numEdges];
  const FieldData3* coordinate = gCenterCoordinates + index;
  GLfloat fX = coordinate->x[0];
  GLfloat fY = coordinate->x[1];
  GLfloat fZ = coordinate->x[2];
  
  for(unsigned iEdge = 0; iEdge < numEdges; iEdge++)
  {
    //if there is an intersection on this edge
    if(iEdgeFlags & (1 << iEdge))
    {
      GLfloat fOffset = fGetOffset(afCubeValue[ a2iEdgeConnection[iEdge][0] ],
                           afCubeValue[ a2iEdgeConnection[iEdge][1] ], fTargetValue);
      asEdgeVertex[iEdge].fX = fX + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][0]  +  fOffset * a2fEdgeDirection[iEdge][0]) * gScaleX;
      asEdgeVertex[iEdge].fY = fY + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][1]  +  fOffset * a2fEdgeDirection[iEdge][1]) * gScaleY;
      asEdgeVertex[iEdge].fZ = fZ + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][2]  +  fOffset * a2fEdgeDirection[iEdge][2]) * gScaleZ;

#if 1//DEBUG ONLY
asEdgeVertex[iEdge].fX *= 20.0;
asEdgeVertex[iEdge].fY *= 20.0;
asEdgeVertex[iEdge].fZ *= 20.0;
#endif

#if DEBUG
      std::cout << "vertex " << asEdgeVertex[iEdge].fX << " " << asEdgeVertex[iEdge].fY << " " << asEdgeVertex[iEdge].fZ << std::endl;
#endif
    }
  }
  
  unsigned numVertices = 0;
  GLvector vertex[3];
  unsigned edge[3];
  for(unsigned iEdge = 0; iEdge < numEdges; iEdge++)
  {
    //if there is an intersection on this edge
    if(iEdgeFlags & (1 << iEdge))
    {
      GLvector v = { asEdgeVertex[iEdge].fX, asEdgeVertex[iEdge].fY, asEdgeVertex[iEdge].fX };
      edge[numVertices % 3] = iEdge;
      vertex[numVertices++ % 3] = v;
      GLvector v01, v02;
      v01.fX = vertex[1].fX - vertex[0].fX;
      v01.fY = vertex[1].fY - vertex[0].fY;
      v01.fZ = vertex[1].fZ - vertex[0].fZ;
      v02.fX = vertex[2].fX - vertex[0].fX;
      v02.fY = vertex[2].fY - vertex[0].fY;
      v02.fZ = vertex[2].fZ - vertex[0].fZ;
      GLvector normal;
      normal.fX = v01.fY * v02.fZ - v01.fZ * v02.fY;
      normal.fY = v01.fZ * v02.fX - v01.fX * v02.fZ;
      normal.fZ = v01.fX * v02.fY - v01.fY * v02.fX;
      GLfloat dot = normal.fX * gCameraLookAt[0] + normal.fY * gCameraLookAt[1] + normal.fZ * gCameraLookAt[2];
      GLfloat dotNegative = -normal.fX * gCameraLookAt[0] - normal.fY * gCameraLookAt[1] - normal.fZ * gCameraLookAt[2];
      GLfloat norm = sqrt(normal.fX * normal.fX + normal.fY * normal.fY + normal.fZ * normal.fZ);
      normal.fX /= norm;
      normal.fY /= norm;
      normal.fZ /= norm;
#if DEBUG
      std::cout << "dot " << dot << " dotNegative " << dotNegative << " normal " << normal.fX << " " << normal.fY << " " << normal.fZ << std::endl;
#endif
      if(dot > dotNegative) {
        normal.fX *= -1.0;
        normal.fY *= -1.0;
        normal.fZ *= -1.0;
      }
      
    }
  }
  
  //Draw the triangles that were found.  There can be up to five per cube
  const unsigned maxTriangles = 5;
  glBegin(GL_TRIANGLES);
  for(unsigned iTriangle = 0; iTriangle < maxTriangles; iTriangle++)
  {
    if(a2iTriangleConnectionTable[iFlagIndex][3*iTriangle] < 0) {
      break;
    }
    const unsigned numTriVertices = 3;
    for(unsigned iCorner = 0; iCorner < numTriVertices; iCorner++)
    {
      int iVertex = a2iTriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];
      
      // get color by scale
      FieldData isosurfaceScalar = getFluidSample(index);
      GLfloat color[4];
      scaledTemperatureToColor(isosurfaceScalar, color, gIsosurfaceScale);
      glColor4fv(color);
      glMaterialfv(GL_FRONT, GL_AMBIENT,   color);
      glMaterialfv(GL_FRONT, GL_DIFFUSE,   color);
      glMaterialfv(GL_FRONT, GL_SPECULAR,  afSpecularWhite);
      
      glNormal3f(asEdgeNorm[iVertex].fX, asEdgeNorm[iVertex].fY,   asEdgeNorm[iVertex].fZ);
      glVertex3f(asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);
      gDrawnTriangles++;
    }

  }
  glEnd();
}



#if 1
void drawParticle(GLUquadricObj* qobj, const FieldData3* position, FieldData density, FieldData particleTemperature, float particleSize);
#endif

//vMarchingCubes iterates over the entire dataset, calling vMarchCube on each cube
GLvoid vMarchingCubes()
{
  gDrawnTriangles = 0;
  GLint iX, iY, iZ;

  for(iX = 0; iX < gNumFluidX; iX++)
    for(iY = 0; iY < gNumFluidY; iY++)
      for(iZ = 0; iZ < gNumFluidZ; iZ++)
      {
        int index = iX + gNumFluidX * iY + gNumFluidX * gNumFluidY * iZ;
#define RENDER_CELL_CENTERS 0
#if RENDER_CELL_CENTERS
        GLUquadricObj *qobj = gluNewQuadric();
        double scale = 1.0e-3;
        const FieldData3* coordinate = gCenterCoordinates + index;
        drawParticle(qobj, coordinate, 500, 300, scale);
        gluDeleteQuadric(qobj);
#else
        vMarchCube(index);
return;
#endif
      }
  std::cout << "drew " << gDrawnTriangles << " triangles" << std::endl;
}


// For any edge, if one vertex is inside of the surface and the other is outside of the surface
//  then the edge intersects the surface
// For each of the 4 vertices of the tetrahedron can be two possible states : either inside or outside of the surface
// For any tetrahedron the are 2^4=16 possible sets of vertex states
// This table lists the edges intersected by the surface for all 16 possible vertex states
// There are 6 edges.  For each entry in the table, if edge #n is intersected, then bit #n is set to 1

GLint aiTetrahedronEdgeFlags[16]=
{
  0x00, 0x0d, 0x13, 0x1e, 0x26, 0x2b, 0x35, 0x38, 0x38, 0x35, 0x2b, 0x26, 0x1e, 0x13, 0x0d, 0x00,
};


// For each of the possible vertex states listed in aiTetrahedronEdgeFlags there is a specific triangulation
// of the edge intersection points.  a2iTetrahedronTriangles lists all of them in the form of
// 0-2 edge triples with the list terminated by the invalid value -1.
//
// I generated this table by hand

GLint a2iTetrahedronTriangles[16][7] =
{
  {-1, -1, -1, -1, -1, -1, -1},
  { 0,  3,  2, -1, -1, -1, -1},
  { 0,  1,  4, -1, -1, -1, -1},
  { 1,  4,  2,  2,  4,  3, -1},
  
  { 1,  2,  5, -1, -1, -1, -1},
  { 0,  3,  5,  0,  5,  1, -1},
  { 0,  2,  5,  0,  5,  4, -1},
  { 5,  4,  3, -1, -1, -1, -1},
  
  { 3,  4,  5, -1, -1, -1, -1},
  { 4,  5,  0,  5,  2,  0, -1},
  { 1,  5,  0,  5,  3,  0, -1},
  { 5,  2,  1, -1, -1, -1, -1},
  
  { 3,  4,  2,  2,  4,  1, -1},
  { 4,  1,  0, -1, -1, -1, -1},
  { 2,  3,  0, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1},
};

// For any edge, if one vertex is inside of the surface and the other is outside of the surface
//  then the edge intersects the surface
// For each of the 8 vertices of the cube can be two possible states : either inside or outside of the surface
// For any cube the are 2^8=256 possible sets of vertex states
// This table lists the edges intersected by the surface for all 256 possible vertex states
// There are 12 edges.  For each entry in the table, if edge #n is intersected, then bit #n is set to 1

GLint aiCubeEdgeFlags[256]=
{
  0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
};

//  For each of the possible vertex states listed in aiCubeEdgeFlags there is a specific triangulation
//  of the edge intersection points.  a2iTriangleConnectionTable lists all of them in the form of
//  0-5 edge triples with the list terminated by the invalid value -1.
//  For example: a2iTriangleConnectionTable[3] list the 2 triangles formed when corner[0]
//  and corner[1] are inside of the surface, but the rest of the cube is not.
//
//  I found this table in an example program someone wrote long ago.  It was probably generated by hand

GLint a2iTriangleConnectionTable[256][16] =
{
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
  {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
  {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
  {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
  {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
  {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
  {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
  {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
  {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
  {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
  {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
  {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
  {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
  {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
  {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
  {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
  {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
  {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
  {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
  {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
  {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
  {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
  {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
  {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
  {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
  {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
  {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
  {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
  {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
  {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
  {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
  {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
  {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
  {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
  {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
  {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
  {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
  {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
  {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
  {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
  {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
  {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
  {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
  {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
  {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
  {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
  {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
  {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
  {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
  {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
  {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
  {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
  {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
  {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
  {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
  {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
  {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
  {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
  {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
  {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
  {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
  {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
  {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
  {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
  {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
  {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
  {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
  {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
  {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
  {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
  {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
  {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
  {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
  {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
  {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
  {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
  {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
  {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
  {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
  {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
  {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
  {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
  {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
  {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
  {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
  {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
  {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
  {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
  {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
  {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
  {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
  {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
  {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
  {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
  {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
  {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
  {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
  {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
  {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
  {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
  {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
  {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
  {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
  {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
  {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
  {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
  {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
  {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
  {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
  {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
  {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
  {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
  {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
  {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
  {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
  {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
  {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
  {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
  {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
  {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
  {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
  {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
  {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
  {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
  {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
  {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
  {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
  {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};


