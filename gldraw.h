#ifndef GLDRAW_H
#define GLDRAW_H

#include <vector>
#include <cmath>
#include <math.h>
#include <locale.h>
#include <string.h>

/*
#include <QGLWidget>
#include <QtOpenGL>
#include <QGLFunctions>
#include <QOpenGLBuffer>
#include <QGLShaderProgram>
#include <QOpenGLFunctions>
#include <QBasicTimer>
#include <QDebug>
*/

#include <fstream>
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <iomanip>

#include "point.h"
#include "line.h"

#define MU                  0.0
#define SIGMA               1.0

#define UPPERBOUND          1.0
#define LOWERBOUND          -1.0

#define BALLNUM             100
#define RADIUS              0.01
#define LINEWIDTH           0.01

#define B_VERTICES          0
#define B_COLOR             1
#define B_CENTER            2
#define B_RADIUS            3
#define B_LINES             4
#define B_COOR              5
#define B_POLYGONS          6
#define B_TEX               7
#define B_LATTRI            8

#define C_COM_COOR          0
#define C_POS_COOR          1

#define T_TEX               0
#define T_POS_1             1
#define T_POS_2             2

#define F_POS               0
#define F_TEX               1

#define PI                  3.1415926535897
#define ATTRACTION          0.1
#define ELECTRO             1

#define WIDTH               800
#define HEIGHT              800

#define EPS                 1.0e-6
#define COM                 0.68
#define SMOOTH              0.5

#define SEGMENT             1

#define S_INITIAL           0.0004       // init. distance to move points
#define P_INITIAL           1            // init. subdivision number
#define P_RATE              2            // subdivision rate increase
#define CYCLE               6           // number of cycles to perform
#define I_INITIAL           50           // init. number of iterations for cycle
#define I_RATE              0.6666667    // rate at which iteration number decreases i.e. 2/3

using namespace std;

typedef float GLfloat;
typedef bool GLboolean;
typedef unsigned int GLuint;
typedef int GLint;

class GLDraw
{
    //Q_OBJECT

public:

    explicit GLDraw();

    ~GLDraw(){};

    void bundleEdges();

    string                  m_sDataNode;
    string                  m_sDataEdge;

protected:

    GLfloat CDF(int i, int n);

    GLfloat gaussian(GLfloat X);

    GLfloat getRandomArbitrary(GLfloat min, GLfloat max);

    GLfloat calcDistance(GLfloat *p1, GLfloat *p2);

    GLfloat vectorDotProduct(GLfloat *p1, GLfloat *p2);

    GLfloat compatibilityScore(Line *l1, Line *l2);

    GLfloat angleCompatibility(Line *l1, Line *l2);

    GLfloat scaleCompatibility(Line *l1, Line *l2);

    GLfloat positionCompatibility(Line *l1, Line *l2);

    GLfloat visibilityCompatibility(Line *l1, Line *l2);

    GLfloat edgeVisibility(Line *l1, Line *l2);

    float safeAcos(float x);

    float length(float x, float y);

    bool areCompatible(Line *l1, Line *L2);

    void computeCompatibility();

    void edgeAsVector(Line *l1, GLfloat *result);

    void norVector(GLfloat *n);

    void projectPointOnLine(GLfloat *Qsource, GLfloat *Psource, GLfloat *Psink,
                            GLfloat *I);

    void edge_midpoint(GLfloat *p1, GLfloat *p2, GLfloat *midPoint);

    void gaussianSmooth(int P);

    void calAttraction(GLfloat *p1, GLfloat *p2, GLfloat *force, GLfloat c);

    void calElectrostatic(GLfloat *p1, GLfloat *p2, GLfloat *force, GLfloat l);

    void addForce(GLfloat *p1, GLfloat *p2, GLfloat *force);

    void initCompatibilityDistance();

    void calcCompatibilityDistance();

    void saveTexBuffer();

    void initBuffer();

    void initVertex();

    void initShaders1();

    void initShaders2();

    void initShaders3();

    void initShaders4();

    void initShaders5();

    void initPosShaders();

    void initSubDShaders();

    void initAdjMatrix();

    void initDistance();

    void initVBO();

    void initTextureFrameBuffer();

    void initComFrameBuffer();

    void initPosFrameBuffer(int P);

    void initPosBuffer();

    void initializeGL();

    void paintGL();

    void resizeGL(int w, int h);

    void renderToFrameBuffer();

    void renderToScreen();

    void swapPosFrameBuffer();

    void swapPosBuffer();

    void updatePositionFrameBuffer(int P, float S);

    void updateDivisionFrameBuffer(int P);

    void updateEdgeDivisions(int P);

    void updateLength(int P);

    void updateForcePosition(int P, float S);

    void updateBuffer();

    void saveBuffer();

    void readBuffer();

    double read_timer();

    void turnoffTexture();



private:

    vector<Point*>          m_vPoints;

    Point                   *m_pPoint;

    vector<Line*>           m_vLines;

    Line                    *m_pLine;

    GLint                   *m_pTexLimit;

    GLfloat                 *m_pPosition;

    GLfloat                 *m_pPixels;

    GLfloat                 *m_pMidPoint;

    GLfloat                 *m_pAdjMatrix;

    GLfloat                 *m_pAdjMatrix_output;

    GLfloat                 *m_pVert;

    GLfloat                 *m_pColor;

    GLfloat                 *m_pGL_dat_lines;

    GLfloat                 *m_pGL_dat_points;

    GLfloat                 *m_pGL_dat_colors;

    GLfloat                 *m_pGL_dat_centers;

    GLfloat                 *m_pGL_dat_radius;

    GLfloat                 *m_pGL_dat_attrs;

    GLfloat                 *m_pGL_dat_polygons;

    GLfloat                 *m_pGL_dat_lattri;

    GLfloat                 *m_pSubPoints_1;

    GLfloat                 *m_pSubPoints_2;

    GLfloat                 *m_pCompatibility;

    GLfloat                 *m_pSwapIn;

    GLfloat                 *m_pSwapOut;

    GLfloat                 *m_pSmoothPoints;

    GLfloat                 *m_pOriginPoints;

    GLfloat                 m_fLineWidth;

    GLfloat                 m_fMaxLength;

    GLint                   m_nLineNum;

    GLint                   m_nSegment;

    GLint                   m_nWidth;

    GLint                   m_nHeight;

    GLint                   m_nBallNum;

    GLint                   m_nMax;

    GLuint                  m_pVertexBuffer[9];

    GLuint                  m_pCalBuffer[2];

    GLuint                  m_pTexBuffer[3];

    GLuint                  m_pFrameBuffer[2];

    GLuint                  m_nSwapPosOff;

    GLuint                  m_nSwapPosRender;

    GLboolean               m_bLineWidthSet;

    GLboolean               m_bSegmentSet;

    GLboolean               m_bVertices;

    GLboolean               m_bEdges;

    GLboolean               m_bPolygons;

    GLboolean               m_bCurveShaders;

    GLboolean               m_bNodeShaders;

    GLboolean               m_bPos_2;

    GLboolean               m_bDivEven;

    GLboolean               m_bSmooth;

    GLboolean               m_bData1;

    GLboolean               m_bUpdate;


};

#endif // GLDRAW_H
