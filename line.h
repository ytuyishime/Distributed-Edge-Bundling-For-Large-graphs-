#ifndef LINE_H
#define LINE_H

#include <assert.h>
#include <vector>
#include <math.h>
#include <stddef.h>
//#include <QtOpenGL>
//#include <QPointF>

using namespace std;

#define ATTRACTION          0.1

typedef float GLfloat;
typedef bool GLboolean;
typedef unsigned int GLuint;

class Line
{

public:

    Line(int sgm, GLfloat width);

    Line(GLfloat *p1, GLfloat *p2, GLboolean even, int cycle, int sgm,
         GLfloat width);

    GLfloat                 getWidth();

    GLfloat                 getOriLength();

    GLfloat                 getLength();

    GLfloat                 getConstant();

    GLfloat                 getSlope();

    GLfloat                 *getSource();

    GLfloat                 *getSink();

    GLfloat                 *getForce();

    GLfloat                 *getSegment();

    vector<int>             &getCompatibility();

    void                    setCompatibility(int i);

    void                    setSegment(GLfloat*);

    void                    setWidth(GLfloat);

    void                    setLength(GLfloat);

    void                    setConstant(GLfloat);

    void                    setSegmentNum(int);

    void                    setForce(GLfloat*);

protected:

private:

    void                    CalcSlope();

    void                    CalcSegments();

    void                    setEndpoints();

    void                    norVector(GLfloat *n);

    GLfloat                 CalcDistance(GLfloat *p1, GLfloat *p2);

    GLfloat                 *m_pSource;

    GLfloat                 *m_pSink;

    GLfloat                 *m_pSegment;

    GLfloat                 *m_pForce;

    GLfloat                 m_fConstant;

    GLfloat                 m_fSubLength;

    GLfloat                 m_fLength;

    GLfloat                 m_fOriLength;

    GLfloat                 m_fWidth;

    GLfloat                 m_fSlope;

    vector<int>             m_vCompatibility;

    GLuint                  m_nSegment;

};

#endif // LINE_H
