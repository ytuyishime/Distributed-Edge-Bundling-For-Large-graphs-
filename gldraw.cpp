#include "gldraw.h"

double                  simulation_time;
double                  compatibility_time;

static const float m_pQuadPoints[] = {
    -1.0f, -1.0f,
    1.0f, -1.0f,
    1.0f,  1.0f,
    -1.0f,  1.0f
};

static const float m_pTexPoints[] = {
    0.0f, 0.0f,
    1.0f, 0.0f,
    1.0f, 1.0f,
    0.0f, 1.0f
};

GLDraw::GLDraw()
{    
    m_bVertices     = true;
    m_bEdges        = true;
    m_bPolygons     = false;
    m_bDivEven      = false;
    m_bSmooth       = true;
    m_bCurveShaders = true;
    m_bNodeShaders  = !m_bCurveShaders;
    m_bLineWidthSet = true;
    m_bSegmentSet   = true;
    m_bData1        = true;
    m_nSegment      = 1;
    m_nLineNum      = 0;

    if (m_bDivEven) {
        m_nMax = pow(2, CYCLE) + 2;
    } else {
        m_nMax = pow(2, (CYCLE + 1)) - 1 + 2;
    }

//    m_sDataEdge = m_InputEdgeFile;
//    m_sDataNode = m_InputNodeFile;

//    m_sDataNode = "../CDATA/airlines_index_nodes.csv";
//    m_sDataEdge = "../CDATA/airlines_index_edges.csv";

//    m_sDataNode = "../CDATA/USMigrationFlowLarge_nodes.csv";
//    m_sDataEdge = "../CDATA/USMigrationFlowLarge_edges.csv";

//    m_sDataNode = "../CDATA/migrations_index_nodes.csv";
//    m_sDataEdge = "../CDATA/migrations_index_edges.csv";

//    m_sDataNode = "../CDATA/airlines_test_nodes.csv";
//    m_sDataEdge = "../CDATA/airlines_test_edges.csv";

//    m_sDataNode = "../CDATA/three_nodes.csv";
//    m_sDataEdge = "../CDATA/three_edges.csv";

//    m_sDataNode = "../CDATA/three_index_nodes.csv";
//    m_sDataEdge = "../CDATA/three_index_edges.csv";

//    m_sDataNode = "../CDATA/one_node.csv";
//    m_sDataEdge = "../CDATA/one_edge.csv";

//    m_sDataNode = "../CDATA/two_nodes.csv";
//    m_sDataEdge = "../CDATA/two_edges.csv";

//    m_sDataNode = "../CDATA/test_nodes.csv";
//    m_sDataEdge = "../CDATA/test_edges.csv";

//    m_sDataNode = "../CDATA/one_test_node.csv";
//    m_sDataEdge = "../CDATA/one_test_edge.csv";

}


double GLDraw::read_timer()
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

GLfloat GLDraw::calcDistance(GLfloat *p1, GLfloat *p2)
{
    GLfloat x = p1[0] - p2[0];
    GLfloat y = p1[1] - p2[1];
    GLfloat z = p1[2] - p2[2];
    return sqrt(x * x + y * y + z * z);
}

GLfloat GLDraw::getRandomArbitrary(GLfloat min, GLfloat max)
{
    return ((double) rand() / (RAND_MAX)) * (max - min) + min;
}

GLfloat GLDraw::vectorDotProduct(GLfloat *p1, GLfloat *p2)
{
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

void GLDraw::edgeAsVector(Line *l1, GLfloat *result)
{
    result[0] = l1->getSource()[0] - l1->getSink()[0];
    result[1] = l1->getSource()[1] - l1->getSink()[1];
    result[2] = l1->getSource()[2] - l1->getSink()[2];
}

void GLDraw::norVector(GLfloat *n)
{
    assert(n != NULL);
    GLfloat length = sqrt((n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]));
    n[0] = n[0] / length;
    n[1] = n[1] / length;
    n[2] = n[2] / length;
}

GLfloat GLDraw::angleCompatibility(Line *l1, Line *l2)
{
    GLfloat *tmp1 = new GLfloat[3];
    GLfloat *tmp2 = new GLfloat[3];
    memset(tmp1, 0.0, 3 * sizeof(GLfloat));
    memset(tmp2, 0.0, 3 * sizeof(GLfloat));
    edgeAsVector(l1, tmp1);
    edgeAsVector(l2, tmp2);
    GLfloat result = abs(vectorDotProduct(tmp1, tmp2) /
                        (calcDistance(l1->getSource(), l1->getSink()) *
                         calcDistance(l2->getSource(), l2->getSink())));
    delete[] tmp1;
    delete[] tmp2;

    return result;
}

GLfloat GLDraw::scaleCompatibility(Line *l1, Line *l2)
{
    GLfloat length1 = calcDistance(l1->getSource(), l1->getSink());
    GLfloat length2 = calcDistance(l2->getSource(), l2->getSink());
    GLfloat lavg = (length1 + length2) / 2.0;
    GLfloat result = 2.0 / (lavg / min(length1, length2) +
                                   max(length1, length2) / lavg);
    return result;
}

GLfloat GLDraw::positionCompatibility(Line *l1, Line *l2)
{
    GLfloat *midP = new GLfloat[3];
    GLfloat *midQ = new GLfloat[3];

    GLfloat length1 = calcDistance(l1->getSource(), l1->getSink());
    GLfloat length2 = calcDistance(l2->getSource(), l2->getSink());
    GLfloat lavg = (length1 + length2) / 2.0;

    midP[0] = (l1->getSource()[0] + l1->getSource()[0]) / 2.0;
    midP[1] = (l1->getSource()[1] + l1->getSource()[1]) / 2.0;
    midP[2] = (l1->getSource()[2] + l1->getSource()[2]) / 2.0;

    midQ[0] = (l2->getSource()[0] + l2->getSource()[0]) / 2.0;
    midQ[1] = (l2->getSource()[1] + l2->getSource()[1]) / 2.0;
    midQ[2] = (l2->getSource()[2] + l2->getSource()[2]) / 2.0;

    GLfloat result = lavg / (lavg + calcDistance(midP, midQ));

    delete[] midP;
    delete[] midQ;

    return result;
}

GLfloat GLDraw::visibilityCompatibility(Line *l1, Line *l2)
{
    return min(edgeVisibility(l1, l2), edgeVisibility(l2, l1));
}

GLfloat GLDraw::edgeVisibility(Line *l1, Line *l2)
{
    GLfloat *I0 = new GLfloat[3];
    GLfloat *I1 = new GLfloat[3];

    projectPointOnLine(l2->getSource(), l1->getSource(), l1->getSink(), I0);
    projectPointOnLine(l2->getSink(), l1->getSource(), l1->getSink(), I1);

    GLfloat *midI = new GLfloat[3];
    GLfloat *midP = new GLfloat[3];

    midI[0] = (I0[0] + I1[0]) / 2.0;
    midI[1] = (I0[1] + I1[1]) / 2.0;
    midI[2] = (I0[2] + I1[2]) / 2.0;

    midP[0] = ((l1->getSource()[0] + l1->getSink()[0]) / 2.0);
    midP[1] = ((l1->getSource()[1] + l1->getSink()[1]) / 2.0);
    midP[2] = ((l1->getSource()[2] + l1->getSink()[2]) / 2.0);

    GLfloat value = 0.0;
    GLfloat result = max(value, 1 - 2 * calcDistance(midP, midI) /
                                        calcDistance(I0, I1));
    delete[] I0;
    delete[] I1;
    delete[] midI;
    delete[] midP;

    return result;
}

// 2D case
void GLDraw::projectPointOnLine(GLfloat *P, GLfloat *Qsource,
                                GLfloat *Qsink, GLfloat *I)
{
    GLfloat x = Qsink[0] - Qsource[0];
    GLfloat y = Qsink[1] - Qsource[1];
    GLfloat z = Qsink[2] - Qsource[2];

    GLfloat L = sqrt(x * x + y * y + z * z);
    GLfloat r = ((Qsource[1] - P[1]) * (Qsource[1] - Qsink[1]) -
            (Qsource[0] - P[0]) * (Qsink[0] - Qsource[0])) / (L * L);
    I[0] = Qsource[0] + r * (Qsink[0] - Qsource[0]);
    I[1] = Qsource[1] + r * (Qsink[1] - Qsource[1]);
    I[2] = Qsource[2] + r * (Qsink[2] - Qsource[2]);
}

GLfloat GLDraw::compatibilityScore(Line *l1, Line *l2)
{
    float result = (angleCompatibility(l1, l2) *
                    scaleCompatibility(l1, l2) *
                    positionCompatibility(l1, l2) *
                    visibilityCompatibility(l1, l2));
    return result;
}

bool GLDraw::areCompatible(Line *l1, Line *l2)
{
    return (compatibilityScore(l1, l2) >= COM);
}

void GLDraw::computeCompatibility()
{
    float result;
    for (int i = 0; i < m_vLines.size(); i++) {
        for (int j = i + 1; j < m_vLines.size(); j++) {
            if (i == j) {
                continue;
            } else {
                //if (areCompatible(m_vLines[i], m_vLines[j])) {
                result = compatibilityScore(m_vLines[i], m_vLines[j]);
                if (result >= COM) {
                    m_vLines[i]->setCompatibility(j);
                    m_vLines[j]->setCompatibility(i);
                    m_pCompatibility[m_nLineNum * i + j] = result;
                    m_pCompatibility[m_nLineNum * j + i] = result;
                }
            }
        }
    }
}

/*
 * calculate the attraction force
 */
void GLDraw::calAttraction(GLfloat *A, GLfloat *B, GLfloat *Force, GLfloat kp)
{
    // Hooke's Law: F = -kx
    GLfloat force = kp;
    GLfloat v_force[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
    //norVector(v_force);

    assert(Force != NULL);
    if (abs(B[0] - A[0]) < EPS || abs(B[1] - A[1]) < EPS) {
        Force[0] = 0.0;
        Force[1] = 0.0;
        Force[2] = 0.0;
    } else {
        Force[0] = force * v_force[0];
        Force[1] = force * v_force[1];
        Force[2] = force * v_force[2];
    }
}

/*
 *calculate the electrostatic force
 */
void GLDraw::calElectrostatic(GLfloat *A, GLfloat *B, GLfloat *Force,
                              GLfloat particleLength) {
    // Coulomb's Law: F = k(Qq/r^2)
    GLfloat force = (ELECTRO / pow(particleLength, 1));
    GLfloat v_force[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
    //norVector(v_force);

    assert(Force != NULL);
    if (abs(B[0] - A[0]) < EPS || abs(B[1] - A[1]) < EPS) {
        Force[0] = 0.0;
        Force[1] = 0.0;
        Force[2] = 0.0;
    } else {
        Force[0] = force * v_force[0];
        Force[1] = force * v_force[1];
        Force[2] = force * v_force[2];
    }
}

void GLDraw::addForce(GLfloat *p1, GLfloat *p2, GLfloat *Force) {
    GLfloat xForce = p1[0] + p2[0];
    GLfloat yForce = p1[1] + p2[1];
    GLfloat zForce = p1[2] + p2[2];
    Force[0] = xForce;
    Force[1] = yForce;
    Force[2] = zForce;
}

void GLDraw::edge_midpoint(GLfloat *p1, GLfloat *p2, GLfloat *middle) {
    GLfloat vtr[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
    GLfloat length = calcDistance(p1, p2);

    norVector(vtr);

    vtr[0] = vtr[0] * (length / 2);
    vtr[1] = vtr[1] * (length / 2);
    vtr[2] = vtr[2] * (length / 2);

    middle[0] = p1[0] + vtr[0];
    middle[1] = p1[1] + vtr[1];
    middle[2] = p1[2] + vtr[2];
}

float GLDraw::safeAcos(float x)
{
    if (x < -1.0) {
        x = -1.0;
    } else if (x > 1.0) {
        x = 1.0;
    }
    return acos (x);
}

float GLDraw::length(float vec_x , float vec_y)
{
    return sqrt(vec_x * vec_x + vec_y * vec_y);
}

/*
 *Edge Division Calculation Methods
 */
void GLDraw::updateEdgeDivisions(int P) {
    GLfloat *midPoint = new GLfloat[3];

    for (int j = 0; j < m_vLines.size(); j++) {
        if (P == 1) {
            // source
            m_pSwapIn[m_nMax * j * 3 + 0] = m_vLines[j]->getSource()[0];
            m_pSwapIn[m_nMax * j * 3 + 1] = m_vLines[j]->getSource()[1];
            m_pSwapIn[m_nMax * j * 3 + 2] = m_vLines[j]->getSource()[2];

            // mid point
            edge_midpoint(m_vLines[j]->getSource(),
                          m_vLines[j]->getSink(),
                          midPoint);
            m_pSwapIn[m_nMax * j * 3 + 1 * 3 + 0] = midPoint[0];
            m_pSwapIn[m_nMax * j * 3 + 1 * 3 + 1] = midPoint[1];
            m_pSwapIn[m_nMax * j * 3 + 1 * 3 + 2] = midPoint[2];

            // target
            m_pSwapIn[m_nMax * j * 3 + 2 * 3 + 0] = m_vLines[j]->getSink()[0];
            m_pSwapIn[m_nMax * j * 3 + 2 * 3 + 1] = m_vLines[j]->getSink()[1];
            m_pSwapIn[m_nMax * j * 3 + 2 * 3 + 2] = m_vLines[j]->getSink()[2];
        } else {
            GLfloat divided_edge_length = m_vLines[j]->getLength();
            GLfloat segment_length  	= divided_edge_length / (P + 1);
            GLfloat current_segment_length = segment_length;

            //source
            m_pSwapOut[m_nMax * j * 3  + 0] = m_vLines[j]->getSource()[0];
            m_pSwapOut[m_nMax * j * 3  + 1] = m_vLines[j]->getSource()[1];
            m_pSwapOut[m_nMax * j * 3  + 2] = m_vLines[j]->getSource()[2];

            //subdivision points
            int size;
            if (m_bDivEven) {
                size = (P / 2) + 2;
            } else {
                size = ((P - 1) / 2) + 2;
            }
            size_t size0 = 1;
            for (int i = 1; i < size; i++) {
                if (m_bDivEven) {
                    GLfloat old_segment_length
                    = calcDistance(&m_pSwapIn[m_nMax * j * 3 + i * 3 + 0],
                             &m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 0]);
                    while (old_segment_length - current_segment_length > 1e-6) {
                        GLfloat percent_position = current_segment_length
                                / old_segment_length;
                        GLfloat new_subdivision_point_x
                             = m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 0];
                        GLfloat new_subdivision_point_y
                             = m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 1];
                        GLfloat new_subdivision_point_z
                             = m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 2];

                        new_subdivision_point_x += percent_position
                            * (m_pSwapIn[m_nMax * j * 3 + i * 3 + 0]
                            - m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 0]);
                        new_subdivision_point_y += percent_position
                            * (m_pSwapIn[m_nMax * j * 3 + i * 3 + 1]
                            - m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 1]);
                        new_subdivision_point_z += percent_position
                            * (m_pSwapIn[m_nMax * j * 3 + i * 3 + 2]
                            - m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 2]);

                        m_pSwapOut[m_nMax * j * 3 + size0 * 3  + 0]
                                = new_subdivision_point_x;
                        m_pSwapOut[m_nMax * j * 3 + size0 * 3  + 1]
                                = new_subdivision_point_y;
                        m_pSwapOut[m_nMax * j * 3 + size0 * 3  + 2]
                                = new_subdivision_point_z;

                        old_segment_length    -= current_segment_length;
                        current_segment_length = segment_length;
                        size0++;
                    }
                    current_segment_length -= old_segment_length;
                } else {
                    GLfloat new_subdivision_point_x
                            = (m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 0]
                            + m_pSwapIn[m_nMax * j * 3 + i * 3 + 0]) / 2;
                    GLfloat new_subdivision_point_y
                            = (m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 1]
                            + m_pSwapIn[m_nMax * j * 3 + i * 3 + 1]) / 2;
                    GLfloat new_subdivision_point_z
                            = (m_pSwapIn[m_nMax * j * 3 + (i - 1) * 3 + 2]
                            + m_pSwapIn[m_nMax * j * 3 + i * 3 + 2]) / 2;

                    m_pSwapOut[m_nMax * j * 3 + size0 * 3 + 0]
                            = new_subdivision_point_x;
                    m_pSwapOut[m_nMax * j * 3 + size0 * 3 + 1]
                            = new_subdivision_point_y;
                    m_pSwapOut[m_nMax * j * 3 + size0 * 3 + 2]
                            = new_subdivision_point_z;
                    size0++;

                    m_pSwapOut[m_nMax * j * 3 + size0 * 3 + 0]
                            = m_pSwapIn[m_nMax * j * 3 + i * 3 + 0];
                    m_pSwapOut[m_nMax * j * 3 + size0 * 3 + 1]
                            = m_pSwapIn[m_nMax * j * 3 + i * 3 + 1];
                    m_pSwapOut[m_nMax * j * 3 + size0 * 3 + 2]
                            = m_pSwapIn[m_nMax * j * 3 + i * 3 + 2];
                    size0++;
                }
            }

            //target
            m_pSwapOut[m_nMax * j * 3  + (P + 1) * 3 + 0]
                    = m_vLines[j]->getSink()[0];
            m_pSwapOut[m_nMax * j * 3  + (P + 1) * 3 + 1]
                    = m_vLines[j]->getSink()[1];
            m_pSwapOut[m_nMax * j * 3  + (P + 1) * 3 + 2]
                    = m_vLines[j]->getSink()[2];
        }
    }

    delete[] midPoint;
    midPoint = NULL;
}

/*
 * update the line length
 */
void GLDraw::updateLength(int P)
{
    GLfloat length = 0;
    for (int i = 0; i < m_vLines.size(); i++) {
        if (abs(m_vLines[i]->getSource()[0] - m_vLines[i]->getSink()[0]) < EPS
         && abs(m_vLines[i]->getSource()[1] - m_vLines[i]->getSink()[1]) < EPS)
        {
            m_vLines[i]->setLength(EPS);
            m_vLines[i]->setConstant(ATTRACTION);
        } else {
            length = 0;
            for (int j = 1; j <= P + 1; j++) {
                GLfloat segment_length;
                segment_length
                        = calcDistance(&m_pSwapOut[m_nMax * i * 3 + j * 3 + 0],
                        &m_pSwapOut[m_nMax * i * 3 + (j - 1) * 3 + 0]);
                length += segment_length;
            }
            m_vLines[i]->setLength(length);
            m_vLines[i]->setConstant(ATTRACTION / (m_vLines[i]->getOriLength()
                                                   * (P + 1)));
            //m_vLines[i]->setConstant(ATTRACTION / m_vLines[i]->getLength());
        }
    }
}

/*
 * update the forces and positions of the subdivision points
 */
void GLDraw::updateForcePosition(int P, float S)
{

    GLfloat *attr_left_force = new GLfloat[3];
    GLfloat *attr_right_force = new GLfloat[3];
    GLfloat *attr_force   = new GLfloat[3];

    GLfloat *electro_force_acc = new GLfloat[3];
    GLfloat *electro_force  = new GLfloat[3];
    GLfloat *total = new GLfloat[3];

    memset(electro_force, 0.0, 3 * sizeof(GLfloat));
    memset(electro_force_acc, 0.0, 3 * sizeof(GLfloat));
    memset(attr_left_force, 0.0, 3 * sizeof(GLfloat));
    memset(attr_right_force, 0.0, 3 * sizeof(GLfloat));
    memset(attr_force, 0.0, 3 * sizeof(GLfloat));
    memset(total, 0.0, 3 * sizeof(GLfloat));

    for (int i = 0; i < m_vLines.size(); i++) {
        // source
        m_pSwapOut[m_nMax * i * 3 + 0] = m_pSwapIn[m_nMax * i * 3 + 0];
        m_pSwapOut[m_nMax * i * 3 + 1] = m_pSwapIn[m_nMax * i * 3 + 1];
        m_pSwapOut[m_nMax * i * 3 + 2] = m_pSwapIn[m_nMax * i * 3 + 2];
        for (int j = 1; j < P + 1; j++) {
            calAttraction(&m_pSwapIn[m_nMax * i * 3 + j * 3 + 0],
                          &m_pSwapIn[m_nMax * i * 3 + (j - 1) * 3 + 0],
                          attr_left_force, m_vLines[i]->getConstant());

            calAttraction(&m_pSwapIn[m_nMax * i * 3 + j * 3 + 0],
                          &m_pSwapIn[m_nMax * i * 3 + (j + 1) * 3 + 0],
                          attr_right_force, m_vLines[i]->getConstant());

            addForce(attr_left_force, attr_right_force, attr_force);

            memset(electro_force, 0.0, 3 * sizeof(GLfloat));
//            for (int k = 0; k < m_vLines.size(); k++) {
            for (int k = 0; k < m_vLines[i]->getCompatibility().size(); k++) {
                int target = m_vLines[i]->getCompatibility()[k];
//                int target = k;
//                if (target != i &&
                if (m_pCompatibility[m_nLineNum * i + target] >= COM) {
                    calElectrostatic(
                        &m_pSwapIn[m_nMax * i * 3 + j * 3 + 0],
                        &m_pSwapIn[m_nMax * target * 3 + j * 3 + 0],
                        electro_force_acc, calcDistance(
                            &m_pSwapIn[m_nMax * i * 3 + j * 3 + 0],
                            &m_pSwapIn[m_nMax * target * 3 + j * 3 + 0]));
                }
                electro_force[0] += electro_force_acc[0];
                electro_force[1] += electro_force_acc[1];
                electro_force[2] += electro_force_acc[2];
            }
            addForce(attr_force, electro_force, total);

            // subdivision
            m_pSwapOut[m_nMax * i * 3 + j * 3 + 0]
                    = m_pSwapIn[m_nMax * i * 3 + j * 3 + 0] + S * total[0];
            m_pSwapOut[m_nMax * i * 3 + j * 3 + 1]
                    = m_pSwapIn[m_nMax * i * 3 + j * 3 + 1] + S * total[1];
            m_pSwapOut[m_nMax * i * 3 + j * 3 + 2]
                    = m_pSwapIn[m_nMax * i * 3 + j * 3 + 2] + S * total[2];

        }

        // target
        m_pSwapOut[m_nMax * i * 3 + (P + 1) * 3 + 0]
                = m_pSwapIn[m_nMax * i * 3 + (P + 1) * 3 + 0];
        m_pSwapOut[m_nMax * i * 3 + (P + 1) * 3 + 1]
                = m_pSwapIn[m_nMax * i * 3 + (P + 1) * 3 + 1];
        m_pSwapOut[m_nMax * i * 3 + (P + 1) * 3 + 2]
                = m_pSwapIn[m_nMax * i * 3 + (P + 1) * 3 + 2];

        memset(electro_force, 0.0, 3 * sizeof(GLfloat));
        memset(electro_force_acc, 0.0, 3 * sizeof(GLfloat));
        memset(attr_left_force, 0.0, 3 * sizeof(GLfloat));
        memset(attr_right_force, 0.0, 3 * sizeof(GLfloat));
        memset(attr_force, 0.0, 3 * sizeof(GLfloat));
        memset(total, 0.0, 3 * sizeof(GLfloat));
    }

    delete[] attr_left_force;
    attr_left_force = NULL;

    delete[] attr_right_force;
    attr_right_force = NULL;

    delete[] attr_force;
    attr_force = NULL;

    delete[] electro_force;
    electro_force = NULL;

    delete[] electro_force_acc;
    electro_force_acc = NULL;

    delete[] total;
    total = NULL;

    GLfloat length = 0;
    for (int i = 0; i < m_vLines.size(); i++) {
        if (abs(m_vLines[i]->getSource()[0] - m_vLines[i]->getSink()[0]) < EPS
         && abs(m_vLines[i]->getSource()[1] - m_vLines[i]->getSink()[1]) < EPS)
        {
            m_vLines[i]->setLength(EPS);
            m_vLines[i]->setConstant(ATTRACTION);
        } else {
            length = 0;
            for (int j = 1; j <= P + 1; j++) {
                GLfloat segment_length;
                segment_length
                        = calcDistance(&m_pSwapOut[m_nMax * i * 3 + j * 3 + 0],
                        &m_pSwapOut[m_nMax * i * 3 + (j - 1) * 3 + 0]);
                length += segment_length;
            }
            m_vLines[i]->setLength(length);
            m_vLines[i]->setConstant(ATTRACTION / (m_vLines[i]->getOriLength()
                                                   * (P + 1)));
            //m_vLines[i]->setConstant(ATTRACTION / m_vLines[i]->getLength());
        }
    }
}

/*
 * edge bundling main function
 */
void GLDraw::updateBuffer()
{
    float   S     = S_INITIAL;    //distance to move points
    int     I     = I_INITIAL;    //number of iterations for cycle
    int     P     = P_INITIAL;    //subdivision number

    initPosBuffer();
    updateEdgeDivisions(P);

    compatibility_time = read_timer();

    computeCompatibility();

    compatibility_time = read_timer() - compatibility_time;

    simulation_time = read_timer();

    double iteration_time = 0.0;
    for (int cycle = 0; cycle < CYCLE; cycle++) {
        for (int i = 0; i < I; i++) {
            updateForcePosition(P, S);
            swapPosBuffer();
        }
        //prepare for next cycle
        S = S / 2;
        if (m_bDivEven) {
            P = P * P_RATE;
        } else {
            P = P * P_RATE + 1;
        }
        I = I * I_RATE;
        m_nSegment = P;

        updateEdgeDivisions(P);
        swapPosBuffer();
    }

    simulation_time = read_timer() - simulation_time;

    if (m_bSmooth && CYCLE != 0) {
        gaussianSmooth(P);
    }
    m_nSegment++;

//    calcCompatibilityDistance();
}

void GLDraw::initCompatibilityDistance()
{
    GLfloat f_dis = 0.0;
    for (int i = 0; i < m_vLines.size(); i++) {
        for (int k = 0; k < m_vLines[i]->getCompatibility().size(); k++) {
            int target = m_vLines[i]->getCompatibility()[k];
            if (target != i) {
                for (int j = 0; j < m_nMax; j++) {
                    f_dis += calcDistance(&m_pOriginPoints[m_nMax * i * 3 +
                                                           j * 3 + 0],
                                          &m_pOriginPoints[m_nMax * target * 3 +
                                                           j * 3 + 0]);
                }
            }
        }
    }
    cout << "before bundled distance: " << f_dis << endl;
}

void GLDraw::calcCompatibilityDistance()
{
    GLfloat f_dis = 0.0;
    for (int i = 0; i < m_vLines.size(); i++) {
        for (int k = 0; k < m_vLines[i]->getCompatibility().size(); k++) {
            int target = m_vLines[i]->getCompatibility()[k];
            if (target != i) {
                for (int j = 0; j < m_nMax; j++) {
                    f_dis += calcDistance(&m_pSmoothPoints[m_nMax * i * 3 +
                                                           j * 3 + 0],
                                          &m_pSmoothPoints[m_nMax * target * 3 +
                                                           j * 3 + 0]);
                }
            }
        }
    }
    cout << "after bundled distance: " << f_dis << endl;
}

void GLDraw::saveBuffer()
{
    string file = "../DATA";
    ofstream output(file.c_str());

    float   attraction  = ATTRACTION;
    float   electro     = ELECTRO;

    float   eps         = EPS;

    int     segment     = SEGMENT;

    float   s_initial   = S_INITIAL;
    float   i_rate      = I_RATE;

    int     p_initial   = P_INITIAL;
    int     p_rate      = P_RATE;
    int     cycle       = CYCLE;
    int     i_initial   = I_INITIAL;
    int     size        = 9;
    int     lines       = m_nLineNum;
    int     maximum     = m_nMax;
    int     ballnum     = BALLNUM;

    assert(m_vLines.size() == m_nLineNum);
    GLfloat *subpoints  = new GLfloat[m_nLineNum * m_nMax * 3];

    memset(subpoints, 0, sizeof(GLfloat) * m_nLineNum * m_nMax * 3);

    for (int i = 0; i < m_nLineNum; i++) {
        for (int j = 0; j < size; j++) {
            subpoints[m_nMax * i * 3 + j] = m_vLines[i]->getSegment()[j];
        }
    }

    if (m_vLines.size() == 0) {
        cout << "There is no lines, no buffer is saved to the disk." << endl;
    }

    /* ATTRACTION */
    output.write(reinterpret_cast<char *>(&attraction), sizeof(float));

    /* REPULSION */
    output.write(reinterpret_cast<char *>(&electro), sizeof(float));

    /* EPS */
    output.write(reinterpret_cast<char *>(&eps), sizeof(float));

    /* SEGMENT */
    output.write(reinterpret_cast<char *>(&segment), sizeof(float));

    /* S_INITIAL */
    output.write(reinterpret_cast<char *>(&s_initial), sizeof(float));

    /* P_INITIAL */
    output.write(reinterpret_cast<char *>(&p_initial), sizeof(int));

    /* P_RATE */
    output.write(reinterpret_cast<char *>(&p_rate), sizeof(int));

    /* CYCLE */
    output.write(reinterpret_cast<char *>(&cycle), sizeof(int));

    /* I_INITIAL */
    output.write(reinterpret_cast<char *>(&i_initial), sizeof(int));

    /* I_RATE */
    output.write(reinterpret_cast<char *>(&i_rate), sizeof(float));

    /* LINE NUMBER */
    output.write(reinterpret_cast<char *>(&lines), sizeof(int));

    /* MAXIMUM NUMBER */
    output.write(reinterpret_cast<char *>(&maximum), sizeof(int));

    /* LINE SUBPOINTS */
    float subpoint;
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < maximum * 3; j++) {
            subpoint = subpoints[maximum * i * 3 + j];
            output.write(reinterpret_cast<char *>(&subpoint), sizeof(float));
        }
    }

    /* BALL NUMBER */
    //output.write(reinterpret_cast<char *>(&ballnum), sizeof(int));

    output.close();

    /*
    string check = "../check";
    ofstream check_1(check.c_str());
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < 3 * 3; j++) {
            subpoint = subpoints[maximum * i * 3 + j];
            check_1 << subpoint << endl;
        }
    }
    */
}

void GLDraw::readBuffer()
{
    string file = "../DATA";
    ifstream input(file.c_str());

    float   attraction;
    float   electro;
    float   eps;
    int     segment;

    float   s_initial;
    int     p_initial;
    int     p_rate;
    int     cycle;
    int     i_initial;
    float   i_rate;

    int     lines;
    int     maximum;

    //int     ballnum;

    /* ATTRACTION */
    input.read(reinterpret_cast<char *>(&attraction), sizeof(float));

    /* REPULSION */
    input.read(reinterpret_cast<char *>(&electro), sizeof(float));

    /* EPS */
    input.read(reinterpret_cast<char *>(&eps), sizeof(float));

    /* SEGMENT */
    input.read(reinterpret_cast<char *>(&segment), sizeof(float));

    /* S_INITIAL */
    input.read(reinterpret_cast<char *>(&s_initial), sizeof(float));

    /* P_INITIAL */
    input.read(reinterpret_cast<char *>(&p_initial), sizeof(int));

    /* P_RATE */
    input.read(reinterpret_cast<char *>(&p_rate), sizeof(int));

    /* CYCLE */
    input.read(reinterpret_cast<char *>(&cycle), sizeof(int));

    /* I_INITIAL */
    input.read(reinterpret_cast<char *>(&i_initial), sizeof(int));

    /* I_RATE */
    input.read(reinterpret_cast<char *>(&i_rate), sizeof(float));

    /* LINE NUMBER */
    input.read(reinterpret_cast<char *>(&lines), sizeof(int));

    /* MAXIMUM NUMBER */
    input.read(reinterpret_cast<char *>(&maximum), sizeof(int));

    /* BALL NUMBER */
    //input.read(reinterpret_cast<char *>(&ballnum), sizeof(int));

    /* LINE SUBPOINTS */
    float subpoint;
    float *subpoints = new float[m_nLineNum * m_nMax * 3];
    for (int i = 0; i < m_nLineNum; i++) {
        for (int j = 0; j < m_nMax * 3; j++) {
            input.read(reinterpret_cast<char *>(&subpoint), sizeof(float));
            subpoints[m_nMax * i * 3 + j] = subpoint;
        }
    }

    input.close();

//    qDebug() << attraction;

//    qDebug() << electro;

//    qDebug() << force;

//    qDebug() << distance;

//    qDebug() << s_initial;

//    qDebug() << i_rate;

//    qDebug() << p_initial;

//    qDebug() << p_rate;

//    qDebug() << c;

//    qDebug() << i_initial;

//    qDebug() << segment;

//    qDebug() << size;

//    qDebug() << lines;

//    for (int i = 0; i < n_LineNum; i++) {
//        for (int j = 0; j < size; j++) {
//            qDebug() << subpoints[i * size + j];
//        }
//    }
}

void GLDraw::initBuffer()
{
    //m_pGL_dat_points    = new GLfloat[18 * m_nBallNum];

    //m_pGL_dat_colors    = new GLfloat[24 * m_nBallNum];

    //m_pGL_dat_centers   = new GLfloat[18 * m_nBallNum];

    //m_pGL_dat_radius    = new GLfloat[6 * m_nBallNum];

    //m_pGL_dat_lines     = new GLfloat[6 * m_nSegment * m_vLines.size()];
    //m_pGL_dat_lines     = new GLfloat[3 * m_nMax * m_vLines.size()];
    //m_pGL_dat_lines     = new GLfloat[6 * SEGMENT * (m_vLines.size() + 2)];

    // debug for DFS program
    m_pGL_dat_lattri    = new GLfloat[3 * m_nMax * m_vLines.size()];
//    m_pGL_dat_lattri    = new GLfloat[3 * m_nMax * m_nLineNum + 9];

    m_pGL_dat_centers   = new GLfloat[3 * m_nBallNum];

    m_pGL_dat_colors    = new GLfloat[4 * m_nBallNum];

    m_pGL_dat_polygons  = new GLfloat[6 * (m_nSegment + 1) * m_vLines.size()];

    // debug for DFS program
//    m_pGL_dat_lines     = new GLfloat[3 * m_nMax * m_nLineNum + 9];

    int count = 0;

    size_t size0 = 0;
    size_t size1 = 0;
    size_t size2 = 0;
    size_t size3 = 0;
    for (int i = 0; i < m_nBallNum; i++) {
        for (int j = 0; j < 3; j++) {
            m_pGL_dat_centers[size0++] = m_vPoints[i]->getPosition()[j];
        }
        for (int k = 0; k < 4; k++) {
            m_pGL_dat_colors[size1++] = m_vPoints[i]->getColor()[k];
        }
//        if (m_pTexLimit[i] == 1) {
//            count++;
//            for (int j = 0; j < 18; j++) {
//                m_pGL_dat_points[size0++] = m_vPoints[i]->getVertices()[j];
//            }
//        }

//        for (int j = 0; j < 6; j++) {
//            for (int k = 0; k < 4; k++) {
//                m_pGL_dat_colors[size1++] = m_vPoints[i]->getColor()[k];
//            }
//            for (int l = 0; l < 3; l++) {
//                m_pGL_dat_centers[size2++] = m_vPoints[i]->getPosition()[l];
//            }
//            m_pGL_dat_radius[size3++] = m_vPoints[i]->getRadius();
//        }
    }

    size0 = 0;
    GLfloat *outpoints;
    if (m_bSmooth && CYCLE != 0) {
        outpoints = m_pSmoothPoints;
    } else {
        outpoints = m_pSwapIn;
    }

    // debug for DFS program
//    size_t size = 0;
//    for (int i = 0; i < m_nLineNum * m_nMax * 3; i++) {
//        m_pGL_dat_lines[size++] = outpoints[i];
//    }
//    m_pGL_dat_lines[size++] = -0.5;
//    m_pGL_dat_lines[size++] = 0.0;
//    m_pGL_dat_lines[size++] = 0.0;

//    m_pGL_dat_lines[size++] = 0.0;
//    m_pGL_dat_lines[size++] = -0.1;
//    m_pGL_dat_lines[size++] = 0.0;

//    m_pGL_dat_lines[size++] = 0.5;
//    m_pGL_dat_lines[size++] = 0.0;
//    m_pGL_dat_lines[size++] = 0.0;
    m_pGL_dat_lines = outpoints;

//    float total_angle = 0.0;
//    float indiv_angle = 0.0;

    float f_distortion = 0.0;

    assert(m_nSegment + 1 == m_nMax);
    for (int i = 0; i < m_vLines.size(); i++) {
//        float f_source_x = m_pGL_dat_lines[3 * i * m_nMax + 0];
//        float f_source_y = m_pGL_dat_lines[3 * i * m_nMax + 1];
//        float f_source_z = m_pGL_dat_lines[3 * i * m_nMax + 2];
//        float f_target_x = m_pGL_dat_lines[3 * i * m_nMax + 3 * m_nMax + 0];
//        float f_target_y = m_pGL_dat_lines[3 * i * m_nMax + 3 * m_nMax + 1];
//        float f_target_z = m_pGL_dat_lines[3 * i * m_nMax + 3 * m_nMax + 2];
//        float f_vector_x = f_target_x - f_source_x;
//        float f_vector_y = f_target_y - f_source_y;
//        float f_vector_z = f_target_z - f_source_z;
//        cout << m_pGL_dat_lines[3 * i * m_nMax + 3 * (m_nMax - 1) + 0] << endl;
//        cout << m_vLines[i]->getSegment()[3 * (m_nMax - 1) + 0] << endl;
//        cout << m_vLines[i]->getSink()[0] << endl;
        assert(m_pGL_dat_lines[3 * i * m_nMax + 0] - m_vLines[i]->getSegment()[0] == 0);
        assert(m_pGL_dat_lines[3 * i * m_nMax + 1] - m_vLines[i]->getSegment()[1] == 0);
        assert(m_pGL_dat_lines[3 * i * m_nMax + 2] - m_vLines[i]->getSegment()[2] == 0);
        assert(m_pGL_dat_lines[3 * i * m_nMax + 3 * (m_nMax - 1) + 0] - m_vLines[i]->getSegment()[3 * (m_nMax - 1) + 0] == 0);
        assert(m_pGL_dat_lines[3 * i * m_nMax + 3 * (m_nMax - 1) + 1] - m_vLines[i]->getSegment()[3 * (m_nMax - 1) + 1] == 0);
        assert(m_pGL_dat_lines[3 * i * m_nMax + 3 * (m_nMax - 1) + 2] - m_vLines[i]->getSegment()[3 * (m_nMax - 1) + 2] == 0);
        for (int j = 0; j < m_nSegment + 1; j++) {
            m_pGL_dat_lattri[size0++] = (m_vLines[i]->getSlope() *
                                         (180.0 / PI)) / 180.0;
            m_pGL_dat_lattri[size0++] = m_vLines[i]->getLength();
            m_pGL_dat_lattri[size0++] = (j * 1.0f) / (m_nSegment) * 1.0f;
            float f_source[3];
            float f_target[3];
            f_source[0] = m_pGL_dat_lines[3 * i * m_nMax + 0];
            f_source[1] = m_pGL_dat_lines[3 * i * m_nMax + 1];
            f_source[2] = m_pGL_dat_lines[3 * i * m_nMax + 2];
            f_target[0] = m_pGL_dat_lines[3 * i * m_nMax + 3 * m_nMax + 0];
            f_target[1] = m_pGL_dat_lines[3 * i * m_nMax + 3 * m_nMax + 1];
            f_target[2] = m_pGL_dat_lines[3 * i * m_nMax + 3 * m_nMax + 2];
            f_distortion += calcDistance(f_source + 0, f_target + 0);
//            m_pGL_dat_lines[3 * i * m_nMax + 3 * j + 0] = m_vLines[i]->getSegment()[3 * j + 0];
//            m_pGL_dat_lines[3 * i * m_nMax + 3 * j + 1] = m_vLines[i]->getSegment()[3 * j + 1];
//            m_pGL_dat_lines[3 * i * m_nMax + 3 * j + 2] = m_vLines[i]->getSegment()[3 * j + 2];
        }
    }

    cout << "total distortion: " << f_distortion / (m_nMax * m_nLineNum * 1.0) << endl;

//    cout << "total maximum angle: " << total_angle << endl;

    // debug for DFS program
//    m_pGL_dat_lattri[size0++] = 1.0;
//    m_pGL_dat_lattri[size0++] = 1.0;
//    m_pGL_dat_lattri[size0++] = 1.0;

//    m_pGL_dat_lattri[size0++] = 1.0;
//    m_pGL_dat_lattri[size0++] = 1.0;
//    m_pGL_dat_lattri[size0++] = 1.0;

//    m_pGL_dat_lattri[size0++] = 1.0;
//    m_pGL_dat_lattri[size0++] = 1.0;
//    m_pGL_dat_lattri[size0++] = 1.0;

    cout << "total line ponits: " << size0 << endl;
    cout << "total line ponits: " << m_nMax * m_vLines.size() << endl;

    /*
    string file = "../groundtruth_output";
    ofstream output(file.c_str());
    for (int i = 0; i < m_nLineNum * m_nMax * 3; i++) {
        output << outpoints[i] << endl;
    }
    output.close();
    */

    m_pAdjMatrix_output = new GLfloat[m_nBallNum * m_nBallNum];
    memset(m_pAdjMatrix_output, 0.0, sizeof(GLfloat) * m_nBallNum * m_nBallNum);

    int n_count = 0;

    GLint id, source, target;
    ifstream infile(m_sDataEdge.c_str());
    if (infile.is_open()) {
        while (!infile.eof()) {
//            if (count >= 2180) {
//                break;
//            }
            infile >> id >> source >> target;
            m_pAdjMatrix_output[m_nBallNum * source + target] = 1.0;
            m_pAdjMatrix_output[m_nBallNum * target + source] = 1.0;

            n_count++;
            if (infile.fail()) {
                break;
            }
        }
    }
    infile.close();

    ofstream outf("../comMatrix_3");
    for (int i = 0; i < m_nBallNum; i++) {
        for (int j = 0; j < m_nBallNum; j++) {
            assert(m_pAdjMatrix_output[i * m_nBallNum + j] ==
                   m_pAdjMatrix_output[j * m_nBallNum + i]);
            outf << m_pAdjMatrix_output[i * m_nBallNum + j];
            if (j != m_nBallNum - 1) {
                outf << ",";
            }
        }
        outf << endl;
    }
    outf.close();


    cout << "Number of effective vertices is "
         << count
         << endl;

    cout << "The total number of element is "
         << m_vLines.size() * m_nMax * 3
         << endl;

    cout << "The compatibility time is "
         << compatibility_time
         << " seconds"
         << endl;

    cout << "The simulation time is "
         << simulation_time
         << " seconds"
         << endl;

    cout << "The total time is "
         << compatibility_time + simulation_time
         << " seconds"
         << endl;


    /*
    outf.open("../dat_lattri");
    if (!outf) {
        cout << "Cannot open file data_lattri" << endl;
    } else {
        cout << "Write data_lattri" << endl;
    }
    outf.write(reinterpret_cast<char *>(m_pGL_dat_lattri), sizeof(GLfloat) * 3 * m_nMax * m_vLines.size());
    outf.close();


    outf.open("../dat_centers");
    if (!outf) {
        cout << "Cannot open file dat_centers" << endl;
    } else {
        cout << "Write dat_centers" << endl;
    }
    outf.write(reinterpret_cast<char *>(m_pGL_dat_centers), sizeof(GLfloat) * 3 * m_nBallNum);
    outf.close();


    outf.open("../dat_colors");
    if (!outf) {
        cout << "Cannot open file dat_colors" << endl;
    } else {
        cout << "Write dat_colors" << endl;
    }
    outf.write(reinterpret_cast<char *>(m_pGL_dat_colors), sizeof(GLfloat) * 4 * m_nBallNum);
    outf.close();

    outf.open("../dat_polygons");
    if (!outf) {
        cout << "Cannot open file dat_polygons" << endl;
    } else {
        cout << "Write dat_polygson" << endl;
    }
    outf.write(reinterpret_cast<char *>(m_pGL_dat_polygons), sizeof(GLfloat) * 6 * (m_nSegment + 1) * m_vLines.size());
    outf.close();


    outf.open("../dat_matrix");
    if (!outf) {
        cout << "Cannot open file dat_matrix" << endl;
    } else {
        cout << "Write dat_matrix" << endl;
    }
    outf.write(reinterpret_cast<char *>(m_pAdjMatrix_output), sizeof(GLfloat) * m_nBallNum * m_nBallNum);
    outf.close();


    outf.open("../data_smoothpoints");
    if (!outf) {
        cout << "Cannot open file data_smoothpoints" << endl;
    } else {
        cout << "Write data_smoothpoints" << endl;
    }
    outf.write(reinterpret_cast<char *>(m_pSmoothPoints), sizeof(GLfloat) * m_nLineNum * m_nMax * 3);
    outf.close();
    */


}

/*
 * initialize the vertex
 */
void GLDraw::initVertex()
{
    GLfloat id, map_x, map_y;

    cout << "Open " << m_sDataNode << endl;
     //cout<<"HERE "<<m_sDataNode<< " lol";
    ifstream infile(m_sDataNode.c_str());
    if (infile.is_open()) {
        while (!infile.eof()) {
            //infile >> id >> x >> y >> map_x >> map_y;
            infile >> id >> map_x >> map_y;
	    //cout<<id<<" " <<map_x<<" "<< map_y<<endl;
            if (infile.fail()) {
                break;
            }
            //cout << id << ", " << map_x << ", " << map_y << endl;
            m_pPoint = new Point(map_x * 0.9, map_y * 0.9, RADIUS);
            //m_pPoint = new Point(map_x, map_y, RADIUS);
            m_vPoints.push_back(m_pPoint);
        }
    } else {
        cout << "Cannot open " << m_sDataNode << endl;

    }
    //cout<<"closed"<<endl;
    infile.close();

    m_nBallNum = m_vPoints.size();

    cout << m_nBallNum << endl;
}

/*
 *initialize the distance matrix
 */
void GLDraw::initAdjMatrix()
{
    m_pAdjMatrix = new GLfloat[m_nBallNum * m_nBallNum];
    memset(m_pAdjMatrix, 0.0, sizeof(GLfloat) * m_nBallNum * m_nBallNum);

    m_pTexLimit = new GLint[m_nBallNum];
    memset(m_pTexLimit, 0.0, sizeof(GLint) * m_nBallNum);

    int count = 0;

    GLint id, source, target;

    ifstream infile(m_sDataEdge.c_str());

    cout << "Open " << m_sDataEdge << endl;

    if (infile.is_open()) {
        while (!infile.eof()) {
//            if (count >= 2180) {
		  //cout<<"count is greater than 2180"<<endl;
//                break;
//            }
            infile >> id >> source >> target;
	  //  cout<<id<<" " <<source<< " " <<target<<endl;
            m_pAdjMatrix[m_nBallNum * source + target] = 1.0;
            //m_pAdjMatrix[m_nBallNum * target + source] = 1.0;

            m_pTexLimit[source] = 1;
            m_pTexLimit[target] = 1;

            count++;
            if (infile.fail()) {
                break;
            }
        }
    } else {
        cout << "Cannot open " << m_sDataEdge << endl;
    }


    infile.close();
}

/*
 * initialize the lines and distance
 */
void GLDraw::initDistance()
{
    m_fMaxLength = 0;
    assert(m_pAdjMatrix != NULL);
    for (int i = 0; i < m_nBallNum; i++) {
        for (int j = 0; j < m_nBallNum; j++) {
            if (m_pAdjMatrix[m_nBallNum * i + j] == 1.0) {
                GLfloat distance = calcDistance(m_vPoints[i]->getPosition(),
                                                m_vPoints[j]->getPosition());
                if (distance != 0.0) {
                    m_pLine = new Line(m_vPoints[i]->getPosition(),
                                       m_vPoints[j]->getPosition(),
                                       //m_bDivEven, CYCLE, P_INITIAL, LINEWIDTH);
                                       m_bDivEven, CYCLE, m_nMax, LINEWIDTH);
                    m_vLines.push_back(m_pLine);
                    m_fMaxLength = m_fMaxLength > m_pLine->getLength() ?
                                   m_fMaxLength : m_pLine->getLength();
                    m_nLineNum++;
                }
            }
        }
    }

    m_pSmoothPoints    = new GLfloat[m_nLineNum * m_nMax * 3];
    memset(m_pSmoothPoints, 0, sizeof(GLfloat) * m_nLineNum * m_nMax * 3);

    m_pOriginPoints    = new GLfloat[m_nLineNum * m_nMax * 3];
    memset(m_pOriginPoints, 0, sizeof(GLfloat) * m_nLineNum * m_nMax * 3);

    cout << "Total number of lines: " << m_nLineNum << endl;
}


void GLDraw::initPosBuffer()
{
    m_pSubPoints_1     = new GLfloat[m_nLineNum * m_nMax * 3];
    m_pSubPoints_2     = new GLfloat[m_nLineNum * m_nMax * 3];
    m_pCompatibility   = new GLfloat[m_nLineNum * m_nLineNum];

    memset(m_pSubPoints_1, 0, sizeof(GLfloat) * m_nLineNum * m_nMax * 3);
    memset(m_pSubPoints_2, 0, sizeof(GLfloat) * m_nLineNum * m_nMax * 3);
    memset(m_pCompatibility, 0, sizeof(GLfloat) * m_nLineNum * m_nLineNum);

    m_bUpdate     = true;
    m_pSwapIn     = m_pSubPoints_1;
    m_pSwapOut    = m_pSubPoints_2;
}

void GLDraw::swapPosBuffer()
{
    if (m_bUpdate) {
        m_bUpdate     = false;
        m_pSwapIn     = m_pSubPoints_2;
        m_pSwapOut    = m_pSubPoints_1;
    } else {
        m_bUpdate     = true;
        m_pSwapIn     = m_pSubPoints_1;
        m_pSwapOut    = m_pSubPoints_2;
    }
}

GLfloat GLDraw::gaussian(GLfloat X)
{
    GLfloat mu = MU;
    GLfloat div = X - mu;
    GLfloat sigma = SIGMA;
    GLfloat pi = M_PI;
    GLfloat mul = 2 * pi * sigma * sigma;
    GLfloat m = div / sigma;
    GLfloat s = sqrt(mul);
    GLfloat e = exp(-0.5f * m * m);

    GLfloat result = (1 / s) * e;
    return result;
}

GLfloat GLDraw::CDF(int i, int n)
{
    int num         = 1000;

    GLfloat MIN     = -5;
    GLfloat MAX     = 5;
    GLfloat d       = (MAX - MIN) / n;

    GLfloat min     = MIN + (i * d);
    GLfloat max     = MIN + ((i + 1) * d);
    GLfloat subd    = d / num;

    GLfloat x = 0.0;
    GLfloat y = 0.0;

    GLfloat result = 0.0;

    for (int j = 0; j < 1000; j++) {
        assert(max >= (min + j * subd));
        x = min + j * subd;
        result += gaussian(x) * subd;
    }

    return result;
}

void GLDraw::gaussianSmooth(int P)
{
    GLfloat r_x     = 0.0;
    GLfloat r_y     = 0.0;
    GLfloat r_z     = 0.0;

    GLfloat tmp_x   = 0.0;
    GLfloat tmp_y   = 0.0;
    GLfloat tmp_z   = 0.0;

    GLint cycle     = CYCLE;
    GLint k         = pow(2, cycle) * SMOOTH;

    if (m_bDivEven) {
        k = pow(2, cycle)  * SMOOTH;
    } else {
        k = (pow(2, (cycle + 1)) - 1) * SMOOTH;
    }
    GLint n         = k + 1;
    GLint h         = k / 2.0f;

    vector<GLfloat> *v_weight = new vector<GLfloat>[h];

    GLfloat result = 0.0;
    for (int j = 0; j < h; j++) {
        result = 0.0;
        GLint tmp_h = j + 1;
        GLint tmp_n = 2 * tmp_h + 1;
        for (int i = -tmp_h; i <= tmp_h; i++) {
            v_weight[j].push_back(CDF(tmp_h + i, tmp_n));
            result += v_weight[j][tmp_h + i];
        }
    }

    int size0 = 0;
    int size = P + 2;
    for (int i = 0; i < m_vLines.size(); i++) {
        // source
        m_pSmoothPoints[size0++] = m_pSwapIn[m_nMax * i * 3 + 0];
        m_pSmoothPoints[size0++] = m_pSwapIn[m_nMax * i * 3 + 1];
        m_pSmoothPoints[size0++] = m_pSwapIn[m_nMax * i * 3 + 2];

        for (int j = 0; j < size; j++) {
            if (j == 0 || j == (size - 1)) {
                continue;
            }

            r_x = 0.0;
            r_y = 0.0;
            r_z = 0.0;

            if ((j - h) < 0) {
                GLint tmp_h = j - 0;
                GLint tmp_w = tmp_h - 1;
                for (int m = -tmp_h; m <= tmp_h; m++) {
                    tmp_x = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 0];
                    tmp_y = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 1];
                    tmp_z = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 2];

                    r_x += v_weight[tmp_w][tmp_h + m] * tmp_x;
                    r_y += v_weight[tmp_w][tmp_h + m] * tmp_y;
                    r_z += v_weight[tmp_w][tmp_h + m] * tmp_z;
                }
            } else if ((j + h) > (size - 1)) {
                GLint tmp_h = size - 1 - j;
                GLint tmp_w = tmp_h - 1;
                for (int m = -tmp_h; m <= tmp_h; m++) {
                    tmp_x = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 0];
                    tmp_y = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 1];
                    tmp_z = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 2];

                    r_x += v_weight[tmp_w][tmp_h + m] * tmp_x;
                    r_y += v_weight[tmp_w][tmp_h + m] * tmp_y;
                    r_z += v_weight[tmp_w][tmp_h + m] * tmp_z;
                }
            } else {
                for (int m = -h; m <= h; m++) {
                    GLint tmp_w = h - 1;
                    tmp_x = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 0];
                    tmp_y = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 1];
                    tmp_z = m_pSwapIn[m_nMax * i * 3 + (j + m) * 3 + 2];

                    r_x += v_weight[tmp_w][h + m] * tmp_x;
                    r_y += v_weight[tmp_w][h + m] * tmp_y;
                    r_z += v_weight[tmp_w][h + m] * tmp_z;
                }
            }
            m_pSmoothPoints[size0++] = r_x;
            m_pSmoothPoints[size0++] = r_y;
            m_pSmoothPoints[size0++] = r_z;
        }

        // target
        m_pSmoothPoints[size0++]
                = m_pSwapIn[m_nMax * i * 3  + (P + 1) * 3 + 0];
        m_pSmoothPoints[size0++]
                = m_pSwapIn[m_nMax * i * 3  + (P + 1) * 3 + 1];
        m_pSmoothPoints[size0++]
                = m_pSwapIn[m_nMax * i * 3  + (P + 1) * 3 + 2];
    }

    assert(size0 == m_nLineNum * m_nMax * 3);

    delete[] v_weight;
    v_weight = NULL;
}

void GLDraw::saveTexBuffer()
{
    GLfloat *p_pixel = new GLfloat[4 * WIDTH * HEIGHT];

    int n_count = 0;
    for (int i = 0; i < WIDTH * HEIGHT; i++) {
        if (p_pixel[4 * i + 0] != 0.0 ||
            p_pixel[4 * i + 1] != 0.0 ||
            p_pixel[4 * i + 2] != 0.0) {
            n_count++;
        }
    }

    string file;
    if (m_bCurveShaders) {
//        file = "../ONE_EDGE";
//        file = "../ONE_TEST_EDGE";
        file = "../MULTIPLE_EDGE_800_800";
//        file = "../ENTROPY_GRAPH/AIRLINE_EDGES_FDEB_800";
//        file = "../ENTROPY_GRAPH/MIGRATIONS_EDGES_FDEB_800";
//        file = "../ENTROPY_GRAPH/Edison_Mean_Shift/THREE_EDGE_50";
//        file = "../ENTROPY_GRAPH/Edison_Mean_Shift/TEST_EDGE_50";
//        file = "../DATA_EDGE";
//        file = "../ENTROPY_GRAPH/THREE_EDGES_800";
//        file = "../ENTROPY_GRAPH/TEST_EDGES_800";
//        file = "../ENTROPY_GRAPH/ONE_TEST_EDGES_200";
    } else {
//        file = "../ONE_NODE";
//        file = "../ONE_TEST_NODE";
        file = "../MULTIPLE_NODE_800_800";
//        file = "../ENTROPY_GRAPH/AIRLINE_NODES_FDEB_800";
//        file = "../ENTROPY_GRAPH/MIGRATIONS_NODES_FDEB_800";
//        file = "../ENTROPY_GRAPH/Edison_Mean_Shift/THREE_NODE_50";
//        file = "../ENTROPY_GRAPH/Edison_Mean_Shift/TEST_NODE_50";
//        file = "../DATA_NODE";
//        file = "../ENTROPY_GRAPH/THREE_NODES_800";
//        file = "../ENTROPY_GRAPH/TEST_NODES_800";
//        file = "../ENTROPY_GRAPH/ONE_TEST_NODES_400";
    }

    ofstream output(file.c_str());

    if (!output) {
        cout << "Cannot open " << file << endl;
    }

    int width = WIDTH;
    int height = HEIGHT;

    output.write(reinterpret_cast<char *>(p_pixel),
                 sizeof(float) * 4 * width * height);

    output.close();    

    cout << "Save " << file << endl;
}

void GLDraw::bundleEdges()
{
    initVertex();

    initAdjMatrix();

    initDistance();

    updateBuffer();

    if (m_nSegment == 1) {
        m_nSegment++;
    }

    initBuffer();    
}

