#ifndef MILLERINDEXHEX_H
#define MILLERINDEXHEX_H

#include "Vector3d.h"

enum {MillerHexIndicesDimension = 4};

struct MillerReciprocalHexIndices
{
	MillerReciprocalHexIndices(double h = 0, double k = 0, double i = 0, double l = 0) {H = h; K = k; I = i; L = l;}
    double H, K, I, L;

    bool isCorrect()
    {
        /**property of hexagonal Miller indices**/
        if((H + K + I) != 0)
            return false;
        else
            return true;
    }
};

struct MillerDirectHexIndices
{
	MillerDirectHexIndices(double x = 0, double y = 0, double t = 0, double z = 0) {X = x; Y = y; T = t; Z = z;}
    double X, Y, T, Z;
    bool isCorrect()
    {
        /**property of hexagonal Miller indices**/
        if((X + Y + T) != 0)
            return false;
        else
            return true;
    }
};

class MillerHexIndicesTransformator
{
public:
    MillerHexIndicesTransformator(double a, double c);
    virtual ~MillerHexIndicesTransformator();
    Geometry::Vector3d toVector3d(const MillerReciprocalHexIndices& rec_index);
    Geometry::Vector3d toVector3d(const MillerDirectHexIndices& dir_index);
protected:
    Geometry::Vector3d m_a1, m_a2, m_a3, m_cc;
    Geometry::Vector3d m_A1, m_A2, m_A3, m_CC;
};

#endif // MILLERINDEXHEX_H
