#ifndef MISFITSet_H
#define MISFITSet_H

#include "DislocationClasses.h"
#include "Distribution1d.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>

class MisfitSet
{
public:
    /** Constructor */
    MisfitSet(Distribution1d * distribution, double depth, double rho, const Geometry::Vector3d& normal,
                           const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, double nu);
    /** Destructor */
    virtual ~MisfitSet();
    void update();
    const Geometry::Vector3d& U(const Geometry::Vector3d& r) const;
    size_t getNbDislocations() const {return m_dislocations.size();}
    static double m_width;
    static bool m_isPeriodic;
protected:
    void setDislocations(double depth, double rho, const Geometry::Vector3d& normal,
                           const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, double nu);

    std::vector<DislocationParallel *> m_dislocations;
    mutable Geometry::Vector3d m_u;

    double m_d;

    Distribution1d * m_distribution;
};

#endif // MISFITSet_H
