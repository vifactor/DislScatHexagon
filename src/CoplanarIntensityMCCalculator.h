#ifndef CoplanarIntensityMCCalculator_H
#define CoplanarIntensityMCCalculator_H

#include "MCAbstractClasses.h"

class CoplanarIntensityMCCalculator : public MCCalculator
{
public:
    CoplanarIntensityMCCalculator(MCSample * sample, gsl_rng * rng, const Geometry::Vector3d& Q, double sigmax, double sigmaz);
    void run(MCData * data);
    virtual ~CoplanarIntensityMCCalculator();
	enum {m_idx_reI, m_idx_imI, m_nb_idx_func};
	enum {m_idx_qx, m_idx_qz, m_nb_idx_arg};
protected:
    gsl_rng * m_rng;
    double m_sigma_x, m_sigma_z;

    void setupLaboratoryFrame();
    Geometry::Vector3d m_Q;
    Geometry::Vector3d m_axis_x, m_axis_y, m_axis_z;
    const static double m_epsilon = 1e-10;

    size_t m_steps_to_update;
    void add(MCData * data);
};

#endif // CoplanarIntensityMCCalculator_H
