#include "MisfitSet.h"
//#include <iostream>
double MisfitSet::m_width = 0.0;
bool MisfitSet::m_isPeriodic = true;

using namespace Geometry;

MisfitSet::MisfitSet(Distribution1d * distribution, double depth, double rho,
		const Vector3d& normal, const Vector3d& Burgers, const Vector3d& Line,
		double nu)
{
	m_distribution = distribution;
	setDislocations(depth, rho, normal, Burgers, Line, nu);
}

MisfitSet::~MisfitSet()
{
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        delete m_dislocations[i];
    }
}

const Vector3d& MisfitSet::U(const Vector3d& r) const
{
    m_u.set(0.0, 0.0, 0.0);
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        m_u += m_dislocations[i]->U(r);
    }

    return m_u;
}

void MisfitSet::update()
{
    static double pos;
    static double step, shift;
    static int sx, sy, sz;
    //double angle;

    /*sign of the misfit component does not change*/
    sx = 1;

    /*initial position is shifted by a mean distance to the left*/
    step = m_distribution->getStart();
    pos = -m_width/2 - step;

    //angle = gsl_ran_gaussian(m_distribution->getRng(), M_PI / 75);
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        /** setup new position */
    	if(m_isPeriodic)
    	{
    		shift = 0;
    		step = m_distribution->getMean();
    	}
    	else
    	{
            shift = m_distribution->getShift();
            step = m_distribution->getStep();
    	}
        pos += step;

        /** flip randomly by and bz components */
        sy = 1 - 2 * gsl_rng_uniform_int (m_distribution->getRng(), 2);
        sz = 1 - 2 * gsl_rng_uniform_int (m_distribution->getRng(), 2);

        /*set dislocation in the new position*/
        m_dislocations[i]->moveTo(pos + shift);
        /*flip randomly signs of nonmisfit components*/
        m_dislocations[i]->multBurgers(sx, sy, sz);
        /*rotate dislocation line*/
        //m_dislocations[i]->rotate(angle);
    }
}

void MisfitSet::setDislocations(double depth, double rho, const Vector3d& normal,
                           const Vector3d& Burgers, const Vector3d& Line, double nu)
{
    size_t nb;
    double mean;

    /*nb of dislocations*/

    nb = rho * m_width;
    /*mean distance between dislocations*/
    if(nb > 1)
        mean = m_width / (nb - 1);
    else
        mean = m_width;

    m_distribution->setMean(mean);

    m_dislocations.resize(nb);
    for(size_t i = 0; i < m_dislocations.size(); ++i)
    {
        m_dislocations[i] = new DislocationParallel(Burgers, Line, normal, depth, nu);
    }
}
