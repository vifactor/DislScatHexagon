/*
 * CoplanarCorrelationMCCalculator.h
 *
 *  Created on: 8 nov. 2013
 *      Author: kopp
 */

#ifndef COPLANARCORRELATIONMCCALCULATOR_H_
#define COPLANARCORRELATIONMCCALCULATOR_H_

#include "MCAbstractClasses.h"

class CoplanarCorrelationMCCalculator: public MCCalculator
{
public:
	CoplanarCorrelationMCCalculator(MCSample * sample, const gsl_rng * rng, const Geometry::Vector3d& Q);
	virtual ~CoplanarCorrelationMCCalculator();
	virtual void run(MCData * data);
	enum{m_idx_reG, m_idx_imG,
		m_nb_idx_func};
	enum {m_idx_x, m_idx_z1, m_idx_z2,  m_nb_idx_arg};
protected:
	/*random number generator*/
	const gsl_rng * m_rng;

    void setupLaboratoryFrame();
    Geometry::Vector3d m_Q;
    Geometry::Vector3d m_axis_x, m_axis_y, m_axis_z;
    const static double m_epsilon = 1e-10;
};

#endif /* COPLANARCORRELATIONMCCALCULATOR_H_ */
