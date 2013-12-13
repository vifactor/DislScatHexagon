/*
 * MeanStrainMCCalculator.h
 *
 *  Created on: 24 вер. 2013
 *      Author: kopp
 */

#ifndef MEANSTRAINMCCALCULATOR_H_
#define MEANSTRAINMCCALCULATOR_H_

#include "MCAbstractClasses.h"

class MeanStrainMCCalculator: public MCCalculator
{
public:
	MeanStrainMCCalculator(MCSample * sample, const gsl_rng * rng, double h);
	virtual ~MeanStrainMCCalculator();
	virtual void run(MCData * data);
	enum{m_idx_exx, m_idx_eyy, m_idx_ezz,
		m_idx_exz, m_idx_exy, m_idx_eyz,
		m_idx_e2xx, m_idx_e2yy,	m_idx_e2zz,
		m_idx_e2xz,	m_idx_e2xy,	m_idx_e2yz,
		m_nb_idx_func};
	enum {m_idx_z, m_nb_idx_arg};
protected:
	/*random number generator*/
	const gsl_rng * m_rng;
	/*step in calculation of derivatives*/
	double m_h;
	/*method to calculate local strain due to defects in a sample*/
	virtual void strain(const Geometry::Vector3d& r, double& exx, double& eyy,
				double& ezz, double& exy, double& exz, double& eyz);
};

#endif /* MEANSTRAINMCCALCULATOR_H_ */
