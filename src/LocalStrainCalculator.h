/*
 * LocalStrainCalculator.h
 *
 *  Created on: 20 sep. 2013
 *      Author: kopp
 */

#ifndef LocalStrainCalculator_H_
#define LocalStrainCalculator_H_

#include "MCAbstractClasses.h"

class LocalStrainCalculator : public MCCalculator
{
public:
	LocalStrainCalculator(MCSample * sample, double h);
	virtual ~LocalStrainCalculator();
	virtual void run(MCData * data);
	enum {m_idx_exx, m_idx_eyy, m_idx_ezz, m_idx_exz, m_idx_exy, m_idx_eyz, m_nb_idx_func};
	enum {m_idx_x, m_idx_y, m_idx_z, m_nb_idx_arg};
protected:
	/*step in calculation of derivatives*/
	double m_h;
	/*method to calculate local strain due to defects in a sample*/
	virtual void strain(const Geometry::Vector3d& r, double& exx, double& eyy,
				double& ezz, double& exy, double& exz, double& eyz);
};

#endif /* LocalStrainCalculator_H_ */
