/*
 * ConstantStrainNet.h
 *
 *  Created on: 13 ñ³÷. 2014
 *      Author: kopp
 */

#ifndef CONSTANTSTRAINNET_H_
#define CONSTANTSTRAINNET_H_

#include "MCAbstractClasses.h"

class ConstantStrainNet: public MCInterface
{
public:
	ConstantStrainNet(double eps_xz, double eps_yz, double eps_zz, double d);
	virtual ~ConstantStrainNet();
	virtual void update() {/*field distribution is static*/}
	virtual const Geometry::Vector3d& u(const Geometry::Vector3d& r) const;
	virtual size_t nbDefects() const {return 0;}
protected:
	double m_eps_xz, m_eps_yz, m_eps_zz;
};

#endif /* CONSTANTSTRAINNET_H_ */
