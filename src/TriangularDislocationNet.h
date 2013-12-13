/*
 * TriangularDislocationNet.h
 *
 *  Created on: 30 вер. 2013
 *      Author: kopp
 */

#ifndef TRIANGULARDISLOCATIONNET_H_
#define TRIANGULARDISLOCATIONNET_H_

#include "MCAbstractClasses.h"
#include "MisfitSet.h"

class TriangularDislocationNet: public MCInterface
{
public:
	/*properties of the dislocation net*/
	TriangularDislocationNet(Distribution1d * distribution, double frac, const Geometry::Vector3d & burgers,
			const Geometry::Vector3d& line, double rho, double d);
	virtual ~TriangularDislocationNet();
	virtual size_t nbDefects() const {return m_nbDislocations;}
	virtual void update();
	virtual const Geometry::Vector3d& u(const Geometry::Vector3d& r) const;
protected:
	size_t m_nbDislocations;
	MisfitSet * m_set0, * m_set1,  * m_set2;
	Distribution1d * m_distribution;
	double m_frac;
};

#endif /* TRIANGULARDISLOCATIONNET_H_ */
