/*
 * HexagonalDislocationNet.h
 *
 *  Created on: 23 nov. 2012
 *      Author: kopp
 */

#ifndef HEXAGONALDISLOCATIONNET_H_
#define HEXAGONALDISLOCATIONNET_H_

#include "HexagonalNet.h"
#include "MCAbstractClasses.h"
#include "DislocationClasses.h"
#include <vector>

#include <fstream>

class HexagonalDislocationNet : public MCInterface
{
public:
	/*properties of the dislocation net*/
	HexagonalDislocationNet(HexagonalNet * distribution, const Geometry::Vector3d& Burgers, double d, double nu);
	virtual ~HexagonalDislocationNet();
	virtual size_t nbDefects() const {return m_defects.size();}
	virtual void update();
	virtual const Geometry::Vector3d& u(const Geometry::Vector3d& r) const;
protected:
	HexagonalNet * m_distribution;
	std::vector<DislocationPiShape * > m_defects;
	/*Burgers vector*/
	Geometry::Vector3d m_burgers;
};

#endif /* HEXAGONALDISLOCATIONNET_H_ */
