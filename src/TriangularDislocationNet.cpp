/*
 * TriangularDislocationNet.cpp
 *
 *  Created on: 30 вер. 2013
 *      Author: kopp
 */

#include "TriangularDislocationNet.h"

using namespace Geometry;

TriangularDislocationNet::TriangularDislocationNet(Distribution1d * distribution, double frac,
		const Vector3d & burgers, const Vector3d& line, double rho, double d) :
		MCInterface(d)
{
	Vector3d iLine;

	MisfitSet::m_width = MCInterface::m_lateral_size;
	m_frac = frac;

	m_distribution = distribution;

	m_nbDislocations = 0;

	/*first dislocation set*/
	iLine = line;
	m_set0 = new MisfitSet(distribution, d, rho, MCInterface::m_normal,
			burgers,
			iLine, m_nu);
	m_nbDislocations += m_set0->getNbDislocations();

	/*second dislocation set (rotated by 120 deg)*/
	iLine.RotateAboutAxis(2 * M_PI / 3, MCInterface::m_normal);
	m_set1 = new MisfitSet(distribution, d, rho, MCInterface::m_normal,
			burgers,
			iLine, m_nu);
	m_nbDislocations += m_set1->getNbDislocations();

	/*third dislocation set (rotated by 240 deg)*/
	iLine.RotateAboutAxis(2 * M_PI / 3, MCInterface::m_normal);
	m_set2 = new MisfitSet(distribution, d, rho, MCInterface::m_normal,
			burgers,
			iLine, m_nu);
	m_nbDislocations += m_set2->getNbDislocations();
}

TriangularDislocationNet::~TriangularDislocationNet()
{
	delete m_set0;
	delete m_set1;
	delete m_set2;
}

const Vector3d& TriangularDislocationNet::u(const Vector3d& r) const
{
	//m_u.set(0,0,0);
	m_u = m_set0->U(r);
	m_u += m_set1->U(r);
	m_u += m_set2->U(r);

	return m_u;
}

void TriangularDislocationNet::update()
{
	static double r;

	r = gsl_rng_uniform(m_distribution->getRng());
	if(r < m_frac)
		MisfitSet::m_isPeriodic = true;
	else
		MisfitSet::m_isPeriodic = false;

	m_set0->update();
	m_set1->update();
	m_set2->update();
}

