/*
 * HexagonalDislocationNet.cpp
 *
 *  Created on: 23 nov. 2012
 *      Author: kopp
 */

#include "HexagonalDislocationNet.h"
#include <iostream>

using namespace Geometry;

HexagonalDislocationNet::HexagonalDislocationNet(HexagonalNet * distribution, const Geometry::Vector3d& Burgers, double d, double nu) : MCInterface(d)
{
	size_t i;
	EdgeIter ei, ei_end;

	m_distribution = distribution;

	m_defects.resize(num_edges(*m_distribution));

	i = 0;
	for (tie(ei, ei_end) = edges(*m_distribution); ei != ei_end; ++ei)
	{
		const Vector2d& p1 = (*m_distribution)[ei->m_source].m_disturbedCoord;
		const Vector2d& p2 = (*m_distribution)[ei->m_target].m_disturbedCoord;

		m_defects[i] = new DislocationPiShape(Burgers, p1, p2, m_normal, d, nu);
		++i;
	}
}

HexagonalDislocationNet::~HexagonalDislocationNet()
{
	for(size_t i = 0; i < m_defects.size(); ++i)
	{
		delete m_defects[i];
	}
}

void HexagonalDislocationNet::update()
{
	size_t i;
	EdgeIter ei, ei_end;

	m_distribution->update();
	i = 0;
	for (tie(ei, ei_end) = edges(*m_distribution); ei != ei_end; ++ei)
	{
		const Vector2d& p1 = (*m_distribution)[ei->m_source].m_disturbedCoord;
		const Vector2d& p2 = (*m_distribution)[ei->m_target].m_disturbedCoord;

		m_defects[i]->reset(p1, p2);
		++i;
	}
}

const Vector3d& HexagonalDislocationNet::u(const Vector3d& r) const
{
	m_u.set(0.0, 0.0, 0.0);
	for(size_t i = 0; i < m_defects.size(); ++i)
	{
		m_u += m_defects[i]->U(r);
	}
	return m_u;
}
