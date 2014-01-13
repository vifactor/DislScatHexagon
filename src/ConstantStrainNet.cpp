/*
 * ConstantStrainNet.cpp
 *
 *  Created on: 13 ñ³÷. 2014
 *      Author: kopp
 */

#include "ConstantStrainNet.h"

ConstantStrainNet::ConstantStrainNet(double eps_xz, double eps_yz, double eps_zz, double d) : MCInterface(d)
{
	m_eps_xz = eps_xz;
	m_eps_yz = eps_yz;
	m_eps_zz = eps_zz;
}

ConstantStrainNet::~ConstantStrainNet()
{
}

const Geometry::Vector3d& ConstantStrainNet::u(const Geometry::Vector3d& r) const
{
	/*do not consider other displacement components (symmetric peak)*/
	m_u.x = 0.0;
	m_u.y = 0.0;
	/*linear dependence of displacements from strain fields*/
	m_u.z = 2 * m_eps_xz * r.x + 2 * m_eps_yz * r.y + m_eps_zz * r.z;

	return m_u;
}
