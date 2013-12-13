/*
 * LocalStrainCalculator.cpp
 *
 *  Created on: 20 вер. 2013
 *      Author: kopp
 */

#include "LocalStrainCalculator.h"
#include <iostream>

using namespace Geometry;

LocalStrainCalculator::LocalStrainCalculator(MCSample * sample, double h) : MCCalculator(sample)
{
	m_h = h;
}

LocalStrainCalculator::~LocalStrainCalculator()
{
}

void LocalStrainCalculator::run(MCCalculator::MCData * data)
{
	static Vector3d r;
	m_sample->update();
    for(std::size_t i = 0; i < data->getNbPoints(); ++i)
    {
		r.set(data->arg(i, m_idx_x), data->arg(i, m_idx_y),
				data->arg(i, m_idx_z));

        strain(r, data->func(i, m_idx_exx),
        		data->func(i, m_idx_eyy),
        		data->func(i, m_idx_ezz),
        		data->func(i, m_idx_exy),
        		data->func(i, m_idx_exz),
        		data->func(i, m_idx_eyz));
    }
	data->setNbSteps(1);
}

void LocalStrainCalculator::strain(const Geometry::Vector3d& r, double& exx, double& eyy,
		double& ezz,  double& exy, double& exz, double& eyz)
{
	static Vector3d u, u_xh, u_yh, u_zh;
	static Vector3d r_xh, r_yh, r_zh;

	r_xh = r; r_xh.x += m_h;
	r_yh = r; r_yh.y += m_h;
	r_zh = r; r_zh.z += m_h;

	u = m_sample->u(r);
	u_xh = m_sample->u(r_xh);
	u_yh = m_sample->u(r_yh);
	u_zh = m_sample->u(r_zh);

	exx = (u_xh.x - u.x) / m_h;
	eyy = (u_yh.y - u.y) / m_h;
	ezz = (u_zh.z - u.z) / m_h;
	exy = 0.5 * ((u_yh.x - u.x) / m_h + (u_xh.y - u.y) / m_h);
	exz = 0.5 * ((u_zh.x - u.x) / m_h + (u_xh.z - u.z) / m_h);
	eyz = 0.5 * ((u_zh.y - u.y) / m_h + (u_yh.z - u.z) / m_h);
}
