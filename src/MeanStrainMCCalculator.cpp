/*
 * MeanStrainMCCalculator.cpp
 *
 *  Created on: 24 вер. 2013
 *      Author: kopp
 */

#include "MeanStrainMCCalculator.h"

using namespace Geometry;

MeanStrainMCCalculator::MeanStrainMCCalculator(MCSample * sample, const gsl_rng * rng, double h) : MCCalculator(sample)
{
	m_h = h;
	m_rng = rng;
}

MeanStrainMCCalculator::~MeanStrainMCCalculator()
{
}

void MeanStrainMCCalculator::run(MCData * data)
{
	static double exx, eyy, ezz, exy, exz, eyz;
	static Vector3d r;

	/*equivalence of averaging over ensemble and space*/
	r.x = 0;
	r.y = 0;
	for(unsigned long int istep = 0; istep < data->getNbSteps(); ++istep)
	{
		m_sample->update();
		for(size_t i = 0; i < data->getNbPoints(); ++i)
		{
			r.z = data->arg(i, m_idx_z);

			strain(r, exx, eyy, ezz, exy, exz, eyz);

			data->func(i, m_idx_exx) += exx;
			data->func(i, m_idx_eyy) += eyy;
			data->func(i, m_idx_ezz) += ezz;
			data->func(i, m_idx_exy) += exy;
			data->func(i, m_idx_exz) += exz;
			data->func(i, m_idx_eyz) += eyz;
			data->func(i, m_idx_e2xx) += exx * exx;
			data->func(i, m_idx_e2yy) += eyy * eyy;
			data->func(i, m_idx_e2zz) += ezz * ezz;
			data->func(i, m_idx_e2xy) += exy * exy;
			data->func(i, m_idx_e2xz) += exz * exz;
			data->func(i, m_idx_e2yz) += eyz * eyz;
		}
	}
	/*normalization: division by nb of steps performed*/
	//data->normalize();
}

void MeanStrainMCCalculator::strain(const Geometry::Vector3d& r, double& exx, double& eyy,
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
