/*
 * DisplacementCalculator.cpp
 *
 *  Created on: 11 лист. 2013
 *      Author: kopp
 */

#include "LocalDisplacementCalculator.h"

using namespace Geometry;

LocalDisplacementCalculator::LocalDisplacementCalculator(MCSample * sample) : MCCalculator(sample)
{
}

LocalDisplacementCalculator::~LocalDisplacementCalculator()
{
}

void LocalDisplacementCalculator::run(MCCalculator::MCData * data)
{
	static Vector3d r, U;

	m_sample->update();
    for(std::size_t i = 0; i < data->getNbPoints(); ++i)
    {
		r.set(data->arg(i, m_idx_x), data->arg(i, m_idx_y),
				data->arg(i, m_idx_z));

		U = m_sample->u(r);

		data->func(i, m_idx_ux) = U.x;
		data->func(i, m_idx_uy) = U.y;
		data->func(i, m_idx_uz) = U.z;
    }
	data->setNbSteps(1);
}
