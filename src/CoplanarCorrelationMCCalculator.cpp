/*
 * CoplanarCorrelationMCCalculator.cpp
 *
 *  Created on: 8 nov. 2013
 *      Author: kopp
 */

#include "CoplanarCorrelationMCCalculator.h"

using namespace Geometry;

CoplanarCorrelationMCCalculator::CoplanarCorrelationMCCalculator(MCSample * sample, const gsl_rng * rng, const Geometry::Vector3d& Q) : MCCalculator(sample)
{
	m_rng = rng;
    m_Q = Q;

    /* we need to setup the sample in child constructor,
     *  as the sample is required for laboratory frame initialization
     */
    m_sample = sample;

    setupLaboratoryFrame();
}

CoplanarCorrelationMCCalculator::~CoplanarCorrelationMCCalculator()
{
}

void CoplanarCorrelationMCCalculator::run(MCData * data)
{
	static double x, z1, z2, QU;
	static Vector3d r1, r2, U1, U2;

	for(unsigned long int istep = 0; istep < data->getNbSteps(); ++istep)
	{
		m_sample->update();
		for(size_t i = 0; i < data->getNbPoints(); ++i)
		{
			x = data->arg(i, m_idx_x);
			z1 = data->arg(i, m_idx_z1);
			z2 = data->arg(i, m_idx_z2);

		    /*find coordinates in the lab coord system*/
		    /* minus sign before z-coordinates is due to the chosen laboratory coordinate frame*/
		    r1 = m_axis_x * x / 2 + m_axis_z * z1;
		    r2 = m_axis_x * (-x / 2) + m_axis_z * z2;
		    U1 = m_sample->u(r1);
		    U2 = m_sample->u(r2);

		    QU = m_Q * (U1 - U2);

		    //std::cout << x << "\t" << z1 << "\t" << z2 << "\t" << QU << std::endl;
		    //std::cout << m_Q.x << "\t" << m_Q.y << "\t" << m_Q.z << std::endl;

	        data->func(i, m_idx_reG) += cos(QU);
	        data->func(i, m_idx_imG) += sin(QU);
		}
	}
}

void CoplanarCorrelationMCCalculator::setupLaboratoryFrame()
{
    Geometry::Vector3d Qpar, Qper;

    Qpar = m_sample->normal() * (m_sample->normal() * m_Q);
    Qper = m_Q - Qpar;

    /*z along surface normal*/
    m_axis_z = m_sample->normal();
    m_axis_z.normalize();
    if(Qper.norm() <= m_epsilon)
    {
        /*for symmetric reflection orientation of x-axis is irrelevant*/
        /*find any vector perpendicular to z */
		Geometry::Vector3d vec;

		/*
		 * calculate the cross product [m_axis_z % vec],
		 * where vec is not collinear with m_axis_z
		 */
		if((m_axis_z.y != 0) || (m_axis_z.z != 0))
			m_axis_x = m_axis_z % Geometry::Vector3d(1, 0, 0);
		else
			m_axis_x = m_axis_z % Geometry::Vector3d(0, 1, 0);

		m_axis_x.normalize();
    }
    else
    {
        /*for asymmetric reflection orientation x along Qper*/
        m_axis_x = Qper.normalize();
    }
    /* y perpendicular to both */
    m_axis_y = m_axis_z % m_axis_x;

    /*std::cout << "Qx:\t" << m_Q * m_axis_x << std::endl;
    std::cout << "Qz:\t" << m_Q * m_axis_z << std::endl;*/
}

