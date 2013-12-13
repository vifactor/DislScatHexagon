/*
 * MCAbstractClasses.cpp
 *
 *  Created on: 19 setp. 2013
 *      Author: kopp
 */

#include "MCAbstractClasses.h"

using namespace Geometry;

double MCInterface::m_lateral_size = 0;
double MCInterface::m_nu = 0;
Vector3d MCInterface::m_normal = Vector3d(0.0, 0.0, 1.0);

MCCalculator::MCData::MCData(std::size_t nb_pts, std::size_t dim_func, std::size_t dim_arg)
{
	m_nb_steps = 0;
	m_nb_points = nb_pts;
	m_dim_func = dim_func;
	m_dim_arg = dim_arg;

	m_func_max = new double[dim_func];
	m_func_buffer = new double[m_dim_func * m_nb_points];
    m_arg_buffer = new double[m_dim_arg * m_nb_points];
}

void MCCalculator::MCData::init()
{
    m_nb_steps = 0;
    for(std::size_t ipt = 0; ipt < m_dim_func * m_nb_points; ++ipt)
    {
        m_func_buffer[ipt] = 0.0;
    }
    for(std::size_t idim = 0; idim < m_dim_func; ++idim)
    {
    	m_func_max[idim] = 0.0;
    }
}

void MCCalculator::MCData::normalize()
{
	if(m_nb_steps != 0)
	{
	    for(std::size_t ipt = 0; ipt < m_dim_func * m_nb_points; ++ipt)
	    {
	        m_func_buffer[ipt] /= m_nb_steps;
	    }
	}
	m_nb_steps = 0;
}

void MCCalculator::MCData::transferTo(double * buff)
{
    for(std::size_t ipt = 0; ipt < m_dim_func * m_nb_points; ++ipt)
    {
        buff[ipt] = m_func_buffer[ipt];
        m_func_buffer[ipt] = 0.0;
    }
}

void MCCalculator::MCData::append(const double * buff, unsigned long long nsteps)
{
    static double N0, N1, N;
    static double coef0, coef1;
    N0 = m_nb_steps;
    N1 = nsteps;
    N = N0 + N1;
    coef0 = N0 / N;
    coef1 = 1. / N;

    for(std::size_t idim = 0; idim < m_dim_func; ++idim)
    {
        m_func_max[idim] = coef0 * m_func_buffer[idim * m_nb_points] +
                            coef1 * buff[idim * m_nb_points];
    }


    for(std::size_t ipt = 0; ipt < m_nb_points; ++ipt)
    {
    	for(std::size_t idim = 0; idim < m_dim_func; ++idim)
    	{
            m_func_buffer[ipt + idim * m_nb_points] = coef0 * m_func_buffer[ipt + idim * m_nb_points] +
                            coef1 * buff[ipt + idim * m_nb_points];

            /*renew the max func value*/
            if(fabs(m_func_buffer[ipt + idim * m_nb_points]) > fabs(m_func_max[idim]))
            	m_func_max[idim] = m_func_buffer[ipt + idim * m_nb_points];
    	}
    }
    m_nb_steps = N;
}

MCCalculator::MCData::~MCData()
{
	m_nb_steps = 0;
	m_nb_points = 0;
	delete m_func_max;
	delete m_func_buffer;
    delete m_arg_buffer;
}

MCInterface::MCInterface(double d)
{
	m_depth = d;
}

MCInterface::~MCInterface()
{
}

MCSample::MCSample(double w, double t, const Vector3d & normal)
{
	m_thickness = t;
	m_lateral_size = w;
	m_normal = normal;
}

void MCSample::update()
{
	for(std::size_t i = 0; i < m_interfaces.size(); ++i)
	{
		m_interfaces[i]->update();
	}
}

size_t MCSample::nbDefects() const
{
	static size_t n;

	n = 0;
	for(std::size_t i = 0; i < m_interfaces.size(); ++i)
	{
		n = m_interfaces[i]->nbDefects();
	}

	return n;
}

const Vector3d& MCSample::u(const Geometry::Vector3d& r) const
{
	m_u = Vector3d(0.0, 0.0, 0.0);

	for(std::size_t i = 0; i < m_interfaces.size(); ++i)
	{
		m_u += m_interfaces[i]->u(r);
	}

	return m_u;
}

void MCSample::addInterface(MCInterface * interface)
{
	m_interfaces.push_back(interface);
}
