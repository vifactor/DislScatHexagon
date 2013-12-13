/*
 * Distribution1d.cpp
 *
 *  Created on: 2 זמגע. 2013
 *      Author: kopp
 */

#include "Distribution1d.h"


Distribution1d::Distribution1d(const gsl_rng * rng)
{
	m_rng = rng;
	m_mean = 1;
}

PeriodicDistribution1d::PeriodicDistribution1d(const gsl_rng * rng, double mean): Distribution1d(rng)
{
	m_mean = mean;
}

void PeriodicDistribution1d::setMean(double mean)
{
	m_mean = mean;
}

double PeriodicDistribution1d::getStep() const
{
	return m_mean;
}

double PeriodicDistribution1d::getShift() const
{
	return 0.0;
}

GammaDistribution1d::GammaDistribution1d(const gsl_rng * rng, int gamma, double mean): Distribution1d(rng)
{
	m_gamma = gamma;
	m_mean = mean;
	m_beta = m_mean / m_gamma;
}

void GammaDistribution1d::setMean(double mean)
{
	m_mean = mean;
	m_beta = m_mean / m_gamma;
}

double GammaDistribution1d::getStep() const
{
	return gsl_ran_gamma(m_rng, m_gamma, m_beta);
}

double GammaDistribution1d::getShift() const
{
	return 0.0;
}

GaussDistribution1d::GaussDistribution1d(const gsl_rng * rng, double sigma, double mean): Distribution1d(rng)
{
	m_mean = mean;
	m_sigma = sigma;
	m_variance = m_sigma * m_mean;
}

void GaussDistribution1d::setMean(double mean)
{
	m_mean = mean;
	m_variance = m_sigma * m_mean;
}

double GaussDistribution1d::getStep() const
{
	return m_mean;
}

double GaussDistribution1d::getShift() const
{
	return gsl_ran_gaussian(m_rng, m_variance);
}

GammaGaussDistribution1d::GammaGaussDistribution1d(const gsl_rng * rng, int gamma, double sigma, double mean): Distribution1d(rng)
{
	m_mean = mean;

	m_sigma = sigma;
	m_variance = m_sigma * m_mean;
	m_gamma = gamma;
	m_beta = m_mean / m_gamma;
}

void GammaGaussDistribution1d::setMean(double mean)
{
	m_mean = mean;
	m_beta = m_mean / m_gamma;
	m_variance = m_sigma * m_mean;
}

double GammaGaussDistribution1d::getStep() const
{
	static double step;

	/*atempt to combine gamma-step and gaussian-shift distributions*/
	if(m_gamma <= 0)
	{
		step = m_mean;
	}else
	{
		step = gsl_ran_gamma(m_rng, m_gamma, m_beta);
	}
	return step;
}

double GammaGaussDistribution1d::getShift() const
{
	static double shift;

	shift =	gsl_ran_gaussian(m_rng, m_variance);

	return shift;
}
