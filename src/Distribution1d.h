/*
 * Distribution1d.h
 *
 *  Created on: 2 זמגע. 2013
 *      Author: kopp
 */

#ifndef DISTRIBUTION1D_H_
#define DISTRIBUTION1D_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>

class Distribution1d
{
public:
	Distribution1d(const gsl_rng * rng);
	virtual ~Distribution1d(){}
	virtual double getStep() const = 0;
	virtual double getShift() const  = 0;
	virtual double getMean() const {return m_mean;}
	virtual void setMean(double mean) {m_mean = mean;}
	virtual double getStart() const {return gsl_ran_flat(m_rng, 0, m_mean);}
	virtual const gsl_rng * getRng() const {return m_rng;}
protected:
	double m_mean;
	const gsl_rng * m_rng;
};

class PeriodicDistribution1d: public Distribution1d
{
public:
	PeriodicDistribution1d(const gsl_rng * rng, double mean = 1.0);
	virtual ~PeriodicDistribution1d() {}
	virtual double getStep() const;
	virtual double getShift() const;
	virtual void setMean(double mean);
protected:
};

class GammaDistribution1d: public Distribution1d
{
public:
	GammaDistribution1d(const gsl_rng * rng, int gamma, double mean = 1.0);
	virtual ~GammaDistribution1d() {}
	virtual double getStep() const;
	virtual double getShift() const;
	virtual void setMean(double mean);
protected:
	double m_beta;
	int m_gamma;
};

class GaussDistribution1d: public Distribution1d
{
public:
	GaussDistribution1d(const gsl_rng * rng, double sigma, double mean = 1.0);
	virtual ~GaussDistribution1d() {}
	double getStep() const;
	double getShift() const;
	virtual void setMean(double mean);
protected:
	double m_sigma;
	double m_variance;
};

class GammaGaussDistribution1d: public Distribution1d
{
public:
	GammaGaussDistribution1d(const gsl_rng * rng, int gamma, double sigma, double mean = 1.0);
	virtual ~GammaGaussDistribution1d() {}
	double getStep() const;
	double getShift() const;
	virtual void setMean(double mean);
protected:
	double m_beta;
	int m_gamma;
	double m_sigma;
	double m_variance;
};

#endif /* DISTRIBUTION1D_H_ */
