/*
 * MCAbstractClasses.h
 *
 *  Created on: 19 sep. 2013
 *      Author: kopp
 */

#ifndef MCABSTRACTCLASSES_H_
#define MCABSTRACTCLASSES_H_

#include "Vector3d.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>

class MCInterface
{
public:
	MCInterface(double d);
	virtual ~MCInterface();
	virtual size_t nbDefects() const = 0;
	virtual void update() = 0; /*pure virtual method to update a defect distribution*/
	virtual const Geometry::Vector3d& u(const Geometry::Vector3d& r) const = 0; /*pure virtual method to calculate displacements due to defects*/
	static double m_lateral_size;
	/*poisson ratio*/
	static double m_nu;
	static Geometry::Vector3d m_normal;
protected:
	/*at which depth from the surface the interface is situated*/
	double m_depth;
	/*handy variable to return displacement*/
	mutable Geometry::Vector3d m_u;
};

class MCSample
{
public:
	MCSample(double w, double d, const Geometry::Vector3d& normal = Geometry::Vector3d(0.0, 0.0, 1.0));
	virtual ~MCSample() {}
	virtual void addInterface(MCInterface * interface);
	virtual size_t nbDefects() const;
	virtual double thickness() const {return m_thickness;};
	virtual const Geometry::Vector3d & normal() const {return m_normal;}
	virtual void update(); /*pure virtual method to update a defect distribution*/
	virtual const Geometry::Vector3d& u(const Geometry::Vector3d& r) const; /*method to calculate displacements due to defects*/
protected:
	/*normal to the sample surface*/
	Geometry::Vector3d m_normal;
	/*collection of interfaces*/
	std::vector<MCInterface * > m_interfaces;
	/*spacial parameters of the sample*/
	double m_thickness, m_lateral_size;
	/*handy variable to return displacement*/
	mutable Geometry::Vector3d m_u;
};

class MCCalculator
{
public:
    class MCData
    {
    public:
        MCData(std::size_t nb_pts, std::size_t dim_func = 1, std::size_t dim_arg = 1);
        virtual ~MCData();
        void init();
        void normalize();
        double& func(std::size_t i, std::size_t ifunc = 0) {return m_func_buffer[i + ifunc * m_nb_points];}
        const double& func(std::size_t i, std::size_t ifunc = 0) const {return m_func_buffer[i + ifunc * m_nb_points];}
        double& arg(std::size_t i, std::size_t iarg = 0) {return m_arg_buffer[i + iarg * m_nb_points];}
        const double& arg(std::size_t i, std::size_t iarg = 0) const {return m_arg_buffer[i + iarg * m_nb_points];}
        std::size_t getNbPoints() const {return m_nb_points;}
        std::size_t getFuncDim() const {return m_dim_func;}
        std::size_t getArgDim() const {return m_dim_arg;}
        unsigned long long getNbSteps() const {return m_nb_steps;}
        void setNbSteps(unsigned long long nb) {m_nb_steps = nb;}
        double max(std::size_t ifunc) const {return m_func_max[ifunc];}
        /*performs the MC summation func_res = (func_init + func_add) / (N_init + N_add)*/
        void append(const double * func_buff, unsigned long long nsteps);
        /*copies the content of m_func_buff to func_buff and nullifies m_func_buffer*/
        void transferTo(double * func_buff);
    protected:
        unsigned long long m_nb_steps;
        double * m_func_max;
        std::size_t m_nb_points;
        double * m_func_buffer;
        double * m_arg_buffer;
        std::size_t m_dim_func, m_dim_arg;
    };

    /** Default constructor */
    MCCalculator(MCSample * sample){m_sample = sample;}
    /** Default destructor */
    virtual ~MCCalculator(){}
    /*main calculation loop*/
    virtual void run(MCData * data) = 0;
protected:
    MCSample * m_sample;
};

#endif /* MCABSTRACTCLASSES_H_ */
