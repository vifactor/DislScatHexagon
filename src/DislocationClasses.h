/*
 * Dislocations.h
 *
 *  Created on: 1 august 2012
 *  Modified on: 22 nov. 2012
 *      Author: dreamcatcher
 *
 *      This header file contains interfaces to functions
 *    	which calculate displacements due to various kinds of dislocations
 */

#ifndef DISLOCATIONS_H_
#define DISLOCATIONS_H_

#include <gsl/gsl_math.h>
#include "Vector3d.h"
#include "dislocation_displacement.h"

 class DislocationParallel
{
public:
    /** Default constructor */
    DislocationParallel(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, const Geometry::Vector3d& normal, double d, double nu);
    /** Default destructor */
    virtual ~DislocationParallel();
	void moveTo(const double pos);
	/*rotate dislocation by <angle> around z-axis*/
	void rotate(const double angle);
	void multBurgers(int sx, int sy, int sz);

	const Geometry::Vector3d& U(const Geometry::Vector3d& r) const;
protected:
    void setup(const Geometry::Vector3d& Line, const Geometry::Vector3d& normal);

    /*along Burgers vector*/
	Geometry::Vector3d m_axis_x;
	/*perpendicular to x and z*/
	Geometry::Vector3d m_axis_y;
	/*along surface normal/dislocation line*/
	Geometry::Vector3d m_axis_z;

	/*position on the interface*/
	double m_pos;

	/*Burgers vector with respect to dislocation coord frame*/
	double m_bx, m_by, m_bz;
	/*depth*/
	double m_d;
	/*Poisson ratio*/
	double m_nu;

	mutable Geometry::Vector3d m_u;
};

 class DislocationAngularRight
 {
 public:
 	DislocationAngularRight(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, const Geometry::Vector3d& normal, double d, double nu);
     /** Default destructor */
     virtual ~DislocationAngularRight();
 	void moveTo(const Geometry::Vector2d& pos);
 	void reset(const Geometry::Vector3d& Line);
 	void multBurgers(int sx, int sy, int sz);

 	const Geometry::Vector3d& U(const Geometry::Vector3d& r) const;
 protected:
     void setup(const Geometry::Vector3d& Line, const Geometry::Vector3d& normal);
     /*parallel to the surface and in the plane formed by dislocation shoulders*/
 	Geometry::Vector3d m_axis_x;
 	/*parallel to the surface and perpendicular to the plane formed by dislocation shoulders*/
 	Geometry::Vector3d m_axis_y;
 	/*perpendicular to the free surface and along one of the shoulders of the dislocation*/
 	Geometry::Vector3d m_axis_z;

 	/*position on the interface*/
 	Geometry::Vector3d m_pos;

 	/*Burgers vector with respect to dislocation coord frame*/
 	double m_bx, m_by, m_bz;
 	/*depth*/
 	double m_d;
 	/*Poisson ratio*/
 	double m_nu;

 	mutable Geometry::Vector3d m_u;
 };

class DislocationPiShape
{
public:
	DislocationPiShape(const Geometry::Vector3d& Burgers,
			const Geometry::Vector2d& p1, const Geometry::Vector2d& p2,
			const Geometry::Vector3d& normal, double d, double nu);
	/** Default destructor */
	virtual ~DislocationPiShape();
	void reset(const Geometry::Vector2d& p1, const Geometry::Vector2d& p2);
	void multBurgers(int sx, int sy, int sz);

	const Geometry::Vector3d& U(const Geometry::Vector3d& r) const;
protected:

	/*pi-shape consists of two angular dislocations*/

	DislocationAngularRight * m_dislocation1, * m_dislocation2;

	mutable Geometry::Vector3d m_u;
};

#endif /* DISLOCATIONS_H_ */
