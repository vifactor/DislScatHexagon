/*
 * DislocationClasses.cpp
 *
 *  Created on: 1 august 2012
 *  Modified on: 22 nov. 2012
 *  Modified on: 30 sep. 2013
 *      Author: dreamcatcher
 *
 *      This source file contains implementations of function interfaces
 *      from file Dislocations.h
 */

#include "DislocationClasses.h"

using namespace Geometry;

DislocationParallel::DislocationParallel(const Geometry::Vector3d& Burgers,
		const Geometry::Vector3d& Line, const Geometry::Vector3d& normal,
		double d, double nu)
{
    m_pos = 0.0;
    m_d = d;
    m_nu = nu;

    setup(Line, normal);

	/*components of the Burgers vector*/
	m_bx = Burgers.x;
	m_by = Burgers.y;
	m_bz = Burgers.z;
}

DislocationParallel::~DislocationParallel()
{
}

const Geometry::Vector3d& DislocationParallel::U(const Geometry::Vector3d& r) const
{
    static double x, z;
    static double u, v, w;
    static Geometry::Vector3d dr;

    dr = r - m_axis_x * m_pos;

    x = dr * m_axis_x;
    //y = dr * m_axis_y;
    z = dr * m_axis_z;

    straight_dislocation_parallel_hs(x, z, m_bx, m_by, m_bz, m_nu, m_d, u, v, w);

    m_u = m_axis_x * u + m_axis_y * v + m_axis_z * w;

    return m_u;
}

void DislocationParallel::moveTo(const double pos)
{
    m_pos = pos;
}

void DislocationParallel::setup(const Geometry::Vector3d& Line, const Geometry::Vector3d& normal)
{
    /*z-axis opposite surface normal/dislocation line*/
	m_axis_z = -normal;
	m_axis_z.normalize();
	/*y-axis along dislocation line*/
	m_axis_y = Line;
	m_axis_y.normalize();
    /*x-axis perpendicular to both y and z (cross product)*/
	m_axis_x = m_axis_y%m_axis_z;
}

void DislocationParallel::multBurgers(int sx, int sy, int sz)
{
    m_bx *= sx;
    m_by *= sy;
    m_bz *= sz;
}

void DislocationParallel::rotate(const double angle)
{
	m_axis_x.RotateAboutAxis(angle, m_axis_z);
	m_axis_y.RotateAboutAxis(angle, m_axis_z);
}

DislocationAngularRight::DislocationAngularRight(const Vector3d& Burgers,
		const Vector3d& Line, const Vector3d& normal, double d, double nu)
{
	m_nu = nu;
	m_pos = Vector3d(0, 0, -d);
	m_d = d;
	setup(Line, normal);

	/*setup components of the Burgers vector
	 *
	 * x and y components are replaced on purpose, because it is more common to think about Bx component
	 * as the one perpendicular to dislocation line, whereas coordinate frame of angular dislocation
	 * is defined in such a way that component perpendicular to dislocation line is By
	 */
	m_bx = Burgers.y;
	m_by = Burgers.x;
	m_bz = Burgers.z;
}

DislocationAngularRight::~DislocationAngularRight()
{
}

const Geometry::Vector3d& DislocationAngularRight::U(const Geometry::Vector3d& r) const
{
    static double y1, y2, y3;
    static double u, v, w;
    static Geometry::Vector3d dr;

    dr = r - m_pos;

    y1 = dr * m_axis_x;
    y2 = dr * m_axis_y;
    y3 = dr * m_axis_z;

	dislocation_iso_hs_angular90(y1, y2, y3, m_bx, m_by, m_bz, m_nu, m_d, u, v,
			w);

    m_u = m_axis_x * u + m_axis_y * v + m_axis_z * w;

    return m_u;
}

void DislocationAngularRight::moveTo(const Geometry::Vector2d& pos)
{
    m_pos = Vector3d(pos, -m_d);
}

void DislocationAngularRight::setup(const Vector3d& Line, const Vector3d& normal)
{
    /*z-axis opposite surface normal/dislocation line*/
	m_axis_z = -normal;
	m_axis_z.normalize();
	/*x-axis along the shoulder of dislocation line parallel to the surface*/
	m_axis_x = Line;
	m_axis_x.normalize();
    /*y-axis perpendicular to both z and x (cross product)*/
	m_axis_y = m_axis_z%m_axis_x;
}

void DislocationAngularRight::reset(const Geometry::Vector3d& Line)
{
	/*x-axis along the shoulder of dislocation line parallel to the surface*/
	m_axis_x = Line;
	m_axis_x.normalize();
    /*y-axis perpendicular to both z and x (cross product)*/
	m_axis_y = m_axis_z%m_axis_x;
}

void DislocationAngularRight::multBurgers(int sx, int sy, int sz)
{
    m_bx *= sx;
    m_by *= sy;
    m_bz *= sz;
}

DislocationPiShape::DislocationPiShape(const Geometry::Vector3d& Burgers,
		const Geometry::Vector2d& p1, const Geometry::Vector2d& p2,
		const Geometry::Vector3d& normal, double d, double nu)
{
	Vector3d line;

	line = Vector3d((p1 - p2), 0);

	m_dislocation1 = new DislocationAngularRight(Burgers, line, normal, d, nu);
	m_dislocation1->moveTo(p1);
	m_dislocation2 = new DislocationAngularRight(-Burgers, line, normal, d, nu);
	m_dislocation2->moveTo(p2);
}

DislocationPiShape::~DislocationPiShape()
{
	delete m_dislocation1;
	delete m_dislocation2;
}

void DislocationPiShape::reset(const Geometry::Vector2d& p1, const Geometry::Vector2d& p2)
{
	Vector3d line;

	line = Vector3d((p1 - p2), 0);

	m_dislocation1->reset(line);
	m_dislocation1->moveTo(p1);
	m_dislocation2->reset(line);
	m_dislocation2->moveTo(p2);
}

const Geometry::Vector3d& DislocationPiShape::U(const Geometry::Vector3d& r) const
{
	static Vector3d u1, u2;

	u1 = m_dislocation1->U(r);
	u2 = m_dislocation2->U(r);
	m_u = u1 + u2;

	return m_u;
}
