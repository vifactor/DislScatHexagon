#ifndef DISLOCATION_DISPLACEMENT_H_INCLUDED
#define DISLOCATION_DISPLACEMENT_H_INCLUDED

#include <gsl/gsl_math.h>

/*
    displacement due to a straight dislocation in an infinite isotropic media
    Landau & Lifshitz v. VII page 157
*/
extern void straight_dislocation_inf(double x, double y, double bx, double bz, double nu, double& u, double& v, double& w);

/*
    displacement due to a straight dislocation perpendicular to a free surface of of the half-space isotropic media
    'Elastic strain fields and dislocation mobility' edited by Indenbom, Lothe, 1992
*/
extern void straight_dislocation_perpendicular_hs(double x, double y, double z, double bx, double bz, double nu, double& u, double& v, double& w);

/*
    displacement due to a straight dislocation parallel to a free surface of the half-space isotropic media
    'Elastic strain fields and dislocation mobility' edited by Indenbom, Lothe, 1992 p.364
*/
extern void straight_dislocation_parallel_hs(double x, double z, double bx, double by, double bz, double nu, double d, double& u, double& v, double& w);

/* This function provides the displacement u(ux, uy, uz)
 * in an isotropic elastic half space at the point with coordinates (x, y, z)
 * due to slip on an angular dislocation situated at a distance 'a' from the free surface.
 * Dislocation is characterized by Burgers vector (B1, B2, B3). Angle between the dislocation shoulders is beta
 *
 *
 *	Coordinate frame:
 *		y1 - parallel to the surface and in the plane formed by dislocation shoulders
 *		y2 - parallel to the surface and perpendicular to the plane formed by dislocation shoulders
 *		y3 - perpendicular to the free surface and along one of the shoulders of the dislocation
 *
 *
 * 	Reference: Comninou and Dunders, J. of Elasticity, 5(3-4) 1975
 *
 * 		Note: Some of the equations for the B2 and B3 cases have been corrected following Thomas
 * 1993.  The equations are coded in way such that they roughly correspond
 * to each line in original text.  Exceptions have been made where it made
 * more sense because of grouping symbols.
 */
extern void dislocation_iso_hs_angular(const double& y1, const double& y2,
		const double& y3, const double& B1, const double& B2, const double& B3,
		const double& nu, const double& a, const double& beta, double& ux,
		double& uy, double& uz);

 /* This function provides the displacement u(ux, uy, uz) for for the case of angular dislocation
 * in an isotropic elastic half space with angle beta equal to 90-degrees.
 *
 * For more details, see the description of the previous function.
 *
 *
 */
extern void dislocation_iso_hs_angular90(const double& y1, const double& y2,
		const double& y3, const double& B1, const double& B2, const double& B3,
		const double& nu, const double& a, double& ux,
		double& uy, double& uz);

#endif // DISLOCATION_DISPLACEMENT_H_INCLUDED
