#include "dislocation_displacement.h"

void straight_dislocation_inf(double x, double y, double bx, double bz, double nu, double& u, double& v, double& w)
{
    const static double epsilon = 1e-6;
    double ux, vx, wx;
    double uz, vz, wz;
    double r2;

    r2 = x * x + y * y;

    ux = 0.0; vx = 0.0; wx = 0.0;
    if(fabs(bx) > epsilon)
    {
        ux = atan2(y,x) + x*y/(2 * (1.-nu) * r2);
        vx =-((1.-2*nu)/(4.*(1-nu))*log(r2)+x * x/(2.*(1.-nu)*r2));
        wx = 0.0;
    }

    uz = 0.0; vz = 0.0; wz = 0.0;
    if(fabs(bz) > epsilon)
    {
        uz = 0.0;
        vz = 0.0;
        wz = atan2(y,x);
    }

    u = (bx  * ux + bz * uz) / (2 * M_PI);
    v = (bx  * vx + bz * vz) / (2 * M_PI);
    w = (bx  * wx + bz * wz) / (2 * M_PI);
}

void straight_dislocation_perpendicular_hs(double x, double y, double z, double bx, double bz, double nu, double& u, double& v, double& w)
{
    const static double epsilon = 1e-6;
    static double ux, vx, wx;
    static double uz, vz, wz;
    static double r, r2, z2, y2, r_m_z2;

    r2 = x * x + y * y + z * z;
    r = sqrt(r2);
    r_m_z2 = gsl_pow_2(r - z);
    z2 = gsl_pow_2(z);
    y2 = gsl_pow_2(y);

    ux = 0.0; vx = 0.0; wx = 0.0;
    if(fabs(bx) > epsilon)
    {
        /**main terms correspond to Landau&Lifshitz, vol.VII, page 157 **/
        /**relaxation terms are from Lothe, in Indenbom&Lothe, page 367-368**/
        ux =//main term
            atan2(y,x) + x*y/(2 * (1.-nu) * (r2 - z2))
            //relaxation term
            +nu /(2 * (1.- nu)) * (2*x*y*z/(r * r_m_z2) + (1. - 2 * nu) * x*y/r_m_z2);
        vx =-((1.-2*nu)/(4.*(1. - nu)) * log(r2 - z2)+ x * x/(2.*(1.-nu)*(r2 - z2)))
            //relaxation term
            +nu/(2 * (1.0-nu)) * ((1. - 2 * nu) * log(r-z) - (3. - 2 * nu) * z / (r-z)
            +(3.-2*nu) * y2 / r_m_z2 - 2 * y2 / (r * (r-z)));
        wx =//relaxation term
            nu/(1.-nu) * y * (1./r + (1. - 2 * nu) * 1.0/(r-z));
    }

    uz = 0.0; vz = 0.0; wz = 0.0;
    if(fabs(bz) > epsilon)
    {
        uz = y/(r-z);
        vz = -x/(r-z);
        wz = atan2(y,x);
    }

    u = (bx  * ux + bz * uz) / (2 * M_PI);
    v = (bx  * vx + bz * vz) / (2 * M_PI);
    w = (bx  * wx + bz * wz) / (2 * M_PI);
}

void straight_dislocation_parallel_hs(double x, double z, double bx, double by, double bz, double nu, double d, double& u, double& v, double& w)
{
    const static double epsilon = 1e-6;
    static double u1xBx, u2xBx, u3xBx;
    static double u1zBx, u2zBx, u3zBx;
    static double u1xBz, u2xBz, u3xBz;
    static double u1zBz, u2zBz, u3zBz;
    static double u1yBy, u2yBy;
    static double alpha;
    static double cplus, cminus, cplus2;

    alpha = 1.0 / (2 * (1.0 - nu));
    u = 0;
    v = 0;
    w = 0;

    /**usefull variables often met in the expressions below**/
    cplus = x * x + (z + d) * (z + d);
    cminus = x * x + (z - d) * (z - d);
    cplus2 = gsl_pow_2(cplus);

    if(fabs(bx) > epsilon)
    {
        u1xBx = -(atan2(z-d, x) + alpha * x * (z - d)/cminus);
        u1zBx = ((1.0 - alpha)/2 * log(cminus)+alpha * x * x /cminus);

        u2xBx = (atan2(z+d, x) + alpha * x * (z + d)/cplus);
        u2zBx = -((1.0 - alpha)/2 * log(cplus)+alpha * x * x/cplus);

        u3xBx = 2.0 * d * ((1.0-alpha)* x /cplus-2 * alpha * x * z * (z+d)/cplus2);
        u3zBx = -2.0 * d * ((z+d)/cplus + alpha * z * ((z+d) * (z+d)-x * x)/cplus2);

        u += bx / (2 * M_PI)*(u1xBx + u2xBx + u3xBx);
        w += bx / (2 * M_PI)*(u1zBx + u2zBx + u3zBx);
    }

    if(fabs(by) > epsilon)
    {
        u1yBy =  atan2(x, z - d);
        u2yBy = -atan2(x, z + d);

        v += by/(2 * M_PI) * (u1yBy + u2yBy);
    }

    if(fabs(bz) > epsilon)
    {
        u1xBz = -((1.0 - alpha)/2 * log(cminus)+alpha * (z-d)* (z-d)/cminus);
        u1zBz = (atan2(x, z-d) + alpha * x * (z - d)/cminus);

        u2xBz = ((1.0 - alpha)/2*log(cplus)+alpha * (z+d) * (z+d)/cplus);
        u2zBz = -(atan2(x, z+d) + alpha *x * (z + d)/cplus);

        u3xBz = -2 * d *((1-alpha)*(z+d)/cplus+
                       alpha * (2 * x * x * z+d * cplus)/cplus2);
        u3zBz = -2 * d * ((1-alpha) * x/ cplus+2*alpha*x*z*(z+d)/cplus2);

        u += bz/(2 * M_PI)*(u1xBz + u2xBz + u3xBz);
        w += bz/(2 * M_PI)*(u1zBz + u2zBz + u3zBz);
    }
}

/*FIXME :
 *  for the moment to increase the computational speed of the symmetric reflections,
 *  only uz components of displacement vector is preserved
 */
void dislocation_iso_hs_angular90(const double& y1, const double& y2,
		const double& y3, const double& B1, const double& B2, const double& B3,
		const double& nu, const double& a, double& ux,
		double& uy, double& uz)
{
	static double z1, z3, z1bar, z3bar, R2, R, y3bar, R2bar, Rbar, F, Fbar,
			v1InfB1, v2InfB1, v3InfB1, v1CB1, v2CB1, v3CB1, v1B1, v2B1, v3B1,
			v1InfB2, v2InfB2, v3InfB2, v1CB2, v2CB2, v3CB2, v1B2, v2B2, v3B2,
			v1InfB3, v2InfB3, v3InfB3, v1CB3, v2CB3, v3CB3, v1B3, v2B3, v3B3;

	z1 = -y3;
	z3 = y1;
	R2 = y1 * y1 + y2 * y2 + y3 * y3;
	R = sqrt(R2);
	y3bar = y3 + 2 * a;
	z1bar = y3bar;
	z3bar = -y1;
	R2bar = y1 * y1 + y2 * y2 + y3bar * y3bar;
	Rbar = sqrt(R2bar);
	F = -atan2(y2, y1) + atan2(y2, z1) + atan2(y2 * R, y1 * z1);				//initial variant
	Fbar = -atan2(y2, y1) + atan2(y2, z1bar) + atan2(y2 * Rbar, y1 * z1bar);	//initial variant

	//F = atan2(y1, y2) - atan2(z1, y2) - atan2(y1 * z1, y2 * R);				//second variant
	//Fbar = atan2(y1, y2) - atan2(z1bar, y2) - atan2(y1 * z1bar, y2 * Rbar);	//second variant

	v1B1 = 0;
	v2B1 = 0;
	v3B1 = 0;
	if (B1 != 0)
	{
		/* Case I:Burgers vector (B1,0,0)*/
		v1InfB1 = 0.0;
		v1InfB1 = 0.0;
		/*v1InfB1 =
				2 * (1 - nu) * (F + Fbar)
						- y1 * y2
								* (1.0 / (R * (R - y3))
										+ 1.0 / (Rbar * (Rbar + y3bar)));
		v2InfB1 =
				(1 - 2 * nu) * (log(R - y3) + log(Rbar + y3bar))
						- y2 * y2
								* (1.0 / (R * (R - y3))
										+ 1.0 / (Rbar * (Rbar + y3bar)));*/
		v3InfB1 = y2 * (1.0 / R - 1.0 / Rbar);

		v1InfB1 = v1InfB1 / (8 * M_PI * (1.0 - nu));
		v2InfB1 = v2InfB1 / (8 * M_PI * (1.0 - nu));
		v3InfB1 = v3InfB1 / (8 * M_PI * (1.0 - nu));

		v1CB1 = 0.0;
		v2CB1 = 0.0;
		/*v1CB1 = (1 - 2 * nu) * y2 / (Rbar + y3bar)
				* (-y1 / (Rbar + y3bar) * (nu + a / Rbar))
				+ y2 * (y3bar - a) / (Rbar * (Rbar + y3bar))
						* (+y1 / (Rbar + y3bar) * (2 * nu + a / Rbar)
								+ a * y1 / (Rbar * Rbar));
		v2CB1 = (1 - 2 * nu) * (-nu * log(Rbar + y3bar))
				- (1 - 2 * nu) / (Rbar + y3bar)
						* (nu * y3bar - a
								+ (y2 * y2) / (Rbar + y3bar) * (nu + a / Rbar))
				+ (y3bar - a) / (Rbar + y3bar)
						* (-2 * nu + 1.0 / Rbar * (-a)
								+ (y2 * y2) / (Rbar * (Rbar + y3bar))
										* (2 * nu + a / Rbar)
								+ a * (y2 * y2) / (Rbar * Rbar * Rbar));*/
		v3CB1 = 2 * (1 - nu) * ((y2 / (Rbar + y3bar) * (2 * nu + a / Rbar)))
				+ y2 * (y3bar - a) / Rbar
						* (2 * nu / (Rbar + y3bar) + a / (Rbar * Rbar));

		v1CB1 = v1CB1 / (4 * M_PI * (1.0 - nu));
		v2CB1 = v2CB1 / (4 * M_PI * (1.0 - nu));
		v3CB1 = v3CB1 / (4 * M_PI * (1.0 - nu));

		v1B1 = v1InfB1 + v1CB1;
		v2B1 = v2InfB1 + v2CB1;
		v3B1 = v3InfB1 + v3CB1;
	}

	v1B2 = 0;
	v2B2 = 0;
	v3B2 = 0;
	if (B2 != 0)
	{
		/*Case II:Burgers vector (0,B2,0)*/

		v1InfB2 = 0.0;
		v2InfB2 = 0.0;
		/*v1InfB2 = -(1.0 - 2 * nu) * (log(R - y3) + log(Rbar + y3bar))
				+ y1 * y1
						* (1.0 / (R * (R - y3)) + 1.0 / (Rbar * (Rbar + y3bar)))
				+ z1 * (R - y1) / (R * (R - z3))
				+ z1bar * (Rbar - y1) / (Rbar * (Rbar + z3bar));
		v2InfB2 = 2 * (1.0 - nu) * (F + Fbar)
				+ y1 * y2
						* (1.0 / (R * (R - y3)) + 1.0 / (Rbar * (Rbar + y3bar)))
				- y2 * (z1 / (R * (R - z3)) + z1bar / (Rbar * (Rbar + z3bar)));*/
		v3InfB2 = -(1 - 2 * nu) * (log(R - z3) - log(Rbar + z3bar))
				- y1 * (1.0 / R - 1.0 / Rbar) + z1 * (-y3) / (R * (R - z3))
				- z1bar * y3bar / (Rbar * (Rbar + z3bar));

		v1InfB2 = v1InfB2 / (8.0 * M_PI * (1.0 - nu));
		v2InfB2 = v2InfB2 / (8.0 * M_PI * (1.0 - nu));
		v3InfB2 = v3InfB2 / (8.0 * M_PI * (1.0 - nu));

		v1CB2 = 0.0;
		v2CB2 = 0.0;
		/*v1CB2 =
				(1.0 - 2 * nu) * (nu * log(Rbar + y3bar))
						+ (1 - 2 * nu) / (Rbar + y3bar)
								* (nu * y3bar - a
										+ (y1 * y1) / (Rbar + y3bar)
												* (nu + a / Rbar))
						- (1.0 - 2 * nu) / (Rbar + z3bar)
								* (-a * (Rbar - y1) / (Rbar))
						+ (y3bar - a) / (Rbar + y3bar)
								* (2 * nu + a / Rbar
										- (y1 * y1) / (Rbar * (Rbar + y3bar))
												* (2 * nu + a / Rbar)
										- a * (y1 * y1) / (Rbar * Rbar * Rbar))
						+ (y3bar - a) / (Rbar + z3bar)
								* (a * y1 * y3bar / (Rbar * Rbar * Rbar)
										+ (Rbar - y1) / Rbar
												* (-y3bar / (Rbar + z3bar)
														* (a / Rbar)));
		v2CB2 = (1.0 - 2 * nu) * y2 / (Rbar + y3bar)
				* (y1 / (Rbar + y3bar) * (nu + a / Rbar))
				- (1.0 - 2 * nu) * y2 / (Rbar + z3bar) * (a / Rbar)
				+ y2 * (y3bar - a) / (Rbar * (Rbar + y3bar))
						* (-2 * nu * y1 / (Rbar + y3bar)
								- a * y1 / Rbar
										* (1.0 / Rbar + 1.0 / (Rbar + y3bar)))
				+ y2 * (y3bar - a) / (Rbar * (Rbar + z3bar))
						* (y3bar / (Rbar + z3bar) * (a / Rbar)
								+ a * y3bar / ((Rbar * Rbar)));*/
		v3CB2 = -2 * (1 - nu) * y1 / (Rbar + y3bar) * (2 * nu + a / Rbar)
				+ 2 * (1 - nu) * z1bar / (Rbar + z3bar) * (a / Rbar)
				+ (y3bar - a) / Rbar
						* (-2 * nu * y1 / (Rbar + y3bar)
								- a * y1 / (Rbar * Rbar))
				- (y3bar - a) / (Rbar + z3bar)
						* (a / Rbar
								* (1.0 - y3bar * z1bar / (Rbar * Rbar)
										- z1bar * y3bar
												/ (Rbar * (Rbar + z3bar))));

		v1CB2 = v1CB2 / (4 * M_PI * (1 - nu));
		v2CB2 = v2CB2 / (4 * M_PI * (1 - nu));
		v3CB2 = v3CB2 / (4 * M_PI * (1 - nu));

		v1B2 = v1InfB2 + v1CB2;
		v2B2 = v2InfB2 + v2CB2;
		v3B2 = v3InfB2 + v3CB2;
	}

	v1B3 = 0;
	v2B3 = 0;
	v3B3 = 0;
	if (B3 != 0)
	{
		v1InfB3 = 0.0;
		v2InfB3 = 0.0;
		/*Case III:Burgers vector (0,0,B3)*/
		/*v1InfB3 = y2
				* ((R - y1) / (R * (R - z3))
						+ (Rbar - y1) / (Rbar * (Rbar + z3bar)));
		v2InfB3 =
				(1.0 - 2 * nu) * (log(R - z3) + log(Rbar + z3bar))
						- (y2 * y2)
								* (1.0 / (R * (R - z3))
										+ 1.0 / (Rbar * (Rbar + z3bar)));*/
		v3InfB3 = 2 * (1.0 - nu) * (F - Fbar)
				+ y2 * (-y3 / (R * (R - z3)) - y3bar / (Rbar * (Rbar + z3bar)));

		v1InfB3 = v1InfB3 / (8 * M_PI * (1 - nu));
		v2InfB3 = v2InfB3 / (8 * M_PI * (1 - nu));
		v3InfB3 = v3InfB3 / (8 * M_PI * (1 - nu));

		v1CB3 = 0.0;
		v2CB3 = 0.0;
		/*v1CB3 = (1.0 - 2 * nu) * (y2 / (Rbar + y3bar) * (1 + a / Rbar))
				- y2 * (y3bar - a) / Rbar
						* (a / (Rbar * Rbar) + 1.0 / (Rbar + y3bar));
		v2CB3 = (1.0 - 2 * nu)
				* (-log(Rbar + z3bar) - y1 / (Rbar + y3bar) * (1 + a / Rbar)
						+ z1bar / (Rbar + z3bar) * (a / Rbar))
				+ y1 * (y3bar - a) / Rbar
						* (a / (Rbar * Rbar) + 1.0 / (Rbar + y3bar))
				- (y3bar - a) / (Rbar + z3bar)
						* ((-a / Rbar)
								+ z1bar / Rbar * (1 + a * y3bar / (Rbar * Rbar))
								- 1.0 / (Rbar * (Rbar + z3bar))
										* (-a * z1bar / Rbar * y3bar));*/
		v3CB3 = 2 * (1 - nu) * Fbar
				+ 2 * (1 - nu) * (y2 / (Rbar + z3bar) * (a / Rbar))
				+ y2 * (y3bar - a) / (Rbar * (Rbar + z3bar))
						* (1 + y3bar / (Rbar + z3bar) * (a / Rbar)
								+ a * y3bar / (Rbar * Rbar));

		v1CB3 = v1CB3 / (4 * M_PI * (1.0 - nu));
		v2CB3 = v2CB3 / (4 * M_PI * (1.0 - nu));
		v3CB3 = v3CB3 / (4 * M_PI * (1.0 - nu));

		v1B3 = v1InfB3 + v1CB3;
		v2B3 = v2InfB3 + v2CB3;
		v3B3 = v3InfB3 + v3CB3;
	}

	/*Sum the for each slip component */
	ux = B1 * v1B1 + B2 * v1B2 + B3 * v1B3;
	uy = B1 * v2B1 + B2 * v2B2 + B3 * v2B3;
	uz = B1 * v3B1 + B2 * v3B2 + B3 * v3B3;
}
