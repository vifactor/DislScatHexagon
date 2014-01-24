#include "CoplanarIntensityMCCalculator.h"
#include <iostream>

CoplanarIntensityMCCalculator::CoplanarIntensityMCCalculator(MCSample * sample,
		gsl_rng * rng, const Geometry::Vector3d& Q, double sigmax,
		double sigmaz) :
		MCCalculator(sample)
{
    m_rng = rng;
    m_sigma_x = sigmax;
    m_sigma_z = sigmaz;
    m_Q = Q;

    /* we need to setup the sample in child constructor,
     *  as the sample is required for laboratory frame initialization
     */
    m_sample = sample;

    /*
     * normally, the number of MC steps to perform
     * before the next configuration of the sample is generated
     * is equal to number of defects in the sample.
     *
     * for the samples with hexagonal interface the number of defects is so high,
     * that sample is updated very rarely and each output if just an average
     * over one configuration.
     *
     * to avoid this kind of issue, the maximum steps before update is set to nb_max_to_update
     */
    const size_t nb_max_to_update = 5000;
    if(m_sample->nbDefects() == 0) //when sample does not contain defects (testing purposes, resolution calculation)
    	m_steps_to_update = 1;
    else if(m_sample->nbDefects() <= nb_max_to_update)
    	m_steps_to_update = m_sample->nbDefects();
    else
    	m_steps_to_update = nb_max_to_update;

    setupLaboratoryFrame();
}

CoplanarIntensityMCCalculator::~CoplanarIntensityMCCalculator()
{
}

void CoplanarIntensityMCCalculator::setupLaboratoryFrame()
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

void CoplanarIntensityMCCalculator::add(MCData * data)
{
    static double QU, qr, arg;
    static double x, z1, z2, z, x1, x2, y1 /* == y2 */;
    static Geometry::Vector3d r1, r2, U1, U2;

    /*generate two random points*/
    z1 = -gsl_ran_flat(m_rng, 0.0, m_sample->thickness());
    z2 = -gsl_ran_flat(m_rng, 0.0, m_sample->thickness());
    z = z1 - z2;
    //do
    //{
    //    z = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_z);
    //    z2 = z1 - z;
    //} while((z2 < -m_sample->thickness()) || (z2 > 0.0));

    /*choose point randomly in the plane*/
    x1 = gsl_ran_flat(m_rng, -m_sigma_x, m_sigma_x);
    /*choose x-separation between points*/
    x = gsl_ran_gaussian_ziggurat(m_rng, m_sigma_x);
    /*since x = x1 - x2*/
    x2 = x1 - x;

    /*for coplanar geometry y1 = y2; y = 0  */
    y1 = gsl_ran_flat(m_rng, -m_sigma_x, m_sigma_x);

    /*find coordinates in the lab coord system*/
    /* minus sign before z-coordinates is due to the chosen laboratory coordinate frame*/
    r1 = m_axis_x * x1 + m_axis_z * z1 + m_axis_y * y1;
    r2 = m_axis_x * x2 + m_axis_z * z2 + m_axis_y * y1;

    U1 = m_sample->u(r1);
    U2 = m_sample->u(r2);

    QU = m_Q * (U1 - U2);

    /*std::cout << "x, z1, z2, z:\t" <<  x << "\t" << z1 << "\t" << z2 << "\t" << z << std::endl;
    std::cout << "r1:\t" << r1.x << "\t" << r1.y << "\t" << r1.z << std::endl;
    std::cout << "r2:\t" << r2.x << "\t" << r2.y << "\t" << r2.z << std::endl;
    std::cout << "U1:\t" << U1.x << "\t" << U1.y << "\t" << U1.z << std::endl;
    std::cout << "U2:\t" << U2.x << "\t" << U2.y << "\t" << U2.z << std::endl;
    std::cout << "Q:\t" << m_Q.x << "\t" << m_Q.y << "\t" << m_Q.z << std::endl;
    std::cout << "QU:\t" << QU << std::endl;*/

    for(size_t i = 0; i < data->getNbPoints(); ++i)
    {
        qr = data->arg(i, m_idx_qx) * x + data->arg(i, m_idx_qz) * z;
        arg = QU + qr;

        data->func(i, m_idx_reI) += cos(arg);
        data->func(i, m_idx_imI) += sin(arg);
    }
}

void CoplanarIntensityMCCalculator::run(MCData * data)
{
    static size_t istep;

    istep = 0;

    while(istep < data->getNbSteps())
    {
        m_sample->update();
        for(size_t idisl = 0; idisl < m_steps_to_update; ++idisl)
        {
            add(data);
            ++istep;
        }
    }
    data->setNbSteps(istep);
    //data->normalize();
}
