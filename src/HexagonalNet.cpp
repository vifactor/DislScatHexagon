/*
 * HexagonalNet.cpp
 *
 *  Created on: 27 вер. 2013
 *      Author: kopp
 */

#include "HexagonalNet.h"

using namespace Geometry;

HexagonalNet::HexagonalNet(const Vector2d & a)
{
	m_a = a.norm();
	m_a1 = a;
	m_a2 = a;
	/*the second basis vector is rotated by 120 degrees*/
	m_a2.Rotate(2 * M_PI / 3);

	/*rho = 2/3 a1 + 1/3 a2*/
	m_rho = m_a1 * 1.0 / 3 + m_a2 * 2.0 / 3;
}

HexagonalNet::~HexagonalNet()
{
}

void HexagonalNet::seed(double xmin, double ymin, double xmax, double ymax)
{
	int Imax, Imin, Jmax, Jmin;
	double szAlongX, szAlongY;
	Vector2d center;
	std::vector<HexagonalNet::HexVertexIndex> hexagonVertices;
	Index2DescriptorMap index2descriptorMap;
	vertex_t descriptorPrev, descriptorCur;

	/*identify which basis vector gives more contribution along (1, 0) and (0, 1 axis)*/
	if(fabs(m_a1 * Vector2d(1, 0)) >= fabs(m_a2 * Vector2d(1, 0)))
	{
		szAlongX = fabs(m_a1 * Vector2d(1, 0));
		szAlongY = fabs(m_a2 * Vector2d(0, 1));
	}
	else
	{
		szAlongX = fabs(m_a2 * Vector2d(1, 0));
		szAlongY = fabs(m_a1 * Vector2d(0, 1));
	}

	Imin = xmin / szAlongX;
	Imax = xmax / szAlongX;
	Jmin = ymin / szAlongY;
	Jmax = ymax / szAlongY;

	for (int i = Imin; i < Imax; ++i)
	{
		for (int j = Jmin; j < Jmax; ++j)
		{
			center = getHexagonCoord(i, j);
			//center of the hexagon should not be outside the predefined area
			if ((center.x < xmax) && (center.y < ymax) && (center.x > xmin)
					&& (center.y > ymin))
			{
				hexagonVertices = HexagonalNet::getHexagonVertices(i, j);

				descriptorPrev = insertVertex(index2descriptorMap, hexagonVertices[0]);
				for(size_t ivert = 1; ivert < hexagonVertices.size(); ++ivert)
				{
					descriptorCur = insertVertex(index2descriptorMap, hexagonVertices[ivert]);
					boost::add_edge(descriptorPrev, descriptorCur, *this);

					descriptorPrev = descriptorCur;
				}
				//(n-1) -> 0 edge of hexagon
				descriptorCur = insertVertex(index2descriptorMap, hexagonVertices[0]);
				boost::add_edge(descriptorPrev, descriptorCur, *this);
				m_hexagon_centers.push_back(center);
			}
		}
	}
}

HexagonalNet::vertex_t HexagonalNet::insertVertex(Index2DescriptorMap & mp, const HexVertexIndex & index)
{
	vertex_t descriptor;
	Index2DescriptorMap::const_iterator it;

	it = mp.find(index);

	if(it == mp.end())
	{
		descriptor = add_vertex(*this);
		mp[index] = descriptor;
		(*this)[descriptor].m_undisturbedCoord = getVertexCoord(index);
		(*this)[descriptor].m_disturbedCoord = getVertexCoord(index);

	} else
		descriptor = it->second;

	return descriptor;
}

Vector2d HexagonalNet::getVertexCoord(const HexVertexIndex & index) const
{
	return m_rho * index.s + m_a1 * index.i + m_a2 * index.j;
}

Vector2d HexagonalNet::getHexagonCoord(int i, int j) const
{
	/*gives the coordinates of hexagon's center*/
	return m_rho * 0.5 + m_a1 * (i + 0.5) + m_a2  * j;
}

std::vector<HexagonalNet::HexVertexIndex> HexagonalNet::getVertexNeighbours(const HexVertexIndex & index) const
{
	std::vector<HexagonalNet::HexVertexIndex> neibours;

	/*FUNCTION IS NOT USED AND PROBABLY CONTAINS MISTAKES IN THE LIST OF NEIGHBOURS*/
	/*int i, j;

	i = index.i;
	j = index.j;
	if(index.s == 0)
	{
		neibours.push_back(HexVertexIndex(1, i, j));
		neibours.push_back(HexVertexIndex(1, i, j - 1));
		neibours.push_back(HexVertexIndex(1, i - 1, j - 1));
	}else
	{
		neibours.push_back(HexVertexIndex(0, i, j));
		neibours.push_back(HexVertexIndex(0, i, j + 1));
		neibours.push_back(HexVertexIndex(0, i + 1, j + 1));
	}*/

	return neibours;
}

std::vector<HexagonalNet::HexVertexIndex> HexagonalNet::getHexagonVertices(
		int i, int j) const
{
	std::vector<HexagonalNet::HexVertexIndex> vertices;

	vertices.push_back(HexVertexIndex(0,	i,		j));
	vertices.push_back(HexVertexIndex(1,	i,		j));
	vertices.push_back(HexVertexIndex(0,	i + 1,	j + 1));
	vertices.push_back(HexVertexIndex(1,	i + 1,	j));
	vertices.push_back(HexVertexIndex(0,	i + 1,	j));
	vertices.push_back(HexVertexIndex(1,	i,	j - 1));

	return vertices;
}

HexagonalNetRandomShifts::HexagonalNetRandomShifts(const gsl_rng * rng, const Vector2d& center, const Vector2d & a, double f) :
		HexagonalNet(a)
{
	m_rng = rng;
	m_center = center;
	m_f = f;
}

void HexagonalNetRandomShifts::update()
{
	VertexIter vi, vi_end;
	static Vector2d shift;
	static double r, phi;

	/*random origin of the net*/
	//phi = 2 * M_PI * gsl_rng_uniform (m_rng);
	//r =  m_a * sqrt(gsl_rng_uniform (m_rng));
	m_center.set(0.0, 0.0);

	for (tie(vi, vi_end) = vertices(*this); vi != vi_end; ++vi)
	{
		shift = displacement();

		(*this)[*vi].m_disturbedCoord = (*this)[*vi].m_undisturbedCoord + shift + m_center;
	}
}

Vector2d HexagonalNetRandomShifts::displacement()
{
	static Vector2d res;
	static double radius, phi;

	/*uniform distribution in a circle of radius*/
	phi = 2 * M_PI * gsl_rng_uniform (m_rng);
	radius =  m_a * m_f * sqrt(gsl_rng_uniform (m_rng));

	res.x = radius * cos(phi);
	res.y = radius * sin(phi);

	return res;
}

HexagonalNetRandomSources::HexagonalNetRandomSources(
		const gsl_rng * rng, const Vector2d& a, double center_shift, double alpha, int k, double frac) :
		HexagonalNet(a)
{
	m_rng = rng;
	m_center = Vector2d(0, 0);
	m_center_shift = center_shift;
	m_alpha = alpha;
	m_source_frac = frac;
	m_k = k;
}

void HexagonalNetRandomSources::update()
{
	static VertexIter vi, vi_end;
	static Source source;
	static Vector2d r, shift;
	static double rnum;
	static double R, phi;
	/*
	 * the new distortion centers are generated here
	 * the number of sources is equal to nbhexagons * frac
	 * all source are chosen at random among hexagon centers
	 */
	m_sources.clear();

	for(size_t i = 0; i < m_hexagon_centers.size(); ++i)
	{
		rnum = gsl_rng_uniform (m_rng);
		if(rnum < m_source_frac)
		{
			source.m_r0 = m_hexagon_centers[i];
			source.m_amplitude = gsl_ran_flat(m_rng, -m_alpha, m_alpha);
			m_sources.push_back(source);
		}
	}

	/*the net origin is random*/
	phi = 2 * M_PI * gsl_rng_uniform (m_rng);
	R =  m_center_shift * m_a * sqrt(gsl_rng_uniform (m_rng));
	m_center.set(R * cos(phi), R * sin(phi));

	/*update all vertices*/
	for (tie(vi, vi_end) = vertices(*this); vi != vi_end; ++vi)
	{
		r = (*this)[*vi].m_undisturbedCoord;
		shift = displacement(r);

		(*this)[*vi].m_disturbedCoord = r + shift + m_center;
	}
}

Vector2d HexagonalNetRandomSources::displacement(const Vector2d& r)
{
	static Vector2d u, dr, drdir;
	static double drnorm;

	/*sum over all sources of displacements*/
	u = Vector2d(0, 0);
	for(size_t i = 0; i < m_sources.size(); ++i)
	{
		dr = (r - m_sources[i].m_r0);
		drnorm = dr.norm();
		drdir = dr / drnorm;
		/*every source creates the field of type alpha/r^k, where r - is a distance from the source measured in hexagon radii */
		u += drdir / gsl_pow_uint(drnorm / m_a, m_k) * m_sources[i].m_amplitude;
	}

	return u;
}

HexagonalDisturbedNetRandomWaves::HexagonalDisturbedNetRandomWaves(
		const gsl_rng * rng, const Geometry::Vector2d& center, const Vector2d& a, double amplitude,
		double wavevector) :
		HexagonalNet(a)
{
	m_rng = rng;
	m_center = center;
	m_displacement_ampl = amplitude;
	m_wavevector_ampl = wavevector;

	/*
	 * a - direction of the edge of hexagon in the hexagonal lattice
	 *
	 * we want the wave to go in the direction perpendicular a hexagon edge
	 * we want the wave to be longitudinal: amplitude || wavevector
	 */
	m_wavevector_dir = a;
	/*rotate by 30 degrees*/
	m_wavevector_dir.Rotate(M_PI / 6);
	m_wavevector_dir.normalize();
	m_displacement_dir  = m_wavevector_dir;

}

void HexagonalDisturbedNetRandomWaves::update()
{

	Wave wave;
	static double R, phi;
	static Vector2d shift, r;
	double u_ampl, k_ampl;
	Vector2d u_dir, k_dir;

	/*
	 * the new waves are generated here
	 * three direction of waves is possible
	 * amplitude and frequency are random
	 */
	m_waves.clear();
	k_dir = m_wavevector_dir;
	u_dir = m_displacement_dir;
	for(size_t i = 0; i < 3; ++i)
	{
		u_ampl = gsl_ran_gaussian(m_rng, m_displacement_ampl);
		k_ampl = m_a * (1 + gsl_ran_gaussian(m_rng, m_wavevector_ampl));

		/*plane waves extend in 3 directions 120 degrees with respect to each other*/
		k_dir.Rotate(i * 2 * M_PI / 3);
		u_dir.Rotate(i * 2 * M_PI / 3);

		wave.u0 = u_dir * u_ampl;
		wave.k = k_dir * k_ampl;
		m_waves.push_back(wave);
	}

	/*move origin(center of the hexagonal net) */
	/*the net origin is random*/
	phi = 2 * M_PI * gsl_rng_uniform (m_rng);
	R =  m_a * sqrt(gsl_rng_uniform (m_rng));
	m_center.set(R * cos(phi), R * sin(phi));

	/*update all vertices*/
	VertexIter vi, vi_end;
	for (tie(vi, vi_end) = vertices(*this); vi != vi_end; ++vi)
	{
		r = (*this)[*vi].m_undisturbedCoord;
		shift = displacement(r);

		(*this)[*vi].m_disturbedCoord = r + shift + m_center;
	}

}

Vector2d HexagonalDisturbedNetRandomWaves::displacement(const Vector2d& r)
{
	static Vector2d u, dr;

	/*sum over all sources of displacements*/
	u = Vector2d(0, 0);
	for(size_t i = 0; i < m_waves.size(); ++i)
	{
		/*every wave creates the field of type u0 cos(r*k), where r - is a coordinate of a point*/
		u += m_waves[i].u0 * cos(m_waves[i].k * r);
	}
	return u;
}

