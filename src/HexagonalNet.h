/*
 * HexagonalNet.h
 *
 *  Created on: 27 вер. 2013
 *      Author: kopp
 */

#ifndef HEXAGONALNET_H_
#define HEXAGONALNET_H_

#include <boost/graph/adjacency_list.hpp>
#include "Vector2d.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>

struct VertexProperty
{
	Geometry::Vector2d m_disturbedCoord;
	Geometry::Vector2d m_undisturbedCoord;
};

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, VertexProperty> Graph;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIter;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIter;

class HexagonalNet : public Graph
{
public:
	/*dir - direction of a1 basis vector*/
	HexagonalNet(const Geometry::Vector2d & a);
	virtual ~HexagonalNet();
	virtual void update() {};
	virtual void seed(double xmin, double ymin, double xmax, double ymax);
protected:
	struct HexVertexIndex
	{
		HexVertexIndex(size_t ss = 0, size_t ii = 0, size_t jj = 0) {s = ss; i = ii; j = jj;}
		size_t s;
		int i;
		int j;
		friend bool operator<(const HexVertexIndex& lhs,
				const HexVertexIndex& rhs)
		{
			if (lhs.s < rhs.s)
				return true;
			else if ((lhs.s == rhs.s) && (lhs.i < rhs.i))
				return true;
			else if ((lhs.s == rhs.s) && (lhs.i == rhs.i) && (lhs.j < rhs.j))
				return true;
			else
				return false;
		}
	};
	std::vector<Geometry::Vector2d> m_hexagon_centers;
	typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
	typedef std::map<HexagonalNet::HexVertexIndex, vertex_t> Index2DescriptorMap;
	/*hexagonal lattice basis*/
	Geometry::Vector2d m_a1, m_a2, m_rho;
	/*distance between centers of undisturbed hexagons*/
	double m_a;
	Geometry::Vector2d getVertexCoord(const HexVertexIndex & index) const;
	Geometry::Vector2d getHexagonCoord(int i, int j) const;
	std::vector<HexVertexIndex> getVertexNeighbours(const HexVertexIndex & index) const;
	std::vector<HexVertexIndex> getHexagonVertices(int i, int j) const;
	vertex_t insertVertex(Index2DescriptorMap & mp, const HexVertexIndex & index);
};

class HexagonalNetRandomShifts : public HexagonalNet
{
public:
	HexagonalNetRandomShifts(const gsl_rng * rng, const Geometry::Vector2d& center, const Geometry::Vector2d & a, double f);
	virtual ~HexagonalNetRandomShifts() {}

	virtual void update();
private:
	const gsl_rng * m_rng;
	/*origin of hexagonal net*/
	Geometry::Vector2d m_center;
	double m_f;

	Geometry::Vector2d displacement();
};

class HexagonalNetRandomSources : public HexagonalNet
{
public:
	HexagonalNetRandomSources(const gsl_rng * rng,
			const Geometry::Vector2d& center, const Geometry::Vector2d & a,
			double alpha, int k, double frac);
	virtual ~HexagonalNetRandomSources() {}

	virtual void update();
private:
	const gsl_rng * m_rng;
	/*origin of hexagonal net*/
	Geometry::Vector2d m_center;
	struct Source
	{
		/*amplitude of distortion*/
		double m_amplitude;
		Geometry::Vector2d m_r0;
	};
	/*dispersion of amplitudes of distortions*/
	double m_alpha;
	/*power low exponent*/
	int m_k;
	/*fraction of sources*/
	double m_source_frac;
	std::vector<Source> m_sources;

	Geometry::Vector2d displacement(const Geometry::Vector2d& r);
};

class HexagonalDisturbedNetRandomWaves : public HexagonalNet
{
public:
	HexagonalDisturbedNetRandomWaves(const gsl_rng * rng,
			const Geometry::Vector2d& center, const Geometry::Vector2d& a, double amplitude,
			double wavevector);
	virtual ~HexagonalDisturbedNetRandomWaves() {}

	virtual void update();
private:
	const gsl_rng * m_rng;
	/*where the hexagons should be distributed*/
	Geometry::Vector2d m_center;
	struct Wave
	{
		/*the wave vector*/
		Geometry::Vector2d k;
		/*the vector amplitude*/
		Geometry::Vector2d u0;
	};
	/*amplitude of wave*/
	double m_displacement_ampl;
	/*wave vector amplitude*/
	double m_wavevector_ampl;

	Geometry::Vector2d m_wavevector_dir, m_displacement_dir;

	std::vector<Wave> m_waves;

	Geometry::Vector2d displacement(const Geometry::Vector2d& r);
};

#endif /* HEXAGONALNET_H_ */
