/*
 * Engine.cpp
 *
 *  Created on: 20 вер. 2013
 *      Author: kopp
 */

#include "Engine.h"

using namespace Geometry;

Engine::Engine(const ProgramSettings *  settings, int seed)
{
	m_programSettings = settings;
	m_sample = NULL;
	m_calculator = NULL;
	m_rng = NULL;
	m_interface = NULL;
	m_data = NULL;
	m_hexagonalNet = NULL;
	m_straightNet = NULL;

	m_miller_transformator = new MillerHexIndicesTransformator(
			m_programSettings->getSampleSettings().a0,
			m_programSettings->getSampleSettings().c0);

	setupRng(seed);
	setupSample(m_programSettings->getSampleSettings());
	setupCalculator(m_programSettings->getEngineSettings());

	delete m_miller_transformator;
	m_miller_transformator = NULL;
}

Engine::~Engine()
{
	if(m_interface)
		delete m_interface;
	if(m_sample)
		delete m_sample;
	if(m_calculator)
		delete m_calculator;
	if(m_hexagonalNet)
		delete m_hexagonalNet;
	if(m_straightNet)
		delete m_straightNet;

	if(m_rng)
		gsl_rng_free(m_rng);
	if(m_data)
		delete m_data;
}

void Engine::run()
{
	m_calculator->run(m_data);
}

void Engine::setupSample(const ProgramSettings::SampleSettings& stg)
{
	m_sample = new MCSample(stg.width, stg.thickness);
	MCInterface::m_lateral_size = stg.width;
	MCInterface::m_nu = stg.nu;
	MCInterface::m_normal = m_sample->normal();
	if(stg.interfaceType == ProgramSettings::SampleSettings::itfSTRAIGHT_GAMMA)
	{
		m_interface = allocateStrGammaInterface(stg.straightGammaInterfaceSettings);
	} else if(stg.interfaceType == ProgramSettings::SampleSettings::itfSTRAIGHT_GG)
	{
		m_interface = allocateStrGammaGaussInterface(stg.straightGammaGaussInterfaceSettings);
	} else if(stg.interfaceType == ProgramSettings::SampleSettings::itfSTRAIGHT_GAUSS)
	{
		m_interface = allocateStrGaussInterface(stg.straightGaussInterfaceSettings);
	}else if(stg.interfaceType == ProgramSettings::SampleSettings::itfHEXRSH)
	{
		m_interface = allocateHexRandomShiftsInterface(stg.hexRandomShiftsInterfaceSettings);
	}else if(stg.interfaceType == ProgramSettings::SampleSettings::itfHEXRSO)
	{
		m_interface = allocateHexRandomSourcesInterface(stg.hexRandomSourcesInterfaceSettings);
	}else if(stg.interfaceType == ProgramSettings::SampleSettings::itfHEXRW)
	{
		m_interface = allocateHexRandomWavesInterface(stg.hexRandomWavesInterfaceSettings);
	}else if(stg.interfaceType == ProgramSettings::SampleSettings::itfCONSTFIELD)
	{
		m_interface = allocateConstFieldInterface(stg.constFieldInterfaceSettings);
	}
	std::cout << "nb of defects:\t" << m_interface->nbDefects() << std::endl;
	m_sample->addInterface(m_interface);
}

void Engine::setupRng(int seed)
{
	m_rng = gsl_rng_alloc(gsl_rng_default);

	gsl_rng_set (m_rng, time(0) + seed);
}

void Engine::setupCalculator(const ProgramSettings::EngineSettings& stg)
{
	if(stg.calculatorType == ProgramSettings::EngineSettings::calcLDISPL)
	{
		m_calculator = allocateLocalDisplacementCalculator(stg.localDisplacementCalculatorSettings);
		m_data = allocateLocalDisplacementCalculatorData(stg.localDisplacementCalculatorSettings);
	}
	else if(stg.calculatorType == ProgramSettings::EngineSettings::calcLSTRAIN)
	{
		m_calculator = allocateLocalStrainCalculator(stg.localStrainCalculatorSettings);
		m_data = allocateLocalStrainCalculatorData(stg.localStrainCalculatorSettings);
	} else if(stg.calculatorType == ProgramSettings::EngineSettings::calcMSTRAIN)
	{
		m_calculator = allocateMeanStrainCalculator(stg.meanStrainCalculatorSettings);
		m_data = allocateMeanStrainCalculatorData(stg.meanStrainCalculatorSettings);
	} else if(stg.calculatorType == ProgramSettings::EngineSettings::calcCOINTENSITY)
	{
		m_calculator = allocateCoplanarIntensityCalculator(stg.coplanarIntensityCalculatorSettings);
		m_data = allocateCoplanarIntensityCalculatorData(stg.coplanarIntensityCalculatorSettings);
	} else if(stg.calculatorType == ProgramSettings::EngineSettings::calcCOCORRELATION)
	{
		m_calculator = allocateCoplanarCorrelationCalculator(stg.coplanarCorrelationCalculatorSettings);
		m_data = allocateCoplanarCorrelationCalculatorData(stg.coplanarCorrelationCalculatorSettings);
	}
}

MCInterface * Engine::allocateStrGammaGaussInterface(const ProgramSettings::SampleSettings::StraightGammaGaussInterfaceSettings* stg)
{
	Vector3d burgers, line;
	double bx, by, bz;

	burgers = m_miller_transformator->toVector3d(stg->burgers);
	line = m_miller_transformator->toVector3d(stg->line);
	line.normalize();

	bx = burgers * (line % MCInterface::m_normal);
	by = burgers * line;
	bz = burgers * MCInterface::m_normal;

	m_straightNet = new GammaGaussDistribution1d(m_rng, stg->gamma, stg->sigma);

	/*misfit dislocations on hexagonal interface: 1/3<11m20>/<m1100>*/
	return (new TriangularDislocationNet(m_straightNet, stg->frac,
			Vector3d(bx, by, bz),
			line, stg->rho, stg->depth));
}

MCInterface * Engine::allocateStrGammaInterface(const ProgramSettings::SampleSettings::StraightGammaInterfaceSettings* stg)
{
	Vector3d burgers, line;
	double bx, by, bz;

	burgers = m_miller_transformator->toVector3d(stg->burgers);
	line = m_miller_transformator->toVector3d(stg->line);
	line.normalize();

	bx = burgers * (line % MCInterface::m_normal);
	by = burgers * line;
	bz = burgers * MCInterface::m_normal;

	m_straightNet = new GammaDistribution1d(m_rng, stg->gamma);

	/*misfit dislocations on hexagonal interface: 1/3<11m20>/<m1100>*/
	return (new TriangularDislocationNet(m_straightNet, stg->frac,
			Vector3d(bx, by, bz),
			line, stg->rho, stg->depth));
}

MCInterface * Engine::allocateStrGaussInterface(const ProgramSettings::SampleSettings::StraightGaussInterfaceSettings* stg)
{
	Vector3d burgers, line;
	double bx, by, bz;

	burgers = m_miller_transformator->toVector3d(stg->burgers);
	line = m_miller_transformator->toVector3d(stg->line);
	line.normalize();

	bx = burgers * (line % MCInterface::m_normal);
	by = burgers * line;
	bz = burgers * MCInterface::m_normal;

	m_straightNet = new GaussDistribution1d(m_rng, stg->sigma);

	/*misfit dislocations on hexagonal interface: 1/3<11m20>/<m1100>*/
	return (new TriangularDislocationNet(m_straightNet, stg->frac,
			Vector3d(bx, by, bz),
			line, stg->rho, stg->depth));
}

MCInterface * Engine::allocateHexRandomShiftsInterface(const ProgramSettings::SampleSettings::HexRandomShiftsInterfaceSettings* stg)
{
	double a_hex;
	double w;
	Vector3d burgers, line;
	Vector2d proj;
	double bx, by, bz;

	w = MCInterface::m_lateral_size;
	/*this makes equivalence between net of straight dislocations and net of hexagonal loops with the same dislocation density*/
	a_hex = 2.0 / 3.0 / stg->rho;

	burgers = m_miller_transformator->toVector3d(stg->burgers);
	line = m_miller_transformator->toVector3d(stg->line);
	line.normalize();

	bx = burgers * (line % MCInterface::m_normal);
	by = burgers * line;
	bz = burgers * MCInterface::m_normal;

	/*since n = (0, 0, 1)*/
	proj.set(line.x, line.y);
	proj.normalize();
	proj.Rotate(M_PI / 2);
	proj *= a_hex;

	m_hexagonalNet = new HexagonalNetRandomShifts(m_rng,
			Vector2d(0, 0), proj, stg->shift);
	m_hexagonalNet->seed(-w/2, -w/2, w/2, w/2);

	std::cout << "a_hex:\t" << a_hex << std::endl;
	std::cout << "proj:\t" << proj.x << "\t" << proj.y << std::endl;

	return (new HexagonalDislocationNet(m_hexagonalNet, Vector3d(bx, by, bz), stg->depth, MCInterface::m_nu));
}

MCInterface * Engine::allocateConstFieldInterface(const ProgramSettings::SampleSettings::ConstFieldInterfaceSettings* stg)
{
	return (new ConstantStrainNet(stg->eps_xz, stg->eps_yz, stg->eps_zz, stg->depth));
}

MCInterface * Engine::allocateHexRandomSourcesInterface(const ProgramSettings::SampleSettings::HexRandomSourcesInterfaceSettings* stg)
{
	double a_hex;
	double w;
	Vector3d burgers, line;
	Vector2d proj;
	double bx, by, bz;

	w = MCInterface::m_lateral_size;
	/*this makes equivalence between net of straight dislocations and net of hexagonal loops with the same dislocation density*/
	a_hex = 2.0 / 3.0 / stg->rho;

	burgers = m_miller_transformator->toVector3d(stg->burgers);
	line = m_miller_transformator->toVector3d(stg->line);
	line.normalize();

	bx = burgers * (line % MCInterface::m_normal);
	by = burgers * line;
	bz = burgers * MCInterface::m_normal;

	/*since n = (0, 0, 1)*/
	proj.set(line.x, line.y);
	proj.normalize();
	proj.Rotate(M_PI / 2);
	proj *= a_hex;

	m_hexagonalNet = new HexagonalNetRandomSources(m_rng, proj, stg->center_shift, stg->alpha, stg->k, stg->frac);
	m_hexagonalNet->seed(-w/2, -w/2, w/2, w/2);

	std::cout << "a_hex:\t" << a_hex << std::endl;
	std::cout << "proj:\t" << proj.x << "\t" << proj.y << std::endl;

	return (new HexagonalDislocationNet(m_hexagonalNet, Vector3d(bx, by, bz), stg->depth, MCInterface::m_nu));
}

MCInterface * Engine::allocateHexRandomWavesInterface(const ProgramSettings::SampleSettings::HexRandomWavesInterfaceSettings* stg)
{
	double a_hex;
	double w;
	Vector3d burgers, line;
	Vector2d proj;
	double bx, by, bz;

	w = MCInterface::m_lateral_size;
	/*this makes equivalence between net of straight dislocations and net of hexagonal loops with the same dislocation density*/
	a_hex = 2.0 / 3.0 / stg->rho;

	burgers = m_miller_transformator->toVector3d(stg->burgers);
	line = m_miller_transformator->toVector3d(stg->line);
	line.normalize();

	bx = burgers * (line % MCInterface::m_normal);
	by = burgers * line;
	bz = burgers * MCInterface::m_normal;

	/*since n = (0, 0, 1)*/
	proj.set(line.x, line.y);
	proj.normalize();
	proj.Rotate(M_PI / 2);
	proj *= a_hex;

	m_hexagonalNet = new HexagonalDisturbedNetRandomWaves(m_rng,
			Vector2d(0, 0), proj, stg->amplitude, stg->wavevector);

	m_hexagonalNet->seed(-w/2, -w/2, w/2, w/2);

	return (new HexagonalDislocationNet(m_hexagonalNet, Vector3d(bx, by, bz), stg->depth, MCInterface::m_nu));
}

MCCalculator * Engine::allocateLocalDisplacementCalculator(const ProgramSettings::EngineSettings::LocalDisplacementCalculatorSettings* stg)
{
	return (new LocalDisplacementCalculator(m_sample));
}

MCCalculator * Engine::allocateLocalStrainCalculator(const ProgramSettings::EngineSettings::LocalStrainCalculatorSettings* stg)
{
	return (new LocalStrainCalculator(m_sample, stg->hstep));
}

MCCalculator * Engine::allocateMeanStrainCalculator(const ProgramSettings::EngineSettings::MeanStrainCalculatorSettings* stg)
{
	return (new MeanStrainMCCalculator(m_sample, m_rng, stg->hstep));
}

MCCalculator * Engine::allocateCoplanarIntensityCalculator(const ProgramSettings::EngineSettings::CoplanarIntensityCalculatorSettings * stg)
{
	return (new CoplanarIntensityMCCalculator(m_sample, m_rng, m_miller_transformator->toVector3d(stg->Q), stg->sigma_x, stg->sigma_z));
}

MCCalculator * Engine::allocateCoplanarCorrelationCalculator(const ProgramSettings::EngineSettings::CoplanarCorrelationCalculatorSettings * stg)
{
	return (new CoplanarCorrelationMCCalculator(m_sample, m_rng, m_miller_transformator->toVector3d(stg->Q)));
}

MCCalculator::MCData * Engine::allocateLocalDisplacementCalculatorData(const ProgramSettings::EngineSettings::LocalDisplacementCalculatorSettings* stg)
{
	MCCalculator::MCData * data;
	std::size_t nb_points;
	size_t ipt;
	std::vector<double> x, y, z;
	double r[3];
	std::string line;
	std::istringstream is;

	/*initialization of MCData arguments, which are x, y and z*/
	if(!stg->infile.empty())
	{
		std::ifstream fin;
		fin.open(stg->infile.c_str());
		if(!fin)
		{
			std::cout << "Cannot open the file:\t" << stg->infile << std::endl;
			return NULL;
		}
		while (!fin.eof())
		{
			getline(fin, line);

			is.str(line);
			is >> r[0] >> r[1] >> r[2];
			x.push_back(r[0]);
			y.push_back(r[1]);
			z.push_back(r[2]);
			is.clear();
		}
		fin.close();

		nb_points = x.size();
		data = new MCCalculator::MCData(nb_points, LocalDisplacementCalculator::m_nb_idx_func, LocalDisplacementCalculator::m_nb_idx_arg);
		data->init();
		for (size_t i = 0; i < nb_points; ++i)
		{
			data->arg(i, LocalDisplacementCalculator::m_idx_x) = x[i];
			data->arg(i, LocalDisplacementCalculator::m_idx_y) = y[i];
			data->arg(i, LocalDisplacementCalculator::m_idx_z) = z[i];
		}
	}else
	{
		stg->xrange.toVector(x);
		stg->yrange.toVector(y);
		stg->zrange.toVector(z);

		nb_points = stg->xrange.m_sampling * stg->yrange.m_sampling
				* stg->zrange.m_sampling;

		data = new MCCalculator::MCData(nb_points, LocalDisplacementCalculator::m_nb_idx_func, LocalDisplacementCalculator::m_nb_idx_arg);

		data->init();

		ipt = 0;
		for(size_t i = 0; i < x.size(); ++i)
		{
			for(size_t j = 0; j < y.size(); ++j)
			{
				for(size_t k = 0; k < z.size(); ++k)
				{
					data->arg(ipt, LocalDisplacementCalculator::m_idx_x) = x[i];
					data->arg(ipt, LocalDisplacementCalculator::m_idx_y) = y[j];
					data->arg(ipt, LocalDisplacementCalculator::m_idx_z) = z[k];
					++ipt;
				}
			}
		}
	}

	return data;
}

MCCalculator::MCData * Engine::allocateLocalStrainCalculatorData(const ProgramSettings::EngineSettings::LocalStrainCalculatorSettings* stg)
{
	MCCalculator::MCData * data;
	std::size_t nb_points;
	std::vector<double> x, y, z;
	size_t ipt;

	nb_points = stg->xrange.m_sampling * stg->yrange.m_sampling
			* stg->zrange.m_sampling;

	data = new MCCalculator::MCData(nb_points, LocalStrainCalculator::m_nb_idx_func, LocalStrainCalculator::m_nb_idx_arg);

	data->init();
	/*initialization of MCData arguments, which are x, y and z*/
	stg->xrange.toVector(x);
	stg->yrange.toVector(y);
	stg->zrange.toVector(z);

	ipt = 0;
	for(size_t i = 0; i < x.size(); ++i)
	{
		for(size_t j = 0; j < y.size(); ++j)
		{
			for(size_t k = 0; k < z.size(); ++k)
			{
				data->arg(ipt, LocalStrainCalculator::m_idx_x) = x[i];
				data->arg(ipt, LocalStrainCalculator::m_idx_y) = y[j];
				data->arg(ipt, LocalStrainCalculator::m_idx_z) = z[k];
				++ipt;
			}
		}
	}

	return data;
}

MCCalculator::MCData * Engine::allocateMeanStrainCalculatorData(const ProgramSettings::EngineSettings::MeanStrainCalculatorSettings* stg)
{
	MCCalculator::MCData * data;
	std::size_t nb_points;
	std::vector<double> z;
	size_t ipt;

	nb_points = stg->zrange.m_sampling;

	data = new MCCalculator::MCData(nb_points, MeanStrainMCCalculator::m_nb_idx_func, MeanStrainMCCalculator::m_nb_idx_arg);

	data->init();
	/*initialization of MCData arguments, which are x, y and z*/
	stg->zrange.toVector(z);

	ipt = 0;
	for(size_t k = 0; k < z.size(); ++k)
	{
		data->arg(ipt, MeanStrainMCCalculator::m_idx_z) = z[k];
		++ipt;
	}

	return data;
}

MCCalculator::MCData * Engine::allocateCoplanarIntensityCalculatorData(const ProgramSettings::EngineSettings::CoplanarIntensityCalculatorSettings* stg)
{
	MCCalculator::MCData * data;
	std::size_t nb_points;
	std::vector<double> qx, qz;
	size_t ipt;

	nb_points = stg->qxrange.m_sampling * stg->qzrange.m_sampling;

	data = new MCCalculator::MCData(nb_points,
			CoplanarIntensityMCCalculator::m_nb_idx_func,
			CoplanarIntensityMCCalculator::m_nb_idx_arg);

	data->init();
	/*initialization of MCData arguments, which are qx and qz*/
	stg->qxrange.toVector(qx);
	stg->qzrange.toVector(qz);

	ipt = 0;
	for(size_t i = 0; i < qx.size(); ++i)
	{
			for(size_t k = 0; k < qz.size(); ++k)
			{
				data->arg(ipt, CoplanarIntensityMCCalculator::m_idx_qx) = qx[i];
				data->arg(ipt, CoplanarIntensityMCCalculator::m_idx_qz) = qz[k];
				++ipt;
			}
	}

	return data;
}

MCCalculator::MCData * Engine::allocateCoplanarCorrelationCalculatorData(const ProgramSettings::EngineSettings::CoplanarCorrelationCalculatorSettings* stg)
{
	MCCalculator::MCData * data;
	std::size_t nb_points;
	std::vector<double> x, z1, z2;
	size_t ipt;

	nb_points = stg->xrange.m_sampling * stg->z1range.m_sampling
			* stg->z2range.m_sampling;

	data = new MCCalculator::MCData(nb_points,
			CoplanarCorrelationMCCalculator::m_nb_idx_func,
			CoplanarCorrelationMCCalculator::m_nb_idx_arg);

	data->init();
	/*initialization of MCData arguments, which are qx and qz*/
	stg->xrange.toVector(x);
	stg->z1range.toVector(z1);
	stg->z2range.toVector(z2);

	ipt = 0;
	for (size_t i = 0; i < x.size(); ++i)
	{
		for (size_t k1 = 0; k1 < z1.size(); ++k1)
		{
			for (size_t k2 = 0; k2 < z2.size(); ++k2)
			{
				data->arg(ipt, CoplanarCorrelationMCCalculator::m_idx_x) = x[i];
				data->arg(ipt, CoplanarCorrelationMCCalculator::m_idx_z1) =
						z1[k1];
				data->arg(ipt, CoplanarCorrelationMCCalculator::m_idx_z2) =
						z2[k2];
				++ipt;
			}
		}
	}
	return data;
}
