/*
 * RootEngine.cpp
 *
 *  Created on: 21 זמגע. 2013
 *      Author: kopp
 */

#include "RootEngine.h"
using namespace Geometry;

RootEngine::RootEngine(const ProgramSettings *  settings, int seed, int nbproc, double mult) : Engine(settings, seed)
{
	m_mult = mult;
	m_nbproc = nbproc;

	copySettings(settings->getConfigfile(), settings->getEngineSettings().outfile);

	if(settings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcLDISPL)
	{
		m_cummul_data = allocateLocalDisplacementCalculatorData(settings->getEngineSettings().localDisplacementCalculatorSettings);
		m_data->setNbSteps(1);
	}else if(settings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcLSTRAIN)
	{
		m_cummul_data = allocateLocalStrainCalculatorData(settings->getEngineSettings().localStrainCalculatorSettings);
		m_data->setNbSteps(1);
	} else if(settings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcMSTRAIN)
	{
		m_cummul_data = allocateMeanStrainCalculatorData(settings->getEngineSettings().meanStrainCalculatorSettings);
		m_data->setNbSteps(1);
	} else if(settings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcCOCORRELATION)
	{
		m_cummul_data = allocateCoplanarCorrelationCalculatorData(settings->getEngineSettings().coplanarCorrelationCalculatorSettings);
		m_data->setNbSteps(1);
	} else if(settings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcCOINTENSITY)
	{
		double a0, c0;

		const MillerReciprocalHexIndices & Q = settings->getEngineSettings().coplanarIntensityCalculatorSettings->Q;
		a0 = settings->getSampleSettings().a0;
		c0 = settings->getSampleSettings().c0;

		m_cummul_data = allocateCoplanarIntensityCalculatorData(settings->getEngineSettings().coplanarIntensityCalculatorSettings);
		m_data->setNbSteps(settings->getEngineSettings().coplanarIntensityCalculatorSettings->nbsteps);
		m_intensity_scale =
				2 * M_PI * settings->getSampleSettings().thickness * settings->getSampleSettings().thickness
						* settings->getEngineSettings().coplanarIntensityCalculatorSettings->sigma_x;
						//* settings->getEngineSettings().coplanarIntensityCalculatorSettings->sigma_z;

		m_Q = 2 * M_PI * sqrt(2.0 / 3 * (Q.H * Q.H + Q.K * Q.K + Q.I * Q.I)) / a0
				+ 2 * M_PI * Q.L / c0;
	}
}

RootEngine::~RootEngine()
{
	if(m_cummul_data)
		delete m_cummul_data;
}

void RootEngine::saveData() const
{
	saveData(m_cummul_data);
}

void RootEngine::saveData(const MCCalculator::MCData * data) const
{
	if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcLDISPL)
	{
		saveLocalDisplacementData(data);
	} else if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcLSTRAIN)
	{
		saveLocalStrainData(data);
	} else if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcMSTRAIN)
	{
		saveMeanStrainData(data);
	} else if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcCOINTENSITY)
	{
		saveCoplanarIntensityData(data);
	} else if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcCOCORRELATION)
	{
		saveCoplanarCorrelationData(data);
	}
}

void RootEngine::saveLocalDisplacementData(const MCCalculator::MCData * data) const
{
	std::ofstream fout;
	double ux, uy, uz, x, y, z;

	fout.open(m_programSettings->getEngineSettings().outfile.c_str());
	fout << "x\ty\tz\tx + ux\ty + uy\tz + uz" << std::endl;
	for (size_t ipt = 0; ipt < m_data->getNbPoints(); ++ipt)
	{
		x = data->arg(ipt, LocalDisplacementCalculator::m_idx_x);
		y = data->arg(ipt, LocalDisplacementCalculator::m_idx_y);
		z = data->arg(ipt, LocalDisplacementCalculator::m_idx_z);

		ux = data->func(ipt, LocalDisplacementCalculator::m_idx_ux);
		uy = data->func(ipt, LocalDisplacementCalculator::m_idx_uy);
		uz = data->func(ipt, LocalDisplacementCalculator::m_idx_uz);

		fout << x << "\t" << y << "\t" << z << "\t" << x + ux << "\t" << y + uy
				<< "\t" << z + uz << std::endl;
	}
	fout.close();
}

void RootEngine::saveLocalStrainData(const MCCalculator::MCData * data) const
{
	std::ofstream fout;
	double spur;

	fout.open(m_programSettings->getEngineSettings().outfile.c_str());
	fout << "x\ty\tz\texx\teyy\tezz\texy\texz\teyz\tspur" << std::endl;
	for (size_t ipt = 0; ipt < m_data->getNbPoints(); ++ipt)
	{
		spur = data->func(ipt, LocalStrainCalculator::m_idx_exx)
				+ data->func(ipt, LocalStrainCalculator::m_idx_eyy)
				+ data->func(ipt, LocalStrainCalculator::m_idx_ezz);
		fout << data->arg(ipt, LocalStrainCalculator::m_idx_x) << "\t"
				<< data->arg(ipt, LocalStrainCalculator::m_idx_y) << "\t"
				<< data->arg(ipt, LocalStrainCalculator::m_idx_z) << "\t"
				<< data->func(ipt, LocalStrainCalculator::m_idx_exx) << "\t"
				<< data->func(ipt, LocalStrainCalculator::m_idx_eyy) << "\t"
				<< data->func(ipt, LocalStrainCalculator::m_idx_ezz) << "\t"
				<< data->func(ipt, LocalStrainCalculator::m_idx_exy) << "\t"
				<< data->func(ipt, LocalStrainCalculator::m_idx_exz) << "\t"
				<< data->func(ipt, LocalStrainCalculator::m_idx_eyz) << "\t" << spur
				<< std::endl;
	}
	fout.close();
}

void RootEngine::saveMeanStrainData(const MCCalculator::MCData * data) const
{
	std::ofstream fout;

	double z;
	double exx, eyy, ezz, exy, exz, eyz;
	double e2xx, e2yy, e2zz, e2xy, e2xz, e2yz;
	double Dexx, Deyy, Dezz, Dexy, Dexz, Deyz;

	fout.open(m_programSettings->getEngineSettings().outfile.c_str());
	fout << "#z\t<exx>\t<eyy>\t<ezz>\t<exy>\t<exz>\t<eyz>\tD<exx>\tD<eyy>\tD<ezz>\tD<exy>\tD<exz>\tD<eyz>" << std::endl;
	fout << "#nbSteps:\t"<< data->getNbSteps() << std::endl;
	for (size_t ipt = 0; ipt < data->getNbPoints(); ++ipt)
	{
		z = data->arg(ipt, MeanStrainMCCalculator::m_idx_z);

		exx = data->func(ipt, MeanStrainMCCalculator::m_idx_exx);
		eyy = data->func(ipt, MeanStrainMCCalculator::m_idx_eyy);
		ezz = data->func(ipt, MeanStrainMCCalculator::m_idx_ezz);
		exy = data->func(ipt, MeanStrainMCCalculator::m_idx_exy);
		exz = data->func(ipt, MeanStrainMCCalculator::m_idx_exz);
		eyz = data->func(ipt, MeanStrainMCCalculator::m_idx_eyz);

		e2xx = data->func(ipt, MeanStrainMCCalculator::m_idx_e2xx);
		e2yy = data->func(ipt, MeanStrainMCCalculator::m_idx_e2yy);
		e2zz = data->func(ipt, MeanStrainMCCalculator::m_idx_e2zz);
		e2xy = data->func(ipt, MeanStrainMCCalculator::m_idx_e2xy);
		e2xz = data->func(ipt, MeanStrainMCCalculator::m_idx_e2xz);
		e2yz = data->func(ipt, MeanStrainMCCalculator::m_idx_e2yz);

		Dexx = sqrt(e2xx - gsl_pow_2(exx));
		Deyy = sqrt(e2yy - gsl_pow_2(eyy));
		Dezz = sqrt(e2zz - gsl_pow_2(ezz));
		Dexy = sqrt(e2xy - gsl_pow_2(exy));
		Dexz = sqrt(e2xz - gsl_pow_2(exz));
		Deyz = sqrt(e2yz - gsl_pow_2(eyz));

		fout << z << "\t" << exx << "\t" << eyy << "\t" << ezz << "\t" << exy
				<< "\t" << exz << "\t" << eyz << "\t" << Dexx << "\t" << Deyy
				<< "\t" << Dezz << "\t" << Dexy << "\t" << Dexz << "\t" << Deyz
				<< "\t" << std::endl;
	}
	fout.close();
}

void RootEngine::saveCoplanarCorrelationData(const MCCalculator::MCData * data) const
{
	std::ofstream fout;

	double x, z1, z2;
	double reG, imG, phiG, r2G;
	double reT, imT;

	fout.open(m_programSettings->getEngineSettings().outfile.c_str());
	fout << "#x\tz1\tz2\treG\timG\treT\timT" << std::endl;
	fout << "#nbSteps:\t"<< data->getNbSteps() << std::endl;
	for (size_t ipt = 0; ipt < data->getNbPoints(); ++ipt)
	{
		x = data->arg(ipt, CoplanarCorrelationMCCalculator::m_idx_x);
		z1 = data->arg(ipt, CoplanarCorrelationMCCalculator::m_idx_z1);
		z2 = data->arg(ipt, CoplanarCorrelationMCCalculator::m_idx_z2);

		reG = data->func(ipt, CoplanarCorrelationMCCalculator::m_idx_reG);
		imG = data->func(ipt, CoplanarCorrelationMCCalculator::m_idx_imG);

		/* G = r_G * exp(phi_G) = exp(-T) =>
		 * T = ln(r_G) + i phi_G
		 */

		r2G = gsl_pow_2(reG) + gsl_pow_2(imG);
		phiG = atan2(imG , reG);

		reT = - 0.5 * log(r2G);
		imT = phiG;

		fout << x << "\t" << z1 << "\t" << z2 << "\t" << reG << "\t" << imG
				<< "\t" << reT << "\t" << imT << std::endl;
	}
	fout.close();
}

void RootEngine::saveCoplanarIntensityData(const MCCalculator::MCData * data) const
{
	std::ofstream fout;
	double omega;

	fout.open(m_programSettings->getEngineSettings().outfile.c_str());
	fout << "#om[deg]\tqx\tqz\treI\timI" << std::endl;
	fout << "#nbSteps:\t"<< data->getNbSteps() << std::endl;
	for (size_t ipt = 0; ipt < data->getNbPoints(); ++ipt)
	{
		omega = data->arg(ipt, CoplanarIntensityMCCalculator::m_idx_qx) / m_Q;
		fout << omega * 180.0 / M_PI << "\t" << data->arg(ipt, CoplanarIntensityMCCalculator::m_idx_qx) << "\t"
				<< data->arg(ipt, CoplanarIntensityMCCalculator::m_idx_qz) << "\t"
				<< data->func(ipt, CoplanarIntensityMCCalculator::m_idx_reI) * m_intensity_scale << "\t"
				<< data->func(ipt, CoplanarIntensityMCCalculator::m_idx_imI) * m_intensity_scale
				<< std::endl;
	}
	fout.close();
}

void RootEngine::saveHexagons(const HexagonalNet * net) const
{
	std::string filename;
	std::ofstream fout;
	EdgeIter ei, ei_end;

	filename = m_programSettings->getEngineSettings().outfile.c_str();
	stripExtension(filename);
	filename += "_hexagons.dat";

	fout.open(filename.c_str());
	for (tie(ei, ei_end) = edges(*m_hexagonalNet); ei != ei_end; ++ei)
	{
		const Vector2d& p1 = (*m_hexagonalNet)[ei->m_source].m_disturbedCoord;
		const Vector2d& p2 = (*m_hexagonalNet)[ei->m_target].m_disturbedCoord;

		fout << p1.x << "\t" << p1.y << "\t" << p2.x << "\t" << p2.y
				<< std::endl;
	}
	fout.close();
}

bool RootEngine::done()
{
	bool isDone;
	isDone = false;

	if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcLDISPL)
	{
		if(m_cummul_data->getNbSteps() > 0)
		{
			isDone = true;
			m_data->setNbSteps(0);
		}
	} else if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcLSTRAIN)
	{
		if(m_cummul_data->getNbSteps() > 0)
		{
			isDone = true;
			m_data->setNbSteps(0);
		}
	} else if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcMSTRAIN)
	{
		std::cout << m_cummul_data->getNbSteps() << std::endl;
		if(m_cummul_data->getNbSteps() >= m_programSettings->getEngineSettings().meanStrainCalculatorSettings->nbsteps)
		{
			isDone = true;
			m_data->setNbSteps(0);
		}
	} else if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcCOCORRELATION)
	{
		std::cout << m_cummul_data->getNbSteps() << std::endl;
		if(m_cummul_data->getNbSteps() >= m_programSettings->getEngineSettings().coplanarCorrelationCalculatorSettings->nbsteps)
		{
			isDone = true;
			m_data->setNbSteps(0);
		}
	} else if(m_programSettings->getEngineSettings().calculatorType == ProgramSettings::EngineSettings::calcCOINTENSITY)
	{
		double maxReI, maxImI, prec;
		maxImI = m_cummul_data->max(CoplanarIntensityMCCalculator::m_idx_imI);
		maxReI = m_cummul_data->max(CoplanarIntensityMCCalculator::m_idx_reI);
		prec = fabs(maxImI / maxReI);

		std::cout << m_cummul_data->getNbSteps() << "\t" << maxImI << "\t" << maxReI << "\t" << prec << std::endl;
		/*
		 * to finish the calculations the max imaginary intensity
		 *  must be 1 / precision times smaller
		 *   than the max real part of intensity
		 */
		if ((maxImI != 0) && (prec < m_programSettings->getEngineSettings().coplanarIntensityCalculatorSettings->precision))
		{

			isDone = true;
			m_data->setNbSteps(0);
		}

	}
	return isDone;
}

void RootEngine::appendData(const double * buff)
{
	m_cummul_data->append(buff, m_data->getNbSteps() * m_nbproc);
}

unsigned long long RootEngine::getNbSteps() const
{
	static unsigned long long nbSteps;

	nbSteps = m_data->getNbSteps();
	m_data->setNbSteps(nbSteps * m_mult);

	return nbSteps;
}

void RootEngine::copySettings(const std::string& src, const std::string& dest)
{
	std::ofstream fout;
	std::ifstream fin;
	std::string filename, line;

	filename = dest;
	stripExtension(filename);
	filename += ".~cfg";

	fout.open(filename.c_str());
	fin.open(src.c_str());

	fout << fin.rdbuf();

	fin.close();
	fout.close();
}

