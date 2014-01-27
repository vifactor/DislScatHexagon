/*
 * ProgramSettings.cpp
 *
 *  Created on: 11 june 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"

Range readRange(const libconfig::Setting& stg)
{
	Range range;

	range.m_min = stg[0][0];
	range.m_max = stg[0][1];
	range.m_sampling = stg[1];

	return range;
}

std::ostream& operator<<(std::ostream& out, const Range& range)
{
	out << "[" << range.m_min << ", " << range.m_max << "]:" << range.m_sampling;
	return out;
}

MillerDirectHexIndices readMillerDirectHexIndices(const libconfig::Setting& stg)
{
    MillerDirectHexIndices index;
	if(stg.isArray() && stg.getLength() == MillerHexIndicesDimension)
	{
			index.X = stg[0];
			index.Y = stg[1];
			index.T = stg[2];
			index.Z = stg[3];
	}
	else
	{
		throw ProgramSettings::Exception("Check setting: " + toString(stg.getPath()));
	}
	return index;
}

MillerReciprocalHexIndices readMillerReciprocalHexIndices(const libconfig::Setting& stg)
{
    MillerReciprocalHexIndices index;
	if(stg.isArray() && stg.getLength() == MillerHexIndicesDimension)
	{
			index.H = stg[0];
			index.K = stg[1];
			index.I = stg[2];
			index.L = stg[3];
	}
	else
	{
		throw ProgramSettings::Exception("Check setting: " + toString(stg.getPath()));
	}
	return index;
}

std::ostream& operator<<(std::ostream& out, const MillerDirectHexIndices& index)
{
	out << "[" << index.X << ", " << index.Y << ", " << index.T << ", " << index.Z << "]";
	return out;
}

std::ostream& operator<<(std::ostream& out, const MillerReciprocalHexIndices& index)
{
	out << "[" << index.H << ", " << index.K << ", " << index.I << ", " << index.L << "]";
	return out;
}

ProgramSettings::ProgramSettings()
{
}

ProgramSettings::~ProgramSettings()
{
	if (m_sampleSettings.interfaceType == SampleSettings::itfSTRAIGHT_GAMMA)
	{
		delete m_sampleSettings.straightGammaInterfaceSettings;
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfSTRAIGHT_GG)
	{
		delete m_sampleSettings.straightGammaGaussInterfaceSettings;
	}else if (m_sampleSettings.interfaceType
			== SampleSettings::itfSTRAIGHT_GAUSS)
	{
		delete m_sampleSettings.straightGaussInterfaceSettings;
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfHEXRSH)
	{
		delete m_sampleSettings.hexRandomShiftsInterfaceSettings;
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfHEXRSO)
	{
		delete m_sampleSettings.hexRandomSourcesInterfaceSettings;
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfHEXRSO)
	{
		delete m_sampleSettings.hexRandomSourcesInterfaceSettings;
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfCONSTFIELD)
	{
		delete m_sampleSettings.constFieldInterfaceSettings;
	}

	if (m_engineSettings.calculatorType == EngineSettings::calcLDISPL)
	{
		delete m_engineSettings.localDisplacementCalculatorSettings;
	}else if (m_engineSettings.calculatorType == EngineSettings::calcLSTRAIN)
	{
		delete m_engineSettings.localStrainCalculatorSettings;
	}else if (m_engineSettings.calculatorType == EngineSettings::calcMSTRAIN)
	{
		delete m_engineSettings.meanStrainCalculatorSettings;
	}else if (m_engineSettings.calculatorType == EngineSettings::calcCOINTENSITY)
	{
		delete m_engineSettings.coplanarIntensityCalculatorSettings;
	}else if (m_engineSettings.calculatorType == EngineSettings::calcCOCORRELATION)
	{
		delete m_engineSettings.coplanarCorrelationCalculatorSettings;
	}
}

void ProgramSettings::read(const std::string& cfgfile)
{
	libconfig::Config cfg;

	m_cfgfile = cfgfile;
	// Read the file. If there is an error, report it
	try
	{
		cfg.readFile(m_cfgfile.c_str());
		cfg.setAutoConvert(true);
		const libconfig::Setting& root = cfg.getRoot();

		readSampleSettings(root);
		//readCalculatorSettings(root);
		readEngineSettings(root);
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Exception(toString(fioex.what()) + " in\t" + cfgfile);
	} catch (const libconfig::ParseException &pex)
	{
		throw Exception(
				toString(pex.what()) + " in\t" + cfgfile + ":"
						+ toString(pex.getLine()) + " - "
						+ toString(pex.getError()));
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		throw Exception(
				toString(nfex.what()) + "\t" + toString(nfex.getPath())
						+ " in\t" + cfgfile);
	} catch (libconfig::SettingTypeException& tex)
	{
		throw Exception(
				toString(tex.what()) + "\t" + toString(tex.getPath()) + " in\t"
						+ cfgfile);
	}
}

void ProgramSettings::readSampleSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &sample = root["Sample"];

	/*lattice parameters*/
	m_sampleSettings.a0 = sample["a0"];
	m_sampleSettings.c0 = sample["c0"];

	/*Poisson ratio*/
	m_sampleSettings.nu = sample["nu"];

	/*Sample sizes*/
	m_sampleSettings.thickness = sample["thickness"];
	m_sampleSettings.width = sample["width"];

	/*dislocation settings*/
	const libconfig::Setting &interface = sample["interface"];
	m_sampleSettings.interfaceType = defineInterfaceType(interface[0]);

	const libconfig::Setting &interface_settings = interface[1];
    if(m_sampleSettings.interfaceType == SampleSettings::itfSTRAIGHT_GAMMA)
	{

		m_sampleSettings.straightGammaInterfaceSettings = new SampleSettings::StraightGammaInterfaceSettings();
		readStrGammaInterfaceSettings(interface_settings);
	}else if(m_sampleSettings.interfaceType == SampleSettings::itfSTRAIGHT_GAUSS)
	{
		m_sampleSettings.straightGaussInterfaceSettings = new SampleSettings::StraightGaussInterfaceSettings();
		readStrGaussInterfaceSettings(interface_settings);
	}else if(m_sampleSettings.interfaceType == SampleSettings::itfSTRAIGHT_GG)
	{
		m_sampleSettings.straightGammaGaussInterfaceSettings = new SampleSettings::StraightGammaGaussInterfaceSettings();
		readStrGammaGaussInterfaceSettings(interface_settings);
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfHEXRSH)
	{
		m_sampleSettings.hexRandomShiftsInterfaceSettings = new SampleSettings::HexRandomShiftsInterfaceSettings();
		readHexRandomShiftsInterfaceSettings(interface_settings);
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfHEXRSO)
	{
		m_sampleSettings.hexRandomSourcesInterfaceSettings = new SampleSettings::HexRandomSourcesInterfaceSettings();
		readHexRandomSourcesInterfaceSettings(interface_settings);
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfHEXRW)
	{
		m_sampleSettings.hexRandomWavesInterfaceSettings = new SampleSettings::HexRandomWavesInterfaceSettings();
		readHexRandomWavesInterfaceSettings(interface_settings);
	}else if (m_sampleSettings.interfaceType == SampleSettings::itfCONSTFIELD)
	{
		m_sampleSettings.constFieldInterfaceSettings = new SampleSettings::ConstFieldInterfaceSettings();
		readConstFieldInterfaceSettings(interface_settings);
	} else
	{
		throw Exception("Unknown interface type");
	}
}

void ProgramSettings::readEngineSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &engine = root["Engine"];

	m_engineSettings.outfile = engine["outfile"].c_str();

	const libconfig::Setting & calculator = engine["calculator"];
	m_engineSettings.calculatorType = defineCalculatorType(calculator[0]);
	if(m_engineSettings.calculatorType == EngineSettings::calcLDISPL)
	{
		m_engineSettings.localDisplacementCalculatorSettings = new EngineSettings::LocalDisplacementCalculatorSettings();
		readLocalDisplacementCalculatorSettings(calculator[1]);
	}
	else if(m_engineSettings.calculatorType == EngineSettings::calcLSTRAIN)
	{
		m_engineSettings.localStrainCalculatorSettings = new EngineSettings::LocalStrainCalculatorSettings();
		readLocalStrainCalculatorSettings(calculator[1]);
	} else if(m_engineSettings.calculatorType == EngineSettings::calcMSTRAIN)
	{
		m_engineSettings.meanStrainCalculatorSettings = new EngineSettings::MeanStrainCalculatorSettings();
		readMeanStrainCalculatorSettings(calculator[1]);
	} else if(m_engineSettings.calculatorType == EngineSettings::calcCOINTENSITY)
	{
		m_engineSettings.coplanarIntensityCalculatorSettings = new EngineSettings::CoplanarIntensityCalculatorSettings();
		readCoplanarIntensityCalculatorSettings(calculator[1]);
	} else if(m_engineSettings.calculatorType == EngineSettings::calcCOCORRELATION)
	{
		m_engineSettings.coplanarCorrelationCalculatorSettings = new EngineSettings::CoplanarCorrelationCalculatorSettings();
		readCoplanarCorrelationCalculatorSettings(calculator[1]);
	} else
	{
		throw Exception("Unknown calculator type");
	}
}

void ProgramSettings::readLocalDisplacementCalculatorSettings(const libconfig::Setting& stg)
{
	if (stg.exists("input"))
	{
		m_engineSettings.localDisplacementCalculatorSettings->infile = stg["input"].c_str();
	}else
	{
		m_engineSettings.localDisplacementCalculatorSettings->infile = "";
		m_engineSettings.localDisplacementCalculatorSettings->xrange =
				readRange(stg["xrange"]);
		m_engineSettings.localDisplacementCalculatorSettings->yrange =
				readRange(stg["yrange"]);
		m_engineSettings.localDisplacementCalculatorSettings->zrange =
				readRange(stg["zrange"]);
	}
}

void ProgramSettings::readLocalStrainCalculatorSettings(const libconfig::Setting& stg)
{
	m_engineSettings.localStrainCalculatorSettings->hstep = stg["hstep"];
	m_engineSettings.localStrainCalculatorSettings->xrange = readRange(stg["xrange"]);
	m_engineSettings.localStrainCalculatorSettings->yrange = readRange(stg["yrange"]);
	m_engineSettings.localStrainCalculatorSettings->zrange = readRange(stg["zrange"]);
}

void ProgramSettings::readMeanStrainCalculatorSettings(const libconfig::Setting& stg)
{
	m_engineSettings.meanStrainCalculatorSettings->hstep = stg["hstep"];
	m_engineSettings.meanStrainCalculatorSettings->nbsteps = stg["nbsteps"];
	m_engineSettings.meanStrainCalculatorSettings->zrange = readRange(stg["zrange"]);
}

void ProgramSettings::readCoplanarIntensityCalculatorSettings(const libconfig::Setting& stg)
{
	m_engineSettings.coplanarIntensityCalculatorSettings->Q = readMillerReciprocalHexIndices(stg["Q"]);
	m_engineSettings.coplanarIntensityCalculatorSettings->nbsteps = stg["nbsteps"];
	m_engineSettings.coplanarIntensityCalculatorSettings->sigma_x = stg["sigma_x"];
	m_engineSettings.coplanarIntensityCalculatorSettings->sigma_z = stg["sigma_z"];
	m_engineSettings.coplanarIntensityCalculatorSettings->precision = stg["precision"];
	m_engineSettings.coplanarIntensityCalculatorSettings->qxrange = readRange(stg["qxrange"]);
	m_engineSettings.coplanarIntensityCalculatorSettings->qzrange = readRange(stg["qzrange"]);
}

void ProgramSettings::readCoplanarCorrelationCalculatorSettings(const libconfig::Setting& stg)
{
	m_engineSettings.coplanarCorrelationCalculatorSettings->Q = readMillerReciprocalHexIndices(stg["Q"]);
	m_engineSettings.coplanarCorrelationCalculatorSettings->nbsteps = stg["nbsteps"];
	m_engineSettings.coplanarCorrelationCalculatorSettings->xrange = readRange(stg["xrange"]);
	m_engineSettings.coplanarCorrelationCalculatorSettings->z1range = readRange(stg["z1range"]);
	m_engineSettings.coplanarCorrelationCalculatorSettings->z2range = readRange(stg["z2range"]);
}

ProgramSettings::SampleSettings::InterfaceType ProgramSettings::defineInterfaceType(const libconfig::Setting& stg)
{
	std::string interfaceType;

	interfaceType = stg.c_str();
	if (interfaceType.compare("STRAIGHT_GAMMA") == 0)
	{
		return SampleSettings::itfSTRAIGHT_GAMMA;
	}else if(interfaceType.compare("STRAIGHT_GG") == 0)
	{
		return SampleSettings::itfSTRAIGHT_GG;
	}else if(interfaceType.compare("STRAIGHT_GAUSS") == 0)
	{
		return SampleSettings::itfSTRAIGHT_GAUSS;
	}else if(interfaceType.compare("HEXRSH") == 0)
	{
		return SampleSettings::itfHEXRSH;
	}else if(interfaceType.compare("HEXRSO") == 0)
	{
		return SampleSettings::itfHEXRSO;
	}else if(interfaceType.compare("HEXRW") == 0)
	{
		return SampleSettings::itfHEXRW;
	}else if(interfaceType.compare("CONSTFIELD") == 0)
	{
		return SampleSettings::itfCONSTFIELD;
	}else
	{
		return SampleSettings::itfUNKNOWN;
	}
}

ProgramSettings::EngineSettings::CalculatorType ProgramSettings::defineCalculatorType(const libconfig::Setting& stg)
{
	std::string calculatorType;

	calculatorType = stg.c_str();
	if(calculatorType.compare("LOCAL_DISPLACEMENT") == 0)
	{
		return EngineSettings::calcLDISPL;
	}
	else if (calculatorType.compare("LOCAL_STRAIN") == 0)
	{
		return EngineSettings::calcLSTRAIN;
	} else if (calculatorType.compare("MEAN_STRAIN") == 0)
	{
		return EngineSettings::calcMSTRAIN;
	} else if (calculatorType.compare("COPLANAR_INTENSITY") == 0)
	{
		return EngineSettings::calcCOINTENSITY;
	}else if (calculatorType.compare("COPLANAR_CORRELATION") == 0)
	{
		return EngineSettings::calcCOCORRELATION;
	}
	else
	{
		return EngineSettings::calcUNKNOWN;
	}
}

void ProgramSettings::print() const
{
	printEngineSettings();
	printSampleSettings();
}

void ProgramSettings::printLocalDisplacementCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;

	if(!m_engineSettings.localDisplacementCalculatorSettings->infile.empty())
	{
		std::cout << "Points file:\t"
				<< m_engineSettings.localDisplacementCalculatorSettings->infile
				<< std::endl;
	}else
	{
		std::cout << "x range:\t"
				<< m_engineSettings.localDisplacementCalculatorSettings->xrange
				<< std::endl;
		std::cout << "y range:\t"
				<< m_engineSettings.localDisplacementCalculatorSettings->yrange
				<< std::endl;
		std::cout << "z range:\t"
				<< m_engineSettings.localDisplacementCalculatorSettings->zrange
				<< std::endl;
	}
}

void ProgramSettings::printLocalStrainCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;

	std::cout << "hstep:\t" << m_engineSettings.localStrainCalculatorSettings->hstep << std::endl;
	std::cout << "x range:\t" << m_engineSettings.localStrainCalculatorSettings->xrange << std::endl;
	std::cout << "y range:\t" << m_engineSettings.localStrainCalculatorSettings->yrange << std::endl;
	std::cout << "z range:\t" << m_engineSettings.localStrainCalculatorSettings->zrange << std::endl;
}

void ProgramSettings::printMeanStrainCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;

	std::cout << "hstep:\t" << m_engineSettings.meanStrainCalculatorSettings->hstep << std::endl;
	std::cout << "nb MC steps:\t" << m_engineSettings.meanStrainCalculatorSettings->nbsteps << std::endl;
	std::cout << "z range:\t" << m_engineSettings.meanStrainCalculatorSettings->zrange << std::endl;
}

void ProgramSettings::printCoplanarIntensityCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;

	std::cout << "reflection:\t" << m_engineSettings.coplanarIntensityCalculatorSettings->Q << std::endl;
	std::cout << "resolutions:\t"
			<< m_engineSettings.coplanarIntensityCalculatorSettings->sigma_x
			<< ", "
			<< m_engineSettings.coplanarIntensityCalculatorSettings->sigma_z
			<< std::endl;
	std::cout << "initial nb MC steps:\t" << m_engineSettings.coplanarIntensityCalculatorSettings->nbsteps << std::endl;
	std::cout << "required precision:\t" << m_engineSettings.coplanarIntensityCalculatorSettings->precision << std::endl;
	std::cout << "qx range:\t" << m_engineSettings.coplanarIntensityCalculatorSettings->qxrange << std::endl;
	std::cout << "qz range:\t" << m_engineSettings.coplanarIntensityCalculatorSettings->qzrange << std::endl;
}

void ProgramSettings::printCoplanarCorrelationCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;

	std::cout << "reflection:\t" << m_engineSettings.coplanarCorrelationCalculatorSettings->Q << std::endl;
	std::cout << "nb MC steps to perform:\t" << m_engineSettings.coplanarCorrelationCalculatorSettings->nbsteps << std::endl;
	std::cout << "x range:\t" << m_engineSettings.coplanarCorrelationCalculatorSettings->xrange << std::endl;
	std::cout << "z1 range:\t" << m_engineSettings.coplanarCorrelationCalculatorSettings->z1range << std::endl;
	std::cout << "z1 range:\t" << m_engineSettings.coplanarCorrelationCalculatorSettings->z2range << std::endl;
}

void ProgramSettings::printSampleSettings() const
{
	std::cout << "---Sample settings---" << std::endl;
	std::cout << "Lattice parameters: (a0, c0)\t" << m_sampleSettings.a0 << ", "
			<< m_sampleSettings.c0 << std::endl;
	std::cout << "Sample sizes (thickness width):\t" << m_sampleSettings.thickness << "\t"
			<< m_sampleSettings.width << std::endl;
	std::cout << "Poisson ratio:\t" << m_sampleSettings.nu << std::endl;

	std::cout << "Interface settings:" << std::endl;
	if(m_sampleSettings.interfaceType == SampleSettings::itfSTRAIGHT_GG)
	{
		printStrGammaGaussInterfaceSettings();
	}if(m_sampleSettings.interfaceType == SampleSettings::itfSTRAIGHT_GAMMA)
	{
		printStrGammaInterfaceSettings();
	}else if(m_sampleSettings.interfaceType == SampleSettings::itfSTRAIGHT_GAUSS)
	{
		printStrGaussInterfaceSettings();
	}else if(m_sampleSettings.interfaceType == SampleSettings::itfHEXRSH)
	{
		printHexRandomShiftsInterfaceSettings();
	}else if(m_sampleSettings.interfaceType == SampleSettings::itfHEXRSO)
	{
		printHexRandomSourcesInterfaceSettings();
	}else if(m_sampleSettings.interfaceType == SampleSettings::itfHEXRW)
	{
		printHexRandomWavesInterfaceSettings();
	}else if(m_sampleSettings.interfaceType == SampleSettings::itfCONSTFIELD)
	{
		printConstFieldInterfaceSettings();
	}
}

void ProgramSettings::printEngineSettings() const
{
	std::cout << "---Engine settings---" << std::endl;

	std::cout << "Output basename:\t" << m_engineSettings.outfile << std::endl;

	if(m_engineSettings.calculatorType == EngineSettings::calcLSTRAIN)
	{
		std::cout << "Calculation of local strain has been requested." << std::endl;
		printLocalStrainCalculatorSettings();
	} else if(m_engineSettings.calculatorType == EngineSettings::calcMSTRAIN)
	{
		std::cout << "Calculation of mean strain has been requested." << std::endl;
		printMeanStrainCalculatorSettings();
	} else if(m_engineSettings.calculatorType == EngineSettings::calcCOINTENSITY)
	{
		std::cout << "Calculation of coplanar intensity has been requested." << std::endl;
		printCoplanarIntensityCalculatorSettings();
	} else if(m_engineSettings.calculatorType == EngineSettings::calcCOCORRELATION)
	{
		std::cout << "Calculation of coplanar correlation function has been requested." << std::endl;
		printCoplanarCorrelationCalculatorSettings();
	}else if(m_engineSettings.calculatorType == EngineSettings::calcLDISPL)
	{
		std::cout << "Calculation of local displacements has been requested." << std::endl;
		printLocalDisplacementCalculatorSettings();
	}
}

void ProgramSettings::readStrGammaGaussInterfaceSettings(const libconfig::Setting& stg)
{
    m_sampleSettings.straightGammaGaussInterfaceSettings->rho = stg["rho"];
    m_sampleSettings.straightGammaGaussInterfaceSettings->gamma = stg["gamma"];
    m_sampleSettings.straightGammaGaussInterfaceSettings->sigma = stg["sigma"];
    m_sampleSettings.straightGammaGaussInterfaceSettings->frac = stg["frac"];
    m_sampleSettings.straightGammaGaussInterfaceSettings->burgers = readMillerDirectHexIndices(stg["burgers"]);
    m_sampleSettings.straightGammaGaussInterfaceSettings->line = readMillerDirectHexIndices(stg["line"]);
    m_sampleSettings.straightGammaGaussInterfaceSettings->depth = stg["depth"];
}

void ProgramSettings::readStrGammaInterfaceSettings(const libconfig::Setting& stg)
{
    m_sampleSettings.straightGammaInterfaceSettings->rho = stg["rho"];
    m_sampleSettings.straightGammaInterfaceSettings->gamma = stg["gamma"];
    m_sampleSettings.straightGammaInterfaceSettings->burgers = readMillerDirectHexIndices(stg["burgers"]);
    m_sampleSettings.straightGammaInterfaceSettings->line = readMillerDirectHexIndices(stg["line"]);
    m_sampleSettings.straightGammaInterfaceSettings->depth = stg["depth"];
}

void ProgramSettings::readStrGaussInterfaceSettings(const libconfig::Setting& stg)
{
    m_sampleSettings.straightGaussInterfaceSettings->rho = stg["rho"];
    m_sampleSettings.straightGaussInterfaceSettings->sigma = stg["sigma"];
    m_sampleSettings.straightGaussInterfaceSettings->burgers = readMillerDirectHexIndices(stg["burgers"]);
    m_sampleSettings.straightGaussInterfaceSettings->line = readMillerDirectHexIndices(stg["line"]);
    m_sampleSettings.straightGaussInterfaceSettings->depth = stg["depth"];
}

void ProgramSettings::readHexRandomShiftsInterfaceSettings(const libconfig::Setting& stg)
{
    m_sampleSettings.hexRandomShiftsInterfaceSettings->rho = stg["rho"];
    m_sampleSettings.hexRandomShiftsInterfaceSettings->shift = stg["shift"];
    m_sampleSettings.hexRandomShiftsInterfaceSettings->burgers = readMillerDirectHexIndices(stg["burgers"]);
    m_sampleSettings.hexRandomShiftsInterfaceSettings->line = readMillerDirectHexIndices(stg["line"]);
    m_sampleSettings.hexRandomShiftsInterfaceSettings->depth = stg["depth"];
}

void ProgramSettings::readHexRandomSourcesInterfaceSettings(const libconfig::Setting& stg)
{
    m_sampleSettings.hexRandomSourcesInterfaceSettings->rho = stg["rho"];
    m_sampleSettings.hexRandomSourcesInterfaceSettings->alpha = stg["alpha"];
    m_sampleSettings.hexRandomSourcesInterfaceSettings->k = stg["k"];
    m_sampleSettings.hexRandomSourcesInterfaceSettings->frac = stg["frac"];
    m_sampleSettings.hexRandomSourcesInterfaceSettings->burgers = readMillerDirectHexIndices(stg["burgers"]);
    m_sampleSettings.hexRandomSourcesInterfaceSettings->line = readMillerDirectHexIndices(stg["line"]);
    m_sampleSettings.hexRandomSourcesInterfaceSettings->depth = stg["depth"];
}

void ProgramSettings::readHexRandomWavesInterfaceSettings(const libconfig::Setting& stg)
{
    m_sampleSettings.hexRandomWavesInterfaceSettings->rho = stg["rho"];
    m_sampleSettings.hexRandomWavesInterfaceSettings->wavevector = stg["wavevector"];
    m_sampleSettings.hexRandomWavesInterfaceSettings->amplitude = stg["amplitude"];
    m_sampleSettings.hexRandomWavesInterfaceSettings->burgers = readMillerDirectHexIndices(stg["burgers"]);
    m_sampleSettings.hexRandomWavesInterfaceSettings->line = readMillerDirectHexIndices(stg["line"]);
    m_sampleSettings.hexRandomWavesInterfaceSettings->depth = stg["depth"];
}

void ProgramSettings::readConstFieldInterfaceSettings(const libconfig::Setting& stg)
{
	m_sampleSettings.constFieldInterfaceSettings->depth = stg["depth"];
    m_sampleSettings.constFieldInterfaceSettings->eps_xz = stg["eps_xz"];
    m_sampleSettings.constFieldInterfaceSettings->eps_yz = stg["eps_yz"];
    m_sampleSettings.constFieldInterfaceSettings->eps_zz = stg["eps_zz"];
}

void ProgramSettings::printStrGammaGaussInterfaceSettings() const
{
	std::cout << "\tStraight dislocation interface with gamma-gaussian distribution." << std::endl;
	std::cout << "\tDensity:\t" << m_sampleSettings.straightGammaGaussInterfaceSettings->rho
			<< std::endl;
	std::cout << "\tGamma correlation parameter:\t" << m_sampleSettings.straightGammaGaussInterfaceSettings->gamma
			<< std::endl;
	std::cout << "\tGauss correlation parameter:\t" << m_sampleSettings.straightGammaGaussInterfaceSettings->sigma
			<< std::endl;
	std::cout << "\t Fraction of ideal dislocation nets:\t" << m_sampleSettings.straightGammaGaussInterfaceSettings->frac << std::endl;
	std::cout << "\tDepth:\t" << m_sampleSettings.straightGammaGaussInterfaceSettings->depth
			<< std::endl;
	std::cout << "\tBurgers vector / Dislocation line:\t" << m_sampleSettings.straightGammaGaussInterfaceSettings->burgers << "/"
			<< m_sampleSettings.straightGammaGaussInterfaceSettings->line << std::endl;
}

void ProgramSettings::printStrGammaInterfaceSettings() const
{
	std::cout << "\tStraight dislocation interface with gamma distribution." << std::endl;
	std::cout << "\tDensity:\t" << m_sampleSettings.straightGammaInterfaceSettings->rho
			<< std::endl;
	std::cout << "\tCorrelation parameter:\t" << m_sampleSettings.straightGammaInterfaceSettings->gamma
			<< std::endl;
	std::cout << "\tDepth:\t" << m_sampleSettings.straightGammaInterfaceSettings->depth
			<< std::endl;
	std::cout << "\tBurgers vector / Dislocation line:\t" << m_sampleSettings.straightGammaInterfaceSettings->burgers << "/"
			<< m_sampleSettings.straightGammaInterfaceSettings->line << std::endl;
}

void ProgramSettings::printStrGaussInterfaceSettings() const
{
	std::cout << "\tStraight dislocation interface with Gauss shifts." << std::endl;
	std::cout << "\tDensity:\t" << m_sampleSettings.straightGaussInterfaceSettings->rho
			<< std::endl;
	std::cout << "\tCorrelation parameter:\t" << m_sampleSettings.straightGaussInterfaceSettings->sigma
			<< std::endl;
	std::cout << "\tDepth:\t" << m_sampleSettings.straightGaussInterfaceSettings->depth
			<< std::endl;
	std::cout << "\tBurgers vector / Dislocation line:\t" << m_sampleSettings.straightGaussInterfaceSettings->burgers << "/"
			<< m_sampleSettings.straightGaussInterfaceSettings->line << std::endl;
}

void ProgramSettings::printHexRandomShiftsInterfaceSettings() const
{
	std::cout << "\tHexagonal interface with random vertex shifts." << std::endl;
	std::cout << "\tDensity:\t" << m_sampleSettings.hexRandomShiftsInterfaceSettings->rho
			<< std::endl;
	std::cout << "\tCorrelation parameter:\t" << m_sampleSettings.hexRandomShiftsInterfaceSettings->shift
			<< std::endl;
	std::cout << "\tDepth:\t" << m_sampleSettings.hexRandomShiftsInterfaceSettings->depth
			<< std::endl;
	std::cout << "\tBurgers vector / Dislocation line:\t"
			<< m_sampleSettings.hexRandomShiftsInterfaceSettings->burgers << "/"
			<< m_sampleSettings.hexRandomShiftsInterfaceSettings->line
			<< std::endl;
}

void ProgramSettings::printHexRandomSourcesInterfaceSettings() const
{
	std::cout << "\tHexagonal interface with random sources of distortion." << std::endl;
	std::cout << "\tDensity:\t" << m_sampleSettings.hexRandomSourcesInterfaceSettings->rho
			<< std::endl;
	std::cout << "\tSource amplitude:\t" << m_sampleSettings.hexRandomSourcesInterfaceSettings->alpha
			<< std::endl;
	std::cout << "\tSource power:\t" << m_sampleSettings.hexRandomSourcesInterfaceSettings->k
			<< std::endl;
	std::cout << "\tFraction of sources:\t" << m_sampleSettings.hexRandomSourcesInterfaceSettings->frac
			<< std::endl;
	std::cout << "\tDepth:\t" << m_sampleSettings.hexRandomSourcesInterfaceSettings->depth
			<< std::endl;
	std::cout << "\tBurgers vector / Dislocation line:\t"
			<< m_sampleSettings.hexRandomSourcesInterfaceSettings->burgers << "/"
			<< m_sampleSettings.hexRandomSourcesInterfaceSettings->line
			<< std::endl;
}

void ProgramSettings::printHexRandomWavesInterfaceSettings() const
{
	std::cout << "\tHexagonal interface with random waves." << std::endl;
	std::cout << "\tDensity:\t" << m_sampleSettings.hexRandomWavesInterfaceSettings->rho
			<< std::endl;
	std::cout << "\tDepth:\t" << m_sampleSettings.hexRandomWavesInterfaceSettings->depth
			<< std::endl;
	std::cout << "\tBurgers vector / Dislocation line:\t"
			<< m_sampleSettings.hexRandomWavesInterfaceSettings->burgers << "/"
			<< m_sampleSettings.hexRandomWavesInterfaceSettings->line
			<< std::endl;
	std::cout << "\tAmplitude, wavevector:\t"
			<< m_sampleSettings.hexRandomWavesInterfaceSettings->amplitude
			<< ", "
			<< m_sampleSettings.hexRandomWavesInterfaceSettings->wavevector
			<< std::endl;
}

void ProgramSettings::printConstFieldInterfaceSettings() const
{
	std::cout << "\tStatic field interface." << std::endl;
	std::cout << "\tEsp_xz:\t" << m_sampleSettings.constFieldInterfaceSettings->eps_xz
			<< std::endl;
	std::cout << "\tEsp_yz:\t" << m_sampleSettings.constFieldInterfaceSettings->eps_yz
			<< std::endl;
	std::cout << "\tEsp_zz:\t" << m_sampleSettings.constFieldInterfaceSettings->eps_zz
			<< std::endl;
	std::cout << "\tDepth:\t" << m_sampleSettings.constFieldInterfaceSettings->depth
			<< std::endl;
}
