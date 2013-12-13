/*
 * ProgramSettings.h
 *
 *  Created on: 11 june 2013
 *      Author: kopp
 */

#ifndef PROGRAMSETTINGS_H_
#define PROGRAMSETTINGS_H_

#include "StringTools.h"
#include "MillerIndexHex.h"
#include <libconfig.h++>

struct Range
{
	double m_min;
	double m_max;
	size_t m_sampling;

	double getStep() const
	{
		return (m_sampling > 1) ? (m_max - m_min) / (m_sampling - 1) : 0.0;
	}
	void toVector(std::vector<double>& vec) const
	{
		double step = getStep();
		for (size_t i = 0; i < m_sampling; ++i)
		{
			vec.push_back(m_min + step * i);
		}
	}
	bool good() const
	{
		if ((m_min <= m_max) && (m_sampling > 0))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

class ProgramSettings
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m)
		{
			msg = "ProgramSettings::" + m;
		}
		~Exception() throw ()
		{
		}
		const char* what() const throw ()
		{
			return msg.c_str();
		}
	private:
		std::string msg;
	};
	struct SampleSettings
	{
		double nu;
		double thickness;
		double width;
		/*hexagonal lattice parameters*/
		double a0, c0;

        struct StraightGammaGaussInterfaceSettings
		{
		    /*dislocation density*/
			double rho;
			/*order parameter*/
			int gamma;
			/*order parameter*/
			double sigma;
			/*depth at which interface is situated*/
			double depth;
			/*fraction of ideal (periodic) dislocation nets*/
			double frac;
			/*Burgers vector*/
			MillerDirectHexIndices burgers;
			/*Dislocation line*/
			MillerDirectHexIndices line;
		};

        struct StraightGammaInterfaceSettings
		{
		    /*dislocation density*/
			double rho;
			/*order parameter*/
			int gamma;
			/*depth at which interface is situated*/
			double depth;
			/*fraction of ideal (periodic) dislocation nets*/
			double frac;
			/*Burgers vector*/
			MillerDirectHexIndices burgers;
			/*Dislocation line*/
			MillerDirectHexIndices line;
		};

        struct StraightGaussInterfaceSettings
		{
		    /*dislocation density*/
			double rho;
			/*order parameter*/
			double sigma;
			/*fraction of ideal (periodic) dislocation nets*/
			double frac;
			/*depth at which interface is situated*/
			double depth;
			/*Burgers vector*/
			MillerDirectHexIndices burgers;
			/*Dislocation line*/
			MillerDirectHexIndices line;
		};

        struct HexRandomShiftsInterfaceSettings
		{
		    /*dislocation density*/
			double rho;
			/*Burgers vector*/
			MillerDirectHexIndices burgers;
			/*Dislocation line*/
			MillerDirectHexIndices line;
			/*order parameter*/
			double shift;
			/*depth at which interface is situated*/
			double depth;
		};
        struct HexRandomWavesInterfaceSettings
		{
		    /*dislocation density*/
			double rho;
			/*Burgers vector*/
			MillerDirectHexIndices burgers;
			/*Dislocation line*/
			MillerDirectHexIndices line;
			/*depth at which interface is situated*/
			double depth;
			/*distorsion wave parameter*/
			double amplitude, wavevector;
		};
        struct HexRandomSourcesInterfaceSettings
		{
		    /*dislocation density*/
			double rho;
			/*Burgers vector*/
			MillerDirectHexIndices burgers;
			/*Dislocation line*/
			MillerDirectHexIndices line;
			/*depth at which interface is situated*/
			double depth;
			/*distorsion wave parameter*/
			double alpha, frac;
			int k;
			/*how much the center of the net is moved*/
			double center_shift;
		};
		enum InterfaceType
		{
			itfUNKNOWN, itfSTRAIGHT_GG, itfSTRAIGHT_GAMMA, itfSTRAIGHT_GAUSS, itfHEXRSH, itfHEXRSO, itfHEXRW
		} interfaceType;
		StraightGammaGaussInterfaceSettings * straightGammaGaussInterfaceSettings;
        StraightGammaInterfaceSettings * straightGammaInterfaceSettings;
        StraightGaussInterfaceSettings * straightGaussInterfaceSettings;
        HexRandomShiftsInterfaceSettings * hexRandomShiftsInterfaceSettings;
        HexRandomWavesInterfaceSettings * hexRandomWavesInterfaceSettings;
        HexRandomSourcesInterfaceSettings * hexRandomSourcesInterfaceSettings;
	};
	struct EngineSettings
	{
		std::string outfile;

        struct LocalDisplacementCalculatorSettings
		{
			std::string infile;
		    /*calculation ranges*/
			Range xrange;
			Range yrange;
			Range zrange;
		};

        struct LocalStrainCalculatorSettings
		{
			double hstep;
		    /*calculation ranges*/
			Range xrange;
			Range yrange;
			Range zrange;
		};

        struct MeanStrainCalculatorSettings
		{
        	/*nb MC steps to perform*/
        	unsigned long int nbsteps;
			double hstep;
		    /*calculation ranges*/
			Range zrange;
		};

        struct CoplanarIntensityCalculatorSettings
		{
        	/*reflection*/
        	MillerReciprocalHexIndices Q;
        	/*resolutions*/
        	double sigma_x, sigma_z;
        	/*nb of initial MC steps*/
        	unsigned long int nbsteps;
		    /*calculation ranges*/
			Range qxrange, qzrange;
			double precision;
		};

        struct CoplanarCorrelationCalculatorSettings
		{
        	/*reflection*/
        	MillerReciprocalHexIndices Q;
        	/*nb MC steps to perform*/
        	unsigned long int nbsteps;
		    /*calculation ranges*/
			Range xrange, z1range, z2range;
		};

		enum CalculatorType {calcUNKNOWN, calcLDISPL, calcLSTRAIN, calcMSTRAIN,
							calcCOCORRELATION, calcCOINTENSITY} calculatorType;
		LocalDisplacementCalculatorSettings * localDisplacementCalculatorSettings;
		LocalStrainCalculatorSettings * localStrainCalculatorSettings;
		MeanStrainCalculatorSettings * meanStrainCalculatorSettings;
		CoplanarIntensityCalculatorSettings * coplanarIntensityCalculatorSettings;
		CoplanarCorrelationCalculatorSettings * coplanarCorrelationCalculatorSettings;
	};
	const SampleSettings& getSampleSettings() const
	{
		return m_sampleSettings;
	}
	const EngineSettings& getEngineSettings() const
	{
		return m_engineSettings;
	}
	const std::string& getConfigfile() const
	{
		return m_cfgfile;
	}
	ProgramSettings();
	virtual ~ProgramSettings();

	void read(const std::string& cfg);
	void print() const;
protected:
	void readSampleSettings(const libconfig::Setting& root);
	void readEngineSettings(const libconfig::Setting& root);

	void readLocalDisplacementCalculatorSettings(const libconfig::Setting& stg);
	void readLocalStrainCalculatorSettings(const libconfig::Setting& stg);
	void readMeanStrainCalculatorSettings(const libconfig::Setting& stg);
	void readCoplanarIntensityCalculatorSettings(const libconfig::Setting& stg);
	void readCoplanarCorrelationCalculatorSettings(const libconfig::Setting& stg);

	void readStrGammaGaussInterfaceSettings(const libconfig::Setting& stg);
	void readStrGammaInterfaceSettings(const libconfig::Setting& stg);
	void readStrGaussInterfaceSettings(const libconfig::Setting& stg);
	void readHexRandomShiftsInterfaceSettings(const libconfig::Setting& stg);
	void readHexRandomSourcesInterfaceSettings(const libconfig::Setting& stg);
	void readHexRandomWavesInterfaceSettings(const libconfig::Setting& stg);

	SampleSettings::InterfaceType defineInterfaceType(const libconfig::Setting& stg);
	EngineSettings::CalculatorType defineCalculatorType(const libconfig::Setting& stg);

	void printSampleSettings() const;
	void printEngineSettings() const;

	void printLocalDisplacementCalculatorSettings() const;
	void printLocalStrainCalculatorSettings() const;
	void printMeanStrainCalculatorSettings() const;
	void printCoplanarIntensityCalculatorSettings() const;
	void printCoplanarCorrelationCalculatorSettings() const;

	void printStrGammaGaussInterfaceSettings() const;
	void printStrGammaInterfaceSettings() const;
	void printStrGaussInterfaceSettings() const;
	void printHexRandomShiftsInterfaceSettings() const;
	void printHexRandomSourcesInterfaceSettings() const;
	void printHexRandomWavesInterfaceSettings() const;

	SampleSettings m_sampleSettings;
	EngineSettings m_engineSettings;
	std::string m_cfgfile;
};

#endif /* PROGRAMSETTINGS_H_ */
