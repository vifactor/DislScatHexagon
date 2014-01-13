/*
 * Engine.h
 *
 *  Created on: 20 вер. 2013
 *      Author: kopp
 */

#ifndef ENGINE_H_
#define ENGINE_H_

#include "ProgramSettings.h"
#include "HexagonalDislocationNet.h"
#include "TriangularDislocationNet.h"
#include "ConstantStrainNet.h"

#include "LocalDisplacementCalculator.h"
#include "LocalStrainCalculator.h"
#include "MeanStrainMCCalculator.h"
#include "CoplanarIntensityMCCalculator.h"
#include "CoplanarCorrelationMCCalculator.h"
#include "MillerIndexHex.h"
#include "StringTools.h"

class Engine
{
public:
	Engine(const ProgramSettings *  settings, int seed);
	virtual ~Engine();
	virtual void run();
	virtual void setNbSteps(unsigned long long nb) {m_data->setNbSteps(nb);}
	virtual void transferDataTo(double * buff) {m_data->transferTo(buff);}
	virtual size_t getDataSize() const {return m_data->getNbPoints() * m_data->getFuncDim();}
protected:
	virtual void setupSample(const ProgramSettings::SampleSettings& stg);
	virtual void setupRng(int seed);
	virtual void setupCalculator(const ProgramSettings::EngineSettings& stg);

	virtual MCInterface * allocateStrGammaGaussInterface(const ProgramSettings::SampleSettings::StraightGammaGaussInterfaceSettings* stg);
	virtual MCInterface * allocateStrGammaInterface(const ProgramSettings::SampleSettings::StraightGammaInterfaceSettings* stg);
	virtual MCInterface * allocateStrGaussInterface(const ProgramSettings::SampleSettings::StraightGaussInterfaceSettings* stg);
	virtual MCInterface * allocateHexRandomShiftsInterface(const ProgramSettings::SampleSettings::HexRandomShiftsInterfaceSettings* stg);
	virtual MCInterface * allocateHexRandomSourcesInterface(const ProgramSettings::SampleSettings::HexRandomSourcesInterfaceSettings* stg);
	virtual MCInterface * allocateHexRandomWavesInterface(const ProgramSettings::SampleSettings::HexRandomWavesInterfaceSettings* stg);
	virtual MCInterface * allocateConstFieldInterface(const ProgramSettings::SampleSettings::ConstFieldInterfaceSettings* stg);

	virtual MCCalculator * allocateLocalDisplacementCalculator(const ProgramSettings::EngineSettings::LocalDisplacementCalculatorSettings* stg);
	virtual MCCalculator * allocateLocalStrainCalculator(const ProgramSettings::EngineSettings::LocalStrainCalculatorSettings* stg);
	virtual MCCalculator * allocateMeanStrainCalculator(const ProgramSettings::EngineSettings::MeanStrainCalculatorSettings* stg);
	virtual MCCalculator * allocateCoplanarCorrelationCalculator(const ProgramSettings::EngineSettings::CoplanarCorrelationCalculatorSettings* stg);
	virtual MCCalculator * allocateCoplanarIntensityCalculator(const ProgramSettings::EngineSettings::CoplanarIntensityCalculatorSettings* stg);

	virtual MCCalculator::MCData * allocateLocalDisplacementCalculatorData(const ProgramSettings::EngineSettings::LocalDisplacementCalculatorSettings* stg);
	virtual MCCalculator::MCData * allocateLocalStrainCalculatorData(const ProgramSettings::EngineSettings::LocalStrainCalculatorSettings* stg);
	virtual MCCalculator::MCData * allocateMeanStrainCalculatorData(const ProgramSettings::EngineSettings::MeanStrainCalculatorSettings* stg);
	virtual MCCalculator::MCData * allocateCoplanarCorrelationCalculatorData(const ProgramSettings::EngineSettings::CoplanarCorrelationCalculatorSettings* stg);
	virtual MCCalculator::MCData * allocateCoplanarIntensityCalculatorData(const ProgramSettings::EngineSettings::CoplanarIntensityCalculatorSettings* stg);
	/*gsl random number generator*/
	gsl_rng * m_rng;
	/*settings*/
	const ProgramSettings *  m_programSettings;
	/*sample*/
	MCSample * m_sample;
	/*calculator*/
	MCCalculator * m_calculator;
	/*distribution of hexagons*/
	HexagonalNet * m_hexagonalNet;
	/*distribution of straight dislocations*/
	Distribution1d * m_straightNet;
	/*interface: this variable is added in order to obey the rule that the class which allocates a pointer destroys it as well*/
	MCInterface * m_interface;
	/*MCData*/
	MCCalculator::MCData * m_data;
	/*for transformation of Miller indices*/
	MillerHexIndicesTransformator * m_miller_transformator;
};

#endif /* ENGINE_H_ */
