/*
 * RootEngine.h
 *
 *  Created on: 21 זמגע. 2013
 *      Author: kopp
 */

#ifndef ROOTENGINE_H_
#define ROOTENGINE_H_

#include "Engine.h"

class RootEngine: public Engine
{
public:
	RootEngine(const ProgramSettings *  settings, int seed, int nbproc, double mult = 1.2);
	virtual ~RootEngine();
	void saveData() const;
	bool done();
	void appendData(const double * buff);
	unsigned long long getNbSteps() const;
private:
	void saveData(const MCCalculator::MCData * data) const;
	void saveHexagons(const HexagonalNet * net) const;

	void saveLocalDisplacementData(const MCCalculator::MCData * data) const;
	void saveLocalStrainData(const MCCalculator::MCData * data) const;
	void saveMeanStrainData(const MCCalculator::MCData * data) const;
	void saveCoplanarIntensityData(const MCCalculator::MCData * data) const;
	void saveCoplanarCorrelationData(const MCCalculator::MCData * data) const;

	void copySettings(const std::string& src, const std::string& dest);

	MCCalculator::MCData * m_cummul_data;
	/*multiplicator to increase the number of steps on each run*/
	double m_mult;
	int m_nbproc;
	double m_intensity_scale;
	double m_Q;
};

#endif /* ROOTENGINE_H_ */
