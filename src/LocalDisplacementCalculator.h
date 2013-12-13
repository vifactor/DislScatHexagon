/*
 * LocalDisplacementCalculator.h
 *
 *  Created on: 11 nov. 2013
 *      Author: kopp
 */

#ifndef LOCALDISPLACEMENTCALCULATOR_H_
#define LOCALDISPLACEMENTCALCULATOR_H_

#include "MCAbstractClasses.h"

class LocalDisplacementCalculator : public MCCalculator
{
public:
	LocalDisplacementCalculator(MCSample * sample);
	virtual ~LocalDisplacementCalculator();
	virtual void run(MCData * data);
	enum {m_idx_ux, m_idx_uy, m_idx_uz, m_nb_idx_func};
	enum {m_idx_x, m_idx_y, m_idx_z, m_nb_idx_arg};
protected:
};

#endif /* LOCALDISPLACEMENTCALCULATOR_H_ */
