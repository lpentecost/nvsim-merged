/*******************************************************************************
* Copyright (c) 2012-2013, The Microsystems Design Labratory (MDL)
* Department of Computer Science and Engineering, The Pennsylvania State University
* Exascale Computing Lab, Hewlett-Packard Company
* All rights reserved.
* 
* This source code is part of NVSim - An area, timing and power model for both 
* volatile (e.g., SRAM, DRAM) and non-volatile memory (e.g., PCRAM, STT-RAM, ReRAM, 
* SLC NAND Flash). The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
* 
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer.
* 
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Author list: 
*   Cong Xu	    ( Email: czx102 at psu dot edu 
*                     Website: http://www.cse.psu.edu/~czx102/ )
*   Xiangyu Dong    ( Email: xydong at cse dot psu dot edu
*                     Website: http://www.cse.psu.edu/~xydong/ )
*******************************************************************************/


#include "SenseAmp.h"
#include "formula.h"
#include "global.h"

SenseAmp::SenseAmp() {
	// TODO Auto-generated constructor stub
	initialized = false;
	invalid = false;
}

SenseAmp::~SenseAmp() {
	// TODO Auto-generated destructor stub
}

void SenseAmp::Initialize(long long _numColumn, bool _currentSense, double _senseVoltage, double _pitchSenseAmp, bool _mlc, double _numLvl, int _numF) {
	if (initialized)
		cout << "[Sense Amp] Warning: Already initialized!" << endl;

	numColumn = _numColumn;
	currentSense = _currentSense;
	senseVoltage = _senseVoltage;
	pitchSenseAmp = _pitchSenseAmp;
        mlc = _mlc;
        numLvl = _numLvl;
        numF = _numF;

	if (pitchSenseAmp <= tech->featureSize * 2) {
		/* too small, cannot do the layout */
		invalid = true;
	}

	initialized = true;
}

void SenseAmp::CalculateArea() {
	if (!initialized) {
		cout << "[Sense Amp] Error: Require initialization first!" << endl;
	} else if (invalid) {
		height = width = area = 1e41;
	} else {
		height = width = area = 0;
		double tempHeight = 0;
		double tempWidth = 0;

		CalculateGateArea(INV, 1, 0, W_SA_TOP * tech->featureSize,
				pitchSenseAmp, *tech, &tempWidth, &tempHeight);
		width = MAX(width, tempWidth);
		height += tempHeight;
		
		CalculateGateArea(INV, 1, 0, W_SA_P * tech->featureSize,
				pitchSenseAmp, *tech, &tempWidth, &tempHeight);
		width = MAX(width, tempWidth);
		height += 2 * tempHeight;
		
		CalculateGateArea(INV, 1, W_SA_N * tech->featureSize, 0,
				pitchSenseAmp, *tech, &tempWidth, &tempHeight);
		width = MAX(width, tempWidth);
		height += 2 * tempHeight;
		
		CalculateGateArea(INV, 1, W_SA_BOT * tech->featureSize, 0,
				pitchSenseAmp, *tech, &tempWidth, &tempHeight);
		width = MAX(width, tempWidth);
		height += 2 * tempHeight;
        
                /* Scale SA area for MLC programming, or defer to ISLPED modeled values for FeFET sensing */
                if (mlc && (cell->memCellType == MLCCTT || cell->memCellType == MLCRRAM)) {
                    if (numF == 8)
                        width = 2.4E-6;
                    else if (numF == 16)
                        width = 3.8E-6;
                    else if (numF == 24)
                        width = 5.2E-6;
                    else if (numF == 32)
                        width = 6.8E-6;
                    else if (numF == 40)
                        width = 8.2E-6;
                    else if (numF == 48)
                        width = 9.6E-6;
                    else if (numF == 56)
                        width = 11E-6;
            
                    height = 1.4E-6;
                    width = width * (numLvl - 1.0);
                }
		/* transformation so that width meets the pitch */
		height = height * width / pitchSenseAmp;
		width = pitchSenseAmp * numColumn;

                /* Add additional area if IV converter exists */
                height += area / width;
 
		area = height * width;

                /* Override area value for FeFET SAs */
                if (cell->memCellType == FeFET) {
                    area = numColumn*3456;
                } else if (cell->memCellType == MLCFeFET && numLvl == 4){ // 2 bits per cell
                    area = numColumn*41328; 
                } else if (cell->memCellType == MLCFeFET && numLvl == 8){ // 3 bits per cell
                    area = numColumn*113462.4; 
                } // else, given MLC config not supported yet, stick with default SA
	}
}

void SenseAmp::CalculateRC() {
	if (!initialized) {
		cout << "[Sense Amp] Error: Require initialization first!" << endl;
	} else if (invalid) {
		readLatency = writeLatency = 1e41;
	} else {
	    capLoad = CalculateGateCap((W_SA_P + W_SA_N) * tech->featureSize, *tech)
				+ CalculateDrainCap(W_SA_N * tech->featureSize, NMOS, pitchSenseAmp, *tech)
				+ CalculateDrainCap(W_SA_P * tech->featureSize, PMOS, pitchSenseAmp, *tech)
				+ CalculateDrainCap(W_SA_TOP * tech->featureSize, PMOS, pitchSenseAmp, *tech)
				+ CalculateDrainCap(W_SA_BOT * tech->featureSize, NMOS, pitchSenseAmp, *tech);
            /* Scale SA capload for MLC programming, or defer to ISLPED modeled values for FeFET sensing */
            if (mlc && (cell->memCellType == MLCCTT || cell->memCellType == MLCRRAM)) {
                if (numF == 8)
                    capLoad = 831E-18;
                else if (numF == 16)
                    capLoad = 1.092E-15;
                else if (numF == 24)
                    capLoad = 1.351E-15;
                else if (numF == 32)
                    capLoad = 1.609E-15;
                else if (numF == 40)
                    capLoad = 1.868E-15;
                else if (numF == 48)
                    capLoad = 2.125E-15;
                else if (numF == 56)
                    capLoad = 2.383E-15;
                capLoad = capLoad*(numLvl - 1.0);
            }
	}
}

void SenseAmp::CalculateLatency(double _rampInput) {	/* _rampInput is actually no use in SenseAmp */
	if (!initialized) {
		cout << "[Sense Amp] Error: Require initialization first!" << endl;
	} else {
		readLatency = writeLatency = 0;
		
		//Qing: re-model the current S/A
		if (currentSense) {
			/* all the following values achieved from HSPICE */
			if (tech->featureSize >= 119e-9)
				readLatency += 0.49e-9;		/* 120nm */
			else if (tech->featureSize >= 89e-9)
				readLatency += 0.53e-9;		/* 90nm */
			else if (tech->featureSize >= 64e-9)
				readLatency += 0.62e-9;		/* 65nm */
			else if (tech->featureSize >= 44e-9)
				readLatency += 0.80e-9;		/* 45nm */
			else if (tech->featureSize >= 31e-9)
				readLatency += 1.07e-9;		/* 32nm */
			else
			    //use new S/A number
				//readLatency += 1.45e-9;     /* below 22nm */
				readLatency += 1.0e-10;     /* below 22nm */
		}
		else {
		    /* Voltage sense amplifier */
		    double gm = CalculateTransconductance(W_SENSE_N * tech->featureSize, NMOS, *tech)
				+ CalculateTransconductance(W_SENSE_P * tech->featureSize, PMOS, *tech);
                    /* Scale gm for MLC programming, or defer to ISLPED modeled values for FeFET sensing */
                    if (mlc && (cell->memCellType == MLCCTT || cell->memCellType == MLCRRAM)) { 
                        if (numF == 8)     
                            gm = 240.3E-6;
                        else if (numF == 16)
                            gm = 340.5E-6; 
                        else if (numF == 24)
                            gm = 385.9E-6;
                        else if (numF == 32)
                            gm = 408.7E-6;
                        else if (numF == 40)
                            gm = 421.7E-6;
                        else if (numF == 48)
                            gm = 430.0E-6;
                        else if (numF == 56)
                            gm = 435.8E-6;
                    }
	            double tau = capLoad / gm;
                    /* Since this value is going to depend on the sensing voltage 
                     * we can either characterize for all sensing references or 
                     * take an average/pessimistic value: current value is vin = vdd/2 */
	            readLatency += tau * log(tech->vdd / senseVoltage);
		}

                /* Override latency value for FeFET SAs */
                if (cell->memCellType == FeFET) {
                    readLatency = 1.2e-10;
                } else if (cell->memCellType == MLCFeFET && numLvl == 4){ // 2 bits per cell
                    readLatency = 4.5e-10;
                } else if (cell->memCellType == MLCFeFET && numLvl == 8){ // 3 bits per cell
                    readLatency = 7.0e-10;
                } // else, given MLC config not supported yet, stick with default SLC SA
	}
}

void SenseAmp::CalculatePower() {
	if (!initialized) {
		cout << "[Sense Amp] Error: Require initialization first!" << endl;
	} else if (invalid) {
		readDynamicEnergy = writeDynamicEnergy = leakage = 1e41;
	} else {
		readDynamicEnergy = writeDynamicEnergy = 0;
		leakage = 0;
		
		//Qing: re-model the current S/A
		if (currentSense) {
			/* all the following values achieved from HSPICE */
			if (tech->featureSize >= 119e-9) {			/* 120nm */
				readDynamicEnergy += 8.52e-14;	/* Unit: J */
				leakage += 1.40e-8;				/* Unit: W */
			} else if (tech->featureSize >= 89e-9) {	/* 90nm */
				readDynamicEnergy += 8.72e-14;
				leakage += 1.87e-8;
			} else if (tech->featureSize >= 64e-9) {	/* 65nm */
				readDynamicEnergy += 9.00e-14;
				leakage += 2.57e-8;
			} else if (tech->featureSize >= 44e-9) {	/* 45nm */
				readDynamicEnergy += 10.26e-14;
				leakage += 4.41e-9;
			} else if (tech->featureSize >= 31e-9) {	/* 32nm */
				readDynamicEnergy += 12.56e-14;
				leakage += 12.54e-8;
			} else {    /* TO-DO, need calibration below 22nm */
			    //Qing: use new S/A numbers
				//readDynamicEnergy += 15e-14;
				readDynamicEnergy += 8e-15;
				//leakage += 15e-8;
				leakage += 0;
			}
		}
		else {
		    /* Voltage sense amplifier */
		    readDynamicEnergy += capLoad * tech->vdd * tech->vdd;
		    double idleCurrent =  CalculateGateLeakage(INV, 1, W_SENSE_EN * tech->featureSize, 0,
				inputParameter->temperature, *tech) * tech->vdd;
		    if (mlc)
                        idleCurrent = 1.35E-11 * (numLvl - 1.0);
                    leakage += idleCurrent * tech->vdd;

		}

                /* Override energy value for FeFET SAs */
                if (cell->memCellType == FeFET) {
                    readDynamicEnergy = 7.3e-16;
                } else if (cell->memCellType == MLCFeFET && numLvl == 4){ // 2 bits per cell
                    readDynamicEnergy = 1.092e-14;
                } else if (cell->memCellType == MLCFeFET && numLvl == 8){ // 3 bits per cell
                    readDynamicEnergy = 3.823e-14;
                } // else, given MLC config not supported yet, stick with default SLC SA
		
		readDynamicEnergy *= numColumn;
		leakage *= numColumn;
	}
}

void SenseAmp::PrintProperty() {
	cout << "Sense Amplifier Properties:" << endl;
	FunctionUnit::PrintProperty();
}

SenseAmp & SenseAmp::operator=(const SenseAmp &rhs) {
	height = rhs.height;
	width = rhs.width;
	area = rhs.area;
	readLatency = rhs.readLatency;
	writeLatency = rhs.writeLatency;
	readDynamicEnergy = rhs.readDynamicEnergy;
	writeDynamicEnergy = rhs.writeDynamicEnergy;
	resetLatency = rhs.resetLatency;
	setLatency = rhs.setLatency;
	resetDynamicEnergy = rhs.resetDynamicEnergy;
	setDynamicEnergy = rhs.setDynamicEnergy;
	cellReadEnergy = rhs.cellReadEnergy;
	cellSetEnergy = rhs.cellSetEnergy;
	cellResetEnergy = rhs.cellResetEnergy;
	leakage = rhs.leakage;
	initialized = rhs.initialized;
	invalid = rhs.invalid;
	numColumn = rhs.numColumn;
	currentSense = rhs.currentSense;
	senseVoltage = rhs.senseVoltage;
	capLoad = rhs.capLoad;
	pitchSenseAmp = rhs.pitchSenseAmp;

	return *this;
}
