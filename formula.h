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


#ifndef FORMULA_H_
#define FORMULA_H_

#include "Technology.h"
#include "constant.h"
#include <math.h>

#define MAX(a,b) (((a)> (b))?(a):(b))
#define MIN(a,b) (((a)< (b))?(a):(b))

bool isPow2(int n);

/* calculate the gate capacitance */
double CalculateGateCap(double width, Technology tech);

double CalculateGateArea(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double heightTransistorRegion, Technology tech,
		double *height, double *width);

/* calculate date capacitance for FeFET cells */
double CalculateFeFETGateCap(double width, Technology tech);

/* calculate the capacitance of a gate */
void CalculateGateCapacitance(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double heightTransistorRegion, Technology tech,
		double *capInput, double *capOutput);

double CalculateDrainCap(
		double width, int type,
		double heightTransistorRegion, Technology tech);

/* calculate the capacitance of a FBRAM */
double CalculateFBRAMGateCap(double width, double thicknessFactor, Technology tech);

double CalculateFBRAMDrainCap(double width, Technology tech);

double CalculateGateLeakage(
		int gateType, int numInput,
		double widthNMOS, double widthPMOS,
		double temperature, Technology tech);

double CalculateOnResistance(double width, int type, double temperature, Technology tech);

double CalculateTransconductance(double width, int type, Technology tech);

double horowitz(double tr, double beta, double rampInput, double *rampOutput);

double CalculateWireResistance(
		double resistivity, double wireWidth, double wireThickness,
		double barrierThickness, double dishingThickness, double alphaScatter);

double CalculateWireCapacitance(
		double permittivity, double wireWidth, double wireThickness, double wireSpacing,
		double ildThickness, double millarValue, double horizontalDielectric,
		double verticalDielectic, double fringeCap);


#endif /* FORMULA_H_ */
