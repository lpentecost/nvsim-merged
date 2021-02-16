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
* Contributed by: 
*   Qing Guo	    ( qguo@cs.rochester.edu )
* Note:
*   To support subarray differential write
*******************************************************************************/


#include "Buffer.h"
#include "formula.h"
#include "global.h"

Buffer::Buffer() {
	// TODO Auto-generated constructor stub
	initialized = false;
	invalid = false;
}

Buffer::~Buffer() {
	// TODO Auto-generated destructor stub
}

void Buffer::Initialize(long long _numRow, long long _numColumn) {
	if (initialized)
		cout << "[Buffer] Warning: Already initialized!" << endl;

	numRow = _numRow;
	numColumn = _numColumn;
	
	initialized = true;
}

void Buffer::CalculateArea() {
    double cellheight, cellwidth;
    double xorheight, xorwidth;
    
	if (!initialized) {
		cout << "[Buffer] Error: Require initialization first!" << endl;
	} else if (invalid) {
		height = width = area = 1e41;
	} else {
		if (tech->featureSize >= 44e-9) {       /* 45nm */
		    cellheight = 1.52e-6;
		    cellwidth = 1.54375e-6;
		    xorheight = 2.166e-6;
		    xorwidth = 2.166e-6;
		}
		else if (tech->featureSize >= 31e-9) {  /* 32nm */
		    /* scaled from 45nm */
			cellheight = 1.081e-6;
			cellwidth = 1.098e-6;
			xorheight = 1.54e-6;
			xorwidth = 1.54e-6;
		}
		else {                                  /* below 22nm */
		    /* scaled from 45nm */
			cellheight = 0.743e-6;
			cellwidth = 0.755e-6;
			xorheight = 1.059e-6;
			xorwidth = 1.059e-6;
		}
		
		height = (cellheight + xorheight) * numRow;
		//numColumn = num of data bits + num of control signals. Control signals
		//do not need XOR, allowing more room to place XORs.
		width = cellwidth * numColumn;

		area = height * width;
	}
}

void Buffer::CalculateRC() {
	if (!initialized) {
		cout << "[Buffer] Error: Require initialization first!" << endl;
	} else if (invalid) {
		readLatency = writeLatency = xorLatency = 1e41;
	} else {
	    /* nothing to do */
	}
}

void Buffer::CalculateLatency() {
	if (!initialized) {
		cout << "[Buffer] Error: Require initialization first!" << endl;
	} else {
		if (tech->featureSize >= 44e-9) {       /* 45nm */
		    readLatency = 0.043e-9;
		    writeLatency = readLatency;
		    xorLatency = 0.381e-9;
		}
		else if (tech->featureSize >= 31e-9) {  /* 32nm */
		    /* copied from 45nm, need to figure out */
		    readLatency = 0.043e-9;
		    writeLatency = readLatency;
		    xorLatency = 0.381e-9;
		}
		else {                                  /* below 22nm */
		    /* scaled from 45nm */
		    readLatency = 0.024e-9;
		    writeLatency = readLatency;
		    xorLatency = 0.21e-9;
		}
	}
	
	setLatency = resetLatency = writeLatency;
}

void Buffer::CalculatePower() {
	if (!initialized) {
		cout << "[Buffer] Error: Require initialization first!" << endl;
	} else if (invalid) {
		readDynamicEnergy = writeDynamicEnergy = leakage = 1e41;
		xorDynamicEnergy = xorLeakage = 1e41;
	} else {
		if (tech->featureSize >= 44e-9) {       /* 45nm */
		    readDynamicEnergy = 0.0002e-18;
		    writeDynamicEnergy = readDynamicEnergy;
		    leakage = 19.7536e-9;
		    xorDynamicEnergy = 0.0034e-18;
		    xorLeakage = 45.1492e-9;
		}
		else if (tech->featureSize >= 31e-9) {  /* 32nm */
		    /* copied from 45nm, need to figure out */
		    readDynamicEnergy = 0.0002e-18;
		    writeDynamicEnergy = readDynamicEnergy;
		    leakage = 19.7536e-9;
		    xorDynamicEnergy = 0.0034e-18;
		    xorLeakage = 45.1492e-9;
		}
		else {                                  /* below 22nm */
		    /* scaled from 45nm */
		    readDynamicEnergy = 0.00006e-18;
		    writeDynamicEnergy = readDynamicEnergy;
		    leakage = 45.796e-9;
		    xorDynamicEnergy = 0.001e-18;
		    xorLeakage = 104.67e-9;
		}
		
		cellReadEnergy = readDynamicEnergy;
		cellSetEnergy = cellResetEnergy = writeDynamicEnergy;
		
		readDynamicEnergy *= numColumn;
		writeDynamicEnergy = readDynamicEnergy;
		leakage *= numColumn;
		xorDynamicEnergy *= numColumn;
		xorLeakage *= numColumn;
	}
	
	setDynamicEnergy = resetDynamicEnergy = writeDynamicEnergy;
}

void Buffer::PrintProperty() {
	cout << "Subarray Buffer Properties:" << endl;
	FunctionUnit::PrintProperty();
	cout << "XOR Properties:" << endl;
	cout << " -        Latency = " << xorLatency*1e9 << "ns" << endl;
	cout << " - Dynamic Energy = " << xorDynamicEnergy*1e12 << "pJ" << endl;
	cout << " -  Leakage Power = " << xorLeakage*1e3 << "mW" << endl;
}

Buffer & Buffer::operator=(const Buffer &rhs) {
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
	numRow = rhs.numRow;
	numColumn = rhs.numColumn;
	xorLatency = rhs.xorLatency;
	xorDynamicEnergy = rhs.xorDynamicEnergy;
	xorLeakage = rhs.xorLeakage;

	return *this;
}
