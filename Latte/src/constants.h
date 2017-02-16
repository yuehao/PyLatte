/*
 *  constants.h
 *  ERLtrac
 *
 *  Created by Yue Hao on 4/28/11.
 *  Copyright 2011 Brookhaven National Laboratory. All rights reserved.
 *
 */
 
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <cmath>

const double Mass_Electron =510999.0;
const double Mass_Proton =9.38272310e8;
const double Light_Speed=299792458.0;
const double Permittivity=8.854187817e-12;
const double Permeability=1.25663706143592e-6;
const double Unit_Charge=1.602176487e-19;
const double INFZERO=1e-15;
const double Two_Pi=2.0*M_PI;
const int Phase_Dimension=6;

enum PhaseSpaceCoordinate {
	x_, px_, y_, py_, ct_, dE_
};
enum SpaceCoordinate {
	X_, Y_, Z_
};
enum CellType {NOTYPE, DRIFT, DIPOLE ,QUADRUPOLE , SEXTRUPOLE, OCTUPOLE, MULTIPOLE, LINEARMAP, NONLINEARMAP, BEAMBEAM,
	RFCAVITY, MARKER, BPM, SPREADER, MERGER, SOLENOID, CORRECTOR, KICKER};


enum SimplecticInteger {SI_SECOND=2, SI_FOURTH=4};
enum MachineType {
	MT_LINEAR, MT_RING, MT_ERL
};
enum RANDOMTYPE {
	NORMAL, UNIFORM
};



#endif