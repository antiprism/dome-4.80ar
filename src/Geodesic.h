//$Header:$

/*
GEODESIC.H - C++ Header file for geodesic.cpp geodesic dome class

	Copyright (C) 1995 - 2002 Richard J. Bono

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

	 Please direct inquiries, comments and modifications to:
	 Richard J. Bono
	 44 Augusta Rd.
	 Brownsville, TX 78521

	 email: rjbono@applied-synergetics.com

Revision history:

$Log:$
*/

#include<iostream>
#include<iomanip>
#include<fstream>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>

#include "Command.h"
#include "Labels.h"
#include "Points.h"
#include "Chords.h"
#include "Faces.h"
#include "ThetaTruncation.h"
#include "Cartesian.h"

const double RAD_TO_DEG = 57.2957795130824;
const double DEG_TO_RAD = .0174532925199433;

#define BYTE unsigned char
#define WORD unsigned int

//--------------geodesic points class
class CGeodesic	{
private:
	//Various indices, and temporary variables
	long i,	j, k;
	long index;
	double sphi, stheta, ephi, etheta, X, Y, Z, E;
	//actual values used in calculations set dependent on class type
	int status_count;

public:
	long PeakSym;
	long peakrow;
	int ParabolaFlag;
	double ParabolaRadius;
	double ParabolaFocus;
	long frequency, vertex, faces, edges, show_status;
	long sphere_flg, face_flag, ellipse_flag;
	long vertexII, facesII, edgesII;		//class II topology
	long bucky_face, bucky_vertex, bucky_edges, bucky_ball;
	long classtype, polytype;
	long freq_calc, vertex_calc, face_calc, edges_calc;

	//structure containing intermediate coordinates
	struct sphere{
		//vertex coordinates: 0,0,0 equals zenith
		long xprime;
		long yprime;
		long zprime;
	};

	//portability note: huge keyword was needed in order to allocate memory for
	//arrays from the far heap in a MSDOS far memory model.
	//No longer needed in 32-bit flat address space
	sphere  *pnt_calc;			//pointer to spherical point variable array
	CLabels *pnt_label;			//pointer to label array
	CPoints *pntcrd;  			//pointer to coordinate array
	CChords *edgepts;      		//pointer to chord array
	CFaces	*polyface;			//pointer to face array
	CPoints *sphere_pnt;		//pointer to sphere array
	CThetaTruncation *Truncation;	//Pointer to Peak Truncation info
	
	CGeodesic(CCommand &command);	//constructor
	~CGeodesic();		   						//destructor
	void topology(void);	   					//topological abundance function
	void spherical(void);	   				//calculate spherical coordinates
	double chord_length(double, double, double, double); //chord length function
	void chord_factor(void);   				//calculates chord factors
	void face_factor(void);						//calculates face factors
	void icosa_sphere(long);					//Make an icosa face
	void octa_sphere(long);						//Make an octa face
	void tetra_sphere(long);			  		//Make a Tetra face
	long tetra_angle(void);						//Get "A" coordinate to begin correction of bottom tetra face
	void save_dxf(char *);     				//save face data in DXF format
	void save_dxf_wire(char *);				//save DXF wireframe
	void save_buckydxf(char *);				//save bucky chords in DXF format
	void save_ascii(char *);   				//save all coordinate data in ASCII
	void save_POV(char *);						//save data in POV format
	void save_buckypov(char *);				//save buckyball in DXF format
	void save_PRN(char *);						//Save raw data in ASCII
	void save_WRL(char *);						//Save data in VRML 1.0 WRL format - Indexed-face sets
	void save_WRL_2(char *);						//Save data in VRML 2.0 WRL format - Indexed-face sets
	void save_WRL_wire(char *);				//Save data as indexed-line sets
	void save_WRL_wire_2(char *);				//Save data as indexed-line sets
	void save_buckywrl(char *);				//Save buckyball as wire-frame data
	void save_buckywrl_2(char *);				//Save buckyball as wire-frame data
	void save_OFF(char *);						//Save data in OFF format
	void save_buckyoff(char *);				//Save buckyball in OFF format
	void display_data(void);					//Display data during program execution
	double clean_float(double);				//function cleans up triag zeros
	double rotate_phi(double, double, double, double);	//phi rotation function
	double rotate_theta(double, double, double, double); //theta rotation function
	void time_passage(int);						//Show a time passage signal
	void bucky_factor(void);					//Generate Buckyball chord factors
	double ellipse_radius(double, double);	//calculate elliptical radius
	double root_E(double, double);			//root E theta correction for ellipses
	long face_quantity(long, long);			//determine the number of faces required for sphere
	double axial_angle_A(long);				//Determine Axial Angle A given chord index
	double axial_angle_B(long, double);		//Determine Axial Angle B given chord index & Angle A
	double face_angle_A(long);					//Determine Face Angle A given face index
	double face_angle_B(long, double);		//Determine Face Angle B Given face index & Face angle A
	double face_angle_C(double, double);	//Determine Face Angle C given Angle A & B
	double zenith(double, double);			//Calculate zenith height
	double floor(double, double);			//Calculate floor radius 
	double ParabolaTheta(double radius, double focus);
	double parabolic_radius(double theta, double focus);
	int TestChordCrossing(CCartesian &A, CCartesian &B, CCartesian &T);
	int TestParabolicIcosa(CCartesian &A, CCartesian &B, double &TruncTheta, double &ThetaLimit);
	double FindTruncLimit(double Trunc);
};
