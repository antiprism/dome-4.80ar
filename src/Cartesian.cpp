//$Header:$

/*
CARTESIAN.CPP - CARTESIAN Module for Geodesic Class

	Copyright (C) 1998 - 2002 Richard J. Bono

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

// Cartesian.cpp: implementation of the CCartesian class.
//
//////////////////////////////////////////////////////////////////////


#include "Cartesian.h"
#include "Geodesic.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCartesian::CCartesian()
{
	X = 0;
	Y = 0;
	Z = 0;

	phi = 0;
	theta = 0;
	radius = 0;
}

CCartesian::~CCartesian()
{

}

void CCartesian::SphericalToCartesian()
{
	//Convert Spherical Coordinates to Cartesian
	X = radius * clean_float(cos(phi * DEG_TO_RAD) * sin(theta * DEG_TO_RAD));
	Y = radius * clean_float(sin(phi * DEG_TO_RAD) * sin(theta * DEG_TO_RAD));
	Z = radius * clean_float(cos(theta * DEG_TO_RAD));

	return;

}

double CCartesian::clean_float(double fp_number)
{
		//function cleans floating point results close to but not quite zero
		if((fp_number < 0.0) && (fp_number > -0.0000000000001))
			fp_number = fabs(fp_number);

		return fp_number;

}

double CCartesian::chord(CCartesian & A, CCartesian & B)
{
	
	
	double Asqr = pow(B.X-A.X, 2.0);
	double Bsqr = pow(B.Y-A.Y, 2.0);
	double Csqr = pow(B.Z-A.Z, 2.0);

	double result = pow(Asqr+Bsqr+Csqr, 0.5);


	return result;
}
