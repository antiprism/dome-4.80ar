// Cartesian.h: interface for the CCartesian class.
//
///////////////////////////////// /////////////////////////////////////
//$Header:$

/*
CARTESIAN.H - CARTESIAN Module for Geodesic Class

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
#if !defined(AFX_CARTESIAN_H__FD415501_C3B1_11D1_BEB2_006067085732__INCLUDED_)
#define AFX_CARTESIAN_H__FD415501_C3B1_11D1_BEB2_006067085732__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

class CCartesian  
{
public:
	double chord(CCartesian &A, CCartesian &B);
	double clean_float(double fp_number);
	void SphericalToCartesian();

	double X;
	double Y;
	double Z;

	double theta;
	double phi;
	double radius;

	long a;
	long b;

	CCartesian();
	virtual ~CCartesian();

};

#endif // !defined(AFX_CARTESIAN_H__FD415501_C3B1_11D1_BEB2_006067085732__INCLUDED_)
