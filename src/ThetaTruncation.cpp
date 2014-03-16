//$Header:$

/*
THETATRUNCATION.CPP - THETATRUNCATION Module for Geodesic Class

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
// ThetaTruncation.cpp: implementation of the CThetaTruncation class.
//
//////////////////////////////////////////////////////////////////////

#include "ThetaTruncation.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CThetaTruncation::CThetaTruncation()
{
	InitialPeakTheta = 0;
	InitialPeakA = -1;
	InitialPeakB = -1;

	PeakTheta = 0;
	PeakA = -1;
	PeakB = -1;

}

CThetaTruncation::~CThetaTruncation()
{

}
