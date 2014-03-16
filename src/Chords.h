//$Header:$

/*
CHORDS.H - chords Module for Geodesic Class

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
// Chords.h: interface for the CChords class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CHORDS_H__12384A43_B11A_11D1_ACEF_004095137C07__INCLUDED_)
#define AFX_CHORDS_H__12384A43_B11A_11D1_ACEF_004095137C07__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

class CChords  
{
public:
	long	start;		//chord start point
	long 	end;			//chord end point

	CChords();
	virtual ~CChords();

};

#endif // !defined(AFX_CHORDS_H__12384A43_B11A_11D1_ACEF_004095137C07__INCLUDED_)
