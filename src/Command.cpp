//$Header:$

/*
COMMAND.CPP - COMMAND Module for Geodesic Class

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
// Command.cpp: implementation of the CCommand class.
//
//////////////////////////////////////////////////////////////////////

#include "Command.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCommand::CCommand()
{
		freq = 3;				//Frequency of subdivison
		classt = 1;				//Class Type
		polyt = 1;				//Polyhedron Type
		filet = 0;				//File Output Type
		verbose_flag = 0;		//Suppress Topology display
		sphere_flag = 0;		//Generate sphere
		buckyball = 0;			//Generate buckyball
		faceflag = 1;			//Generate data files based on lines or polyfaces
		suppress_status = 0;	//Suppress display of calculation status
		E = 1.0;				//Elliptical eccentricity (0 - 1)
		ParabolicFlag = 0;		//Paraboloid flag
		ParabolicFocus = 0.5;		//Parabolic Focus
		ParabolicRadius = 1.0;	//Parabolic diameter defaults to diameter at focus
}

CCommand::~CCommand()
{

}
