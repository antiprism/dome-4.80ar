//$Header:$

/*
OFFSAVE.CPP - OFF (OFF) Module for Geodesic Class

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
   25/01/2011: Adrian Rossiter <adrian@antiprism.com>
               New file. Modified from Wrlfile.cpp. 

$Log:$

*/

#include "Geodesic.h"

using namespace std;

//----------------------------------------OFF File memeber Functions
//-------------------save OFF file as indexed face set
void CGeodesic::save_OFF(char *filename)
{
	long sX, sY, sZ;
	double dX, dY, dZ;
	long end_face;

	//save chord data for symmetry triangle to OFF file
	//filename should include .OFF extension on MS-DOS systems

	//determine the number of faces required (first face = 0)
	if(sphere_flg)
		end_face = face_quantity(classtype, polytype);
	else
		end_face = 0;

	if(!show_status){
		cout << "Saving Data to File... ";
		status_count = 0;
	}
	
   ofstream OFF(filename);

	//set up OFF file header
   OFF << "OFF\n";
   OFF << vertex_calc*(end_face+1) << " " << face_calc*(end_face+1) << " 0\n" ;

	//Set field widths
	OFF << setiosflags(ios::fixed) << setprecision(14);

	//output OFF data...start with vertex coordinates...
	for(j=0; j<=end_face; j++){
		if(polytype == 1)
			icosa_sphere(j);
		else if(polytype == 2)
			octa_sphere(j);
		else if(polytype == 3)
			tetra_sphere(j);

		for(i=1; i<=vertex_calc; i++){

			if(!show_status){
				time_passage(status_count);
				status_count++;
				if(status_count > 3)
					status_count = 0;
			}
			if(ParabolaFlag){
				double TruncTheta;
				TruncTheta = ParabolaTheta(ParabolaRadius, ParabolaFocus);
				if(sphere_pnt[i].theta > TruncTheta)
					continue;			
			}
			//convert spherical to cartesian
			dX = sphere_pnt[i].radius *
				  clean_float(cos(sphere_pnt[i].phi * DEG_TO_RAD) *
								  sin(sphere_pnt[i].theta * DEG_TO_RAD));
			dY = sphere_pnt[i].radius *
				  clean_float(sin(sphere_pnt[i].phi * DEG_TO_RAD) *
								  sin(sphere_pnt[i].theta * DEG_TO_RAD));
			dZ = sphere_pnt[i].radius *
				  clean_float(cos(sphere_pnt[i].theta * DEG_TO_RAD));

			//Save data
			OFF << dX << " " << dY << " " << dZ << '\n';
		}
	}
	for(j=0; j<=end_face; j++){
		//output OFF face data...
		for(i=1; i<=face_calc; i++){
			if(!show_status){
				time_passage(status_count);
				status_count++;
				if(status_count > 3)
					status_count = 0;
			}
				if(ParabolaFlag){
					double TruncTheta;
					TruncTheta = ParabolaTheta(ParabolaRadius, ParabolaFocus);
					if(sphere_pnt[edgepts[i].start].theta > TruncTheta)
						continue;			
					else if(sphere_pnt[edgepts[i].end].theta > TruncTheta)
						continue;					
				}
			//Data saved is the point coordinate index for the given face
			sX = (polyface[i].cornerA + (vertex_calc * j)) - 1;
			sY = (polyface[i].cornerB + (vertex_calc * j)) - 1;
			sZ = (polyface[i].cornerC + (vertex_calc * j)) - 1;
         OFF << "3 " << sX << " " << sY << " " << sZ << '\n';
		}
   }
}


void CGeodesic::save_buckyoff(char *filename)
{
	double sX, sY, sZ, eX, eY, eZ;
	long end_face;

	//save Buckyball chord data to POV file
	//filename should include .POV extension on MS-DOS systems

	//determine the number of faces required (first face = 0)
	if(sphere_flg)
		end_face = face_quantity(classtype, polytype);
	else
		end_face = 0;

	edges_calc = bucky_edges;

	if(!show_status){
		cout << "Saving Data to File... ";
		status_count = 0;
	}

   ofstream OFF(filename);

	//set up OFF file header
   OFF << "OFF\n";
   int v_cnt = edges_calc*(end_face+1)*2;
   OFF << v_cnt << " " << v_cnt/2 << " 0\n" ;
   
	//Set field widths
	OFF << setiosflags(ios::fixed) << setprecision(14);

	//output OFF data...start with vertex coordinates...
	for(j=0; j<=end_face; j++){

		if(polytype == 1)
			icosa_sphere(j);
		else if(polytype == 2)
			octa_sphere(j);
		else if(polytype == 3)
			tetra_sphere(j);

		for(i=1; i<=edges_calc; i++){
			if(!show_status){
				time_passage(status_count);
				status_count++;
				if(status_count > 3)
					status_count = 0;
			}
			//convert spherical to cartesian
			sX = sphere_pnt[edgepts[i].start].radius *
					clean_float(cos(sphere_pnt[edgepts[i].start].phi*DEG_TO_RAD) *
								   sin(sphere_pnt[edgepts[i].start].theta*DEG_TO_RAD));
			sY = sphere_pnt[edgepts[i].start].radius *
					clean_float(sin(sphere_pnt[edgepts[i].start].phi*DEG_TO_RAD) *
									sin(sphere_pnt[edgepts[i].start].theta*DEG_TO_RAD));
			sZ = sphere_pnt[edgepts[i].start].radius *
					clean_float(cos(sphere_pnt[edgepts[i].start].theta*DEG_TO_RAD));

			eX = sphere_pnt[edgepts[i].end].radius *
					clean_float(cos(sphere_pnt[edgepts[i].end].phi*DEG_TO_RAD) *
									sin(sphere_pnt[edgepts[i].end].theta*DEG_TO_RAD));
			eY = sphere_pnt[edgepts[i].end].radius *
					clean_float(sin(sphere_pnt[edgepts[i].end].phi*DEG_TO_RAD) *
									sin(sphere_pnt[edgepts[i].end].theta*DEG_TO_RAD));
			eZ = sphere_pnt[edgepts[i].end].radius *
					clean_float(cos(sphere_pnt[edgepts[i].end].theta*DEG_TO_RAD));

			//Save data
			OFF << sX << " " << sY << " " << sZ << '\n';
			OFF << eX << " " << eY << " " << eZ << '\n';
		}
	}

   //Save edge data
   for(int i=0; i<v_cnt/2; i++)
      OFF << "2 " << i*2 << " " << i*2+1 << '\n';
}

//End of OFFSAVE
