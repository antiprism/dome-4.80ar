//$Header:$

/*
DXFSAVE.CPP - DXF Module for Geodesic Class

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

#include "Geodesic.h"

using namespace std;

//----------------------------------------DXF File memeber Functions
//-------------------save DXF file
void CGeodesic::save_dxf(char *filename)
{
	long end_face;
   // AR: set to 0.0 to remove warning that may be used unset
	double TruncTheta = 0.0;
   // AR: set to 0.0 to remove warning that may be used unset
	double ThetaLimit = 0.0;

	CCartesian A;
	CCartesian B;
	CCartesian C;

	//save Face data for symmetry triangle to DXF file
	//filename should include .DXF extension on MS-DOS systems
	//For full spheres each face is saved on a different level...

   ofstream DXF(filename);
   ofstream DUMP("parabola.txt");

	//output DXF data
	DXF << "0" << '\n';
	DXF << "SECTION" << '\n';
	DXF << "2" << '\n';
	DXF << "ENTITIES" << '\n';
	//Set field widths
	DXF << setiosflags(ios::fixed) << setw(8) << setprecision(6);
	DUMP << setiosflags(ios::fixed) << setw(8) << setprecision(6);

	//determine the number of faces required (first face = 0)
	if(sphere_flg)
		end_face = face_quantity(classtype, polytype);
	else
		end_face = 0;

	if(ParabolaFlag){
		if(ParabolaRadius == 2 * ParabolaFocus){
			TruncTheta = 90.0;
			ThetaLimit = 90.0;
		}
		else{
			TruncTheta = ParabolaTheta(ParabolaRadius, ParabolaFocus);
			ThetaLimit = FindTruncLimit(TruncTheta);
		}
	}

	if(!show_status){
		cout << "Saving Data to File... ";
		status_count = 0;
	}

	for(j=0; j<=end_face; j++){
		if(ParabolaFlag){
			DUMP << "FACE Type " << j << '\n';
			DUMP << "A" << '\t' << "B" << '\t' << "C" << '\t';;
			DUMP << "Chord a" << '\t' << "Chord b" << '\t' << "Chord c" << '\t' ;
			DUMP << "Face Angle A" << '\t' << "Face Angle B" << '\t' << "Face Angle C" << '\n';
		}

		//Get the vertex data for the face in question
		if(polytype == 1)
			icosa_sphere(j);
		else if(polytype == 2)
			octa_sphere(j);
		else if(polytype == 3)
			tetra_sphere(j);

		for(i=1; i<=face_calc; i++){
			if(!show_status){
				time_passage(status_count);
				status_count++;
				if(status_count > 3)
					status_count = 0;
			}
			//convert spherical to cartesian.
			//Start of chord
			if(ParabolaFlag){
				A.radius = sphere_pnt[polyface[i].cornerA].radius;
				A.phi = sphere_pnt[polyface[i].cornerA].phi;
				A.theta = sphere_pnt[polyface[i].cornerA].theta;
				B.radius = sphere_pnt[polyface[i].cornerB].radius;
				B.phi = sphere_pnt[polyface[i].cornerB].phi;
				B.theta = sphere_pnt[polyface[i].cornerB].theta;
				C.radius = sphere_pnt[polyface[i].cornerC].radius;
				C.phi = sphere_pnt[polyface[i].cornerC].phi;
				C.theta = sphere_pnt[polyface[i].cornerC].theta;

				if(A.theta > ThetaLimit && B.theta > ThetaLimit && C.theta > ThetaLimit)
					continue;

				else if(A.theta > TruncTheta && B.theta <= TruncTheta && C.theta <= TruncTheta){
					A.theta = TruncTheta;
					A.radius = parabolic_radius(A.theta, ParabolaFocus);
				}
				else if(A.theta <= TruncTheta && B.theta > TruncTheta && C.theta <= TruncTheta){
					B.theta = TruncTheta;
					B.radius = parabolic_radius(B.theta, ParabolaFocus);
				}
				else if(A.theta <= TruncTheta && B.theta <= TruncTheta && C.theta > TruncTheta){
					C.theta = TruncTheta;
					C.radius = parabolic_radius(C.theta, ParabolaFocus);
				}
				else if(A.theta > TruncTheta && B.theta > TruncTheta && C.theta <= TruncTheta){
					A.theta = TruncTheta;
					A.radius = parabolic_radius(A.theta, ParabolaFocus);
					B.theta = TruncTheta;
					B.radius = parabolic_radius(B.theta, ParabolaFocus);
				}
				else if(A.theta > TruncTheta && B.theta <= TruncTheta && C.theta > TruncTheta){
					A.theta = TruncTheta;
					A.radius = parabolic_radius(A.theta, ParabolaFocus);
					C.theta = TruncTheta;
					C.radius = parabolic_radius(C.theta, ParabolaFocus);

					A.SphericalToCartesian();
					B.SphericalToCartesian();
					C.SphericalToCartesian();

					if(TestChordCrossing(A, C, B) == 1)
						continue;

				}
				else if(A.theta <= TruncTheta && B.theta > TruncTheta && C.theta > TruncTheta){
					B.theta = TruncTheta;
					B.radius = parabolic_radius(B.theta, ParabolaFocus);
					C.theta = TruncTheta;
					C.radius = parabolic_radius(C.theta, ParabolaFocus);			
				}
				else if(A.theta > TruncTheta && B.theta > TruncTheta && C.theta > TruncTheta)
					continue;



				A.SphericalToCartesian();
				B.SphericalToCartesian();
				C.SphericalToCartesian();


				//Output to Dump File
				double clength = C.chord(A, B);
				double alength = A.chord(B, C);
				double blength = B.chord(C, A);

				double FaceA = acos((blength*blength + clength*clength - alength*alength)/(2*blength*clength)) * RAD_TO_DEG;
 				double FaceB = asin((blength * sin(FaceA * DEG_TO_RAD))/alength) * RAD_TO_DEG;
				double FaceC = 180.0 - FaceA - FaceB;
					
				DUMP << '\"' << pnt_label[polyface[i].cornerA].A << "-" << pnt_label[polyface[i].cornerA].B << '\"' << '\t';
				DUMP << '\"' << pnt_label[polyface[i].cornerB].A << "-" << pnt_label[polyface[i].cornerB].B << '\"' << '\t';
				DUMP << '\"' << pnt_label[polyface[i].cornerC].A << "-" << pnt_label[polyface[i].cornerC].B << '\"' << '\t';
				DUMP << alength << '\t' << blength << '\t' << clength << '\t' << FaceA << '\t' << FaceB << '\t' << FaceC << '\n';

			}
			else{

				//convert spherical to cartesian.
				//Face point A
				A.radius = sphere_pnt[polyface[i].cornerA].radius;
				A.phi = sphere_pnt[polyface[i].cornerA].phi;
				A.theta = sphere_pnt[polyface[i].cornerA].theta;
				A.SphericalToCartesian();

				//Face point B
				B.radius = sphere_pnt[polyface[i].cornerB].radius;
				B.phi = sphere_pnt[polyface[i].cornerB].phi;
				B.theta = sphere_pnt[polyface[i].cornerB].theta;
				B.SphericalToCartesian();

				//Face point C
				C.radius = sphere_pnt[polyface[i].cornerC].radius;
				C.phi = sphere_pnt[polyface[i].cornerC].phi;
				C.theta = sphere_pnt[polyface[i].cornerC].theta;
				C.SphericalToCartesian();
			}
			//save data
			DXF << "0" << '\n';
			DXF << "3DFACE" << '\n';	//Save as a 3D Polyface entity
			DXF << "8" << '\n';
			DXF << (j + 1) << '\n';		//Each face is saved on a different level.
			DXF << "62" << '\n';
			DXF << "3" << '\n';
			DXF << "10" << '\n';
			DXF << A.X << '\n';
			DXF << "20" << '\n';
			DXF << A.Y << '\n';
			DXF << "30" << '\n';
			DXF << A.Z << '\n';
			DXF << "11" << '\n';
			DXF << B.X << '\n';
			DXF << "21" << '\n';
			DXF << B.Y << '\n';
			DXF << "31" << '\n';
			DXF << B.Z << '\n';
			DXF << "12" << '\n';
			DXF << C.X << '\n';
			DXF << "22" << '\n';
			DXF << C.Y << '\n';
			DXF << "32" << '\n';
			DXF << C.Z << '\n';
			DXF << "13" << '\n';
			DXF << C.X << '\n';
			DXF << "23" << '\n';
			DXF << C.Y << '\n';
			DXF << "33" << '\n';
			DXF << C.Z << '\n';
	}

}

	cout << '\r' << "                     " << '\r'; //Clear status signal

	DXF << "0" << '\n';
	DXF << "ENDSEC" << '\n';
	DXF << "0" << '\n';
	DXF << "EOF" << '\n';

	DXF.close();
	DUMP.close();
}


//-------------------save DXF file in Wire-frame mode
void CGeodesic::save_dxf_wire(char *filename)
{
	long end_face;

	CCartesian A;
	CCartesian B;

	//save chord data for symmetry triangle to DXF file
	//filename should include .DXF extension on MS-DOS systems
	//For full spheres each face is saved on a different level...

	ofstream DXF(filename);
	//output DXF data
  	DXF << "0" << '\n';
	DXF << "SECTION" << '\n';
	DXF << "2" << '\n';
	DXF << "ENTITIES" << '\n';
	//Set field widths
	DXF << setiosflags(ios::fixed) << setw(8) << setprecision(6);
	//determine the number of faces required (first face = 0)
	if(sphere_flg)
		end_face = face_quantity(classtype, polytype);
	else
		end_face = 0;


	if(!show_status){
		cout << "Saving Data to File... ";
		status_count = 0;
	}

	for(j=0; j<=end_face; j++){

		//Get the vertex data for the face in question
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
			//convert spherical to cartesian.
			//Start of chord

			A.phi = sphere_pnt[edgepts[i].start].phi;
			B.phi = sphere_pnt[edgepts[i].end].phi;
			A.theta = sphere_pnt[edgepts[i].start].theta;
			B.theta = sphere_pnt[edgepts[i].end].theta;
			A.radius = sphere_pnt[edgepts[i].start].radius;
			B.radius = sphere_pnt[edgepts[i].end].radius;

			A.SphericalToCartesian();
			B.SphericalToCartesian();
			
			//save data
			DXF << "0" << '\n';
			DXF << "LINE" << '\n';	//Save as a 3D Polyface entity
			DXF << "8" << '\n';
			DXF << (j + 1) << '\n';	//Each face is saved on a different level.
			DXF << "62" << '\n';
			DXF << "3" << '\n';		//color number
			DXF << "10" << '\n';
			DXF << A.X << '\n';
			DXF << "20" << '\n';
			DXF << A.Y << '\n';
			DXF << "30" << '\n';
			DXF << A.Z << '\n';
			DXF << "11" << '\n';
			DXF << B.X << '\n';
			DXF << "21" << '\n';
			DXF << B.Y << '\n';
			DXF << "31" << '\n';
			DXF << B.Z << '\n';
		}		

	}
	cout << '\r' << "                     " << '\r'; //Clear status signal

	DXF << "0" << '\n';
	DXF << "ENDSEC" << '\n';
	DXF << "0" << '\n';
	DXF << "EOF" << '\n';

	DXF.close();
}

//-------------------save Buckyball data in DXF format
void CGeodesic::save_buckydxf(char *filename)
{
//	double Ax, Ay, Az, Bx, By, Bz;
	long end_face;

	CCartesian A;
	CCartesian B;

	//save chords data for buckyball chord data to DXF file
	//filename should include .DXF extension on MS-DOS systems

	ofstream DXF(filename);
	//output DXF data
	DXF << "0" << '\n';
	DXF << "SECTION" << '\n';
	DXF << "2" << '\n';
	DXF << "ENTITIES" << '\n';
	//Set field widths
	DXF << setiosflags(ios::fixed) << setw(8) << setprecision(6);



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

	for(j=0; j<=end_face; j++){
		//Get the vertex data for the face in question
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
				//start point of chord
				A.radius = sphere_pnt[edgepts[i].start].radius;
				A.phi = sphere_pnt[edgepts[i].start].phi;
				A.theta = sphere_pnt[edgepts[i].start].theta;
				A.SphericalToCartesian();

				//end point of chord
				B.radius = sphere_pnt[edgepts[i].start].radius;
				B.phi = sphere_pnt[edgepts[i].start].phi;
				B.theta = sphere_pnt[edgepts[i].start].theta;
				B.SphericalToCartesian();

				//save data
				DXF << "0" << '\n';
				DXF << "LINE" << '\n';
				DXF << "8" << '\n';
				DXF << (j + 1) << '\n';
				DXF << "62" << '\n';
				DXF << "3" << '\n';
				DXF << "10" << '\n';
				DXF << A.X << '\n';
				DXF << "20" << '\n';
				DXF << A.Y << '\n';
				DXF << "30" << '\n';
				DXF << A.Z << '\n';
				DXF << "11" << '\n';
				DXF << B.X << '\n';
				DXF << "21" << '\n';
				DXF << B.Y << '\n';
				DXF << "31" << '\n';
				DXF << B.Z << '\n';
		}
	}
	cout << '\r' << "                     " << '\r'; //Clear status signal

	DXF << "0" << '\n';
	DXF << "ENDSEC" << '\n';
	DXF << "0" << '\n';
	DXF << "EOF" << '\n';

	DXF.close();
}

//End of DXFSAVE
