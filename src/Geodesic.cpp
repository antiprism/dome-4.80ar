//$Header:$

/*
GEODESIC.CPP - C++ class for the calculation of geodesic dome properties

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

Acknowledgements & References:
The main reference used in the creation of this code was "Geodesic Math & How
to Use It" by Hugh Kenner, 1976, University of California Press.
ISBN 0-520-02924-0; Library of Congress Catalog Card Number: 74-27292. Many
thanks to Hugh for putting this data in an accessible format.

Also, many thanks to:
	J. F. Nystrom
	My wife & Daughters
	Chris Fearnley
	Kirby Urner
		&
	R. Buckminster Fuller for changing the way I view Universe.

*/

#include "Geodesic.h"

using namespace std;

//--------------geodesic constructor
CGeodesic::CGeodesic(CCommand &command)
{

	//The constructor handles all calculation when an instance is called.
	//The input parameters are the frequency of subdivision, the class
	//type, the polygon type and a flag to determine whether to generate
	//spherical data. This is all that is needed to specify a geodesic dome!
	//(Note: I have not included radius as size is special case and scaling
	//factors can be applied. The program assumes unit radius domes).

	//initialize class variables
	frequency = command.freq;
	classtype = command.classt;
	polytype = command.polyt;
	sphere_flg = command.sphere_flag;
	show_status = command.suppress_status;
	bucky_ball = command.buckyball;
	face_flag = command.faceflag;
	ParabolaFlag = command.ParabolicFlag;
	ParabolaFocus = command.ParabolicFocus;
	ParabolaRadius = command.ParabolicRadius;
	vertex = 0;
	faces = 0;
	edges = 0;

	if(command.E != 1.0){
		//Elliptical Structure formation
		ellipse_flag = 1;
		E = command.E;
	}
	else
		ellipse_flag = 0;

	if(ParabolaFlag == 1){
		if(ParabolaFocus/(2*ParabolaRadius) < 0.25){
			cout << "f/d ratio must be greater than 0.25 --- Execution Terminating" << '\n';
			exit(1);
		}	

		ellipse_flag = 0;
		if(ParabolaRadius == 0.0)
			ParabolaRadius = 2 * ParabolaFocus;
	}

}
//--------------geodesic destructor
CGeodesic::~CGeodesic() {}
//--------------calculate topological abundance of symmetry triangle
void CGeodesic::topology(void)
{

	//This class function calculates the topological abundance
	//of the symmetry triangle.

	//calculate the number of vertexia & faces in symmetry triangle
	//This is based on the class type
	//First the polyhedron face...

	for(i=1; i<=frequency; i++){
		if(i == 1)
			vertex = 3;
		else
			vertex += (i + 1);
	}
	faces = frequency * frequency;
	edges = (vertex + faces) - 1;

	//if class II then calculate the class II symmetry triangle topology.
	//The class II triangle is NOT the polyhedron face!
	if(classtype == 2){
		for(i=1; i<=frequency/2; i++){
			if(i == 1)
				vertexII = 3;
			else
				vertexII += (i + 1);
		}
		facesII = (frequency/2) * (frequency/2);
		edgesII = (vertexII + facesII) - 1;
	}

	//Calculate the abundance of the Icosa face if Buckyball generation is enabled.
	if(bucky_ball){
		bucky_vertex = 0;
		for(i=0; i<=frequency; i+=3)
			bucky_vertex += (2 * i);

		bucky_face = 1;
		for(i=3; i<=frequency; i+=3)
			bucky_face += (frequency - i);

			bucky_edges = (bucky_vertex + bucky_face) - 1;
	}
	//Set the topology dependent on class type
	if(classtype == 1){
		vertex_calc = vertex;
		edges_calc = edges;
		freq_calc = frequency;
		face_calc = faces;
	}
	else if(classtype == 2){
		vertex_calc =vertexII;
		edges_calc = edgesII;
		freq_calc = frequency / 2;
		face_calc = facesII;
	}
}

//-------------calculate spherical coordinates for given polygon type & class
void CGeodesic::spherical(void)
{
	//calculate spherical coordinates of geodesic vertexia
	//given coordinates values
	//care must be taken to avoid generating a trig error

	//Memory is allocated for the full polyhedron face. Create array object for vertexia
	//First allocate memory for label array

	pnt_label = new CLabels[vertex_calc + 1];
	if (pnt_label == NULL){
		cout << "Insufficient memory for label array --- Execution Terminating" << '\n';
		exit(1);
	}

	//now allocate space for the actual vertex coordinate array
	pntcrd = new CPoints[vertex_calc + 1];
	if (pntcrd == NULL){
		cout << "Insufficient memory for coordinate array --- Execution Terminating" << '\n';
		exit(1);
	}

	//create spherical array for one polyhedron face
	sphere_pnt = new CPoints[vertex_calc + 1];
	if (sphere_pnt == NULL){
		cout << "Insufficient memory for face array --- Execution Terminating" << '\n';
		exit(1);
	}

	//now allocate space for temporary variables
	pnt_calc = new sphere[vertex_calc + 1];
	if (pnt_calc == NULL){
		cout << "Insufficient memory for Spherical Array --- Execution Terminating" << '\n';
		exit(1);
	}

	//Initialize vertex labels and coordinates
	index = 1;
	for(i=0; i<=freq_calc; i++){
		for(j=0; j<=i; j++){
			pnt_label[index].A = i;
			pnt_label[index].B = j;
			pnt_calc[index].xprime = j;
			pnt_calc[index].yprime = i - j;
			pnt_calc[index].zprime = frequency - i;
			index++;
		}
	}

	if(!show_status){
		cout << "Calculating Spherical Coordinates... ";
		status_count = 0;
	}

	for(i=1; i<=vertex_calc; i++){
		if(!show_status){
			time_passage(status_count);
			status_count++;
			if(status_count > 3)
				status_count = 0;
		}
		//set-up x,y,z coordinates for spherical calulations
		//See page 75 of Kenner's "Geodesic Math & How to use it"
		//for formulae used below
		if(polytype == 2){
				//Octahedron
				X = pnt_calc[i].xprime;
				Y = pnt_calc[i].yprime;
				Z = pnt_calc[i].zprime;
		}else if(polytype == 1){
				//Icosahedron
				X = pnt_calc[i].xprime * sin(72.0 * DEG_TO_RAD);
				Y = pnt_calc[i].yprime + (pnt_calc[i].xprime * cos(72.0 * DEG_TO_RAD));
				Z = frequency/2.0 + (pnt_calc[i].zprime / (2.0 * cos(36.0 * DEG_TO_RAD)));
		}else if(polytype == 3){
				//Tetrahedron
				X = pnt_calc[i].xprime * pow(3.0, 0.5);
				Y = 2 * pnt_calc[i].yprime - pnt_calc[i].xprime;
				Z = (3.0 * pnt_calc[i].zprime - pnt_calc[i].xprime - pnt_calc[i].yprime) / pow(2.0, 0.5);
		}


		//Calculate phi while avoiding trig errors
		if(Y == 0 && X == 0)
			pntcrd[i].phi = 0.0;
		else if(Y == 0)
			pntcrd[i].phi = 90.0;
		else
			pntcrd[i].phi = atan(X / Y) * RAD_TO_DEG;

		//adjust value to correct quadrant
		if(pntcrd[i].phi < 0.0)
			pntcrd[i].phi += 180.0;

		//now calculate theta...this is class dependent...
		if(Z == 0)
			pntcrd[i].theta = 90.0;
		else if(classtype == 1)
			//All class I types use this equation for theta
			pntcrd[i].theta = atan((pow(pow(X, 2.0) + pow(Y, 2.0), 0.5)) / Z) * RAD_TO_DEG;
		else if(polytype == 2 && classtype == 2)
			//theta for class II octahedra
			pntcrd[i].theta = atan((pow((2.0 * (pow(X, 2.0) + pow(Y, 2.0))), 0.5)) / Z) * RAD_TO_DEG;
		else if(polytype == 1 && classtype == 2)
			//theta for class II Icosahedron
			pntcrd[i].theta = atan((pow(pow(X, 2.0) + pow(Y, 2.0), 0.5)) / (cos(36.0 * DEG_TO_RAD) * Z)) * RAD_TO_DEG;
		else if(polytype == 3 && classtype ==2)
			//theta for class II tetrahedron
			//The formula given in Geodesic math (eq 12.15) appears to be incorrect. The formula
			//below provides results indicated on tables
			pntcrd[i].theta = atan((2 * pow((pow(X, 2.0) + pow(Y, 2.0)), 0.5) / Z)) * RAD_TO_DEG;

		//make sure the right quadrant is used
		if(pntcrd[i].theta < 0.0)
			pntcrd[i].theta += 180.0;

		//Set Radius = 1
		pntcrd[i].radius = 1.0;
	}
	cout << '\r' << "                                     " << '\r'; //Clear status signal
	delete pnt_calc;				//release some memory
}
//-------------------calculate chord factors
void CGeodesic::chord_factor(void)
{
	//Function cycles through each vertex and determines the
	//chord connections. The use of linked lists have simplified this code
	//considerably.

	index = 1;
	if(!show_status){
		cout << "Calculating Chord Factors... ";
		status_count = 0;
	}

	//create array object for edges
	edgepts = new CChords[(edges + 1)];
	if (edgepts == NULL){
		cout << "Insufficient memory for chord array --- Execution Terminating" << '\n';
		exit(1);
	}

	for(i=1; i<vertex_calc; i++){
		if(!show_status){
			time_passage(status_count);
			status_count++;
			if(status_count > 3)
				status_count = 0;
		}
		if(pnt_label[i].A == 0 && pnt_label[i].B == 0){
			//point is the zenith vertex; add 10 & 11
			//line #1 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;

			//line #2 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//end point
			edgepts[index].end =j;

			index++;

		}
		else if(pnt_label[i].A == freq_calc && pnt_label[i].B == 0){
			//point is a left vertex; add 01
			//line #1 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//end point
			edgepts[index].end = j;
			index++;
		}
		else if(pnt_label[i].A == pnt_label[i].B){
			//point is a right edge; add 10 & 11
			//line #1 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;

			//line #2 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;
		}
		else if(pnt_label[i].A == freq_calc){
			//point is a bottom vertex; add 01
			//line #1 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A  && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;
		}
		else if(pnt_label[i].B == 0){
			//point is a right edge; add 01, 10 & 11
			//line #1 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;

			//line #2 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;

			//line #3 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;
		}
		else {
			//point is an interior vertex; add 01, 10 & 11
			//line #1 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;

			//line #2 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;

			//line #3 start point definition
			edgepts[index].start = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//end point
			edgepts[index].end = j;

			index++;
		}
	}
	cout << '\r' << "                             " << '\r'; //Clear status signal
}
//-------------------calculate face factors
void CGeodesic::face_factor(void)
{
	//Function cycles through each vertex and determines the
	//face definitions essentially going CCW about the face

	index = 1;
	if(!show_status){
		cout << "Calculating Face Factors... ";
		status_count = 0;
	}

	//create array object for faces
	polyface = new CFaces[face_calc + 1];
	if(polyface == NULL){
		cout << "Insufficient memory for face array --- Execution Terminating" << '\n';
		exit(1);
	}


	for(i=1; i<=(vertex_calc-(freq_calc+1)); i++){
	//function only interates to the end of line frequency - 1. This is located
	//at vertex#:vertex_calc-(freq_calc+1)
		if(!show_status){
			time_passage(status_count);
			status_count++;
			if(status_count > 3)
				status_count = 0;
		}
		if(pnt_label[i].A == 0 && pnt_label[i].B == 0){
			//point is the zenith vertex; add 10 & 01;
			//Corner A
			polyface[index].cornerA = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
					break;
			}
			//Corner B
			polyface[index].cornerB = j;
			for(k=i; k<=vertex_calc; k++){
				if(pnt_label[k].A == pnt_label[j].A && pnt_label[k].B == pnt_label[j].B + 1)
					break;
			}
			//Corner C
			polyface[index].cornerC = k;

			index++;

		}
		else if(pnt_label[i].B == 0){
			//point is a left edge; add 10 & 01;
			//                      add 11; subtract 10
			//face #1
			//Corner A
			polyface[index].cornerA = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
					break;
			}
			//corner B
			polyface[index].cornerB = j;
			for(k=i; k<=vertex_calc; k++){
				if(pnt_label[k].A == pnt_label[j].A && pnt_label[k].B == pnt_label[j].B + 1)
					break;
			}
			//corner C
			polyface[index].cornerC = k;

			index++;

			//Face #2
			//Corner A
			polyface[index].cornerA = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//Corner B
			polyface[index].cornerB = j;
			for(k=i; k<=vertex_calc; k++){
				if(pnt_label[k].A == pnt_label[j].A - 1 && pnt_label[k].B == pnt_label[j].B)
					break;
			}
			//Corner C
			polyface[index].cornerC = k;

			index++;

		}
		else if(pnt_label[i].A == pnt_label[i].B){
			//point is a right edge; add 10 & 01; subtract 11
			//line #1 start point definition
			//Corner A
			polyface[index].cornerA = i;

			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
					break;
			}
			//Corner B
			polyface[index].cornerB = j;
			for(k=i; k<=vertex_calc; k++){
				if(pnt_label[k].A == pnt_label[j].A && pnt_label[k].B == pnt_label[j].B + 1)
					break;
			}
			//Corner C
			polyface[index].cornerC = k;

			index++;

		}
		else {
			//Interior point; add 10 & 01;
			//                add 11; subtract 10
			//face #1
			//Corner A
			polyface[index].cornerA = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
					break;
			}
			//Corner B
			polyface[index].cornerB = j;
			for(k=i; k<=vertex_calc; k++){
				if(pnt_label[k].A == pnt_label[j].A && pnt_label[k].B == pnt_label[j].B + 1)
					break;
			}
			//Corner C
			polyface[index].cornerC = k;

			index++;

			//Face #2
			//Corner A
			polyface[index].cornerA = i;
			for(j=i; j<=vertex_calc; j++){
				if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
					break;
			}
			//Corner B
			polyface[index].cornerB = j;
			for(k=i; k<=vertex_calc; k++){
				if(pnt_label[k].A == pnt_label[j].A - 1 && pnt_label[k].B == pnt_label[j].B)
					break;
			}
			//Corner C
			polyface[index].cornerC = k;

			index++;
		}

		//faces=index;
	}
	cout << '\r' << "                             " << '\r'; //Clear status signal
}
//--------------------clean floating point numbers
double CGeodesic::clean_float(double fp_number)
{
		//function cleans floating point results close to but not quite zero
		if((fp_number < 0.0) && (fp_number > -0.0000000000001))
			fp_number = fabs(fp_number);

		return fp_number;
}

//-------------------time passage signal
void CGeodesic::time_passage(int count)
{
	//This function displays the twirly to incidate the program is working
	//and that you machine is not hanging somewhere in limbo...
	if(count == 0)
		cout << '\b' << "-";
	else if(count == 1)
		cout << '\b' << "/";
	else if(count == 2)
		cout << '\b' << "|";
	else if(count == 3)
		cout << '\b' << '\\';
}

//-------------------determine number of faces per sphere
long CGeodesic::face_quantity(long class_value, long poly_value)
{
	long return_value = -1; // AR: set to dummy value

	//determine the number of faces required (first face = 0)
	if(class_value == 1){
		if(poly_value == 1)
			return_value = 19;
		else if(poly_value == 2)
			return_value =7;
		else if(poly_value == 3)
			return_value = 3;
	}
	else if(class_value == 2){
		if(poly_value == 1)
			return_value = 59;
		else if(poly_value == 2)
			return_value =39;
		else if(poly_value == 3)
			return_value = 11;
	}

	return return_value;
}

//-------------------save ASCII PRN format
void CGeodesic::save_PRN(char *filename)
{
	double sX, sY, sZ, eX, eY, eZ;

	//This function saves the spherical & cartesian data to an
	//comma-delimited ASCII PRN file. This file can be imported into
	//most spreadsheets and the like.

	ofstream PRN(filename);

	//Set field widths
	PRN << setiosflags(ios::fixed) << setw(8) << setprecision(6);

	PRN << "Spherical Vertexia Coordinates" << '\n';
	//output PRN data...start with vertexia in spherical coordinates...
	for(i=1; i<=vertex_calc; i++)
			PRN << pntcrd[i].phi << ", " << pntcrd[i].theta << '\n';

	PRN << '\n';

	PRN << "Cartesian Vertexia Coordinates" << '\n';
	//...then in XYZ coordinates
	for(i=1; i<=vertex_calc; i++){
			//convert spherical to cartesian
			sX = clean_float(cos(pntcrd[i].phi * DEG_TO_RAD) * sin(pntcrd[i].theta * DEG_TO_RAD));
			sY = clean_float(sin(pntcrd[i].phi * DEG_TO_RAD) * sin(pntcrd[i].theta * DEG_TO_RAD));
			sZ = clean_float(cos(pntcrd[i].theta * DEG_TO_RAD));
			//Save data
			PRN << sX << ", " << sY << ", " << sZ << '\n';
	}
	PRN << '\n';

	PRN << "Spherical Chord Coordinates" << '\n';
	//Now do the chords...first in spherical coordinate pairs...
	for(i=1; i<=edges_calc; i++)
		PRN << pntcrd[edgepts[i].start].phi << ", " << pntcrd[edgepts[i].start].theta << ", "
				<< pntcrd[edgepts[i].end].phi << ", " << pntcrd[edgepts[i].end].theta << '\n';

	PRN << '\n';

	PRN << "Cartesian Chord Coordinates" << '\n';
	//...then in XYZ...
	for(i=1; i<=edges_calc; i++){
		//convert spherical to cartesian
		sX = clean_float(cos(pntcrd[edgepts[i].start].phi * DEG_TO_RAD) *
							  sin(pntcrd[edgepts[i].start].theta * DEG_TO_RAD));
		sY = clean_float(sin(pntcrd[edgepts[i].start].phi * DEG_TO_RAD) *
							  sin(pntcrd[edgepts[i].start].theta * DEG_TO_RAD));
		sZ = clean_float(cos(pntcrd[edgepts[i].start].theta * DEG_TO_RAD));
		eX = clean_float(cos(pntcrd[edgepts[i].end].phi * DEG_TO_RAD) *
							  sin(pntcrd[edgepts[i].end].theta * DEG_TO_RAD));
		eY = clean_float(sin(pntcrd[edgepts[i].end].phi * DEG_TO_RAD) *
							  sin(pntcrd[edgepts[i].end].theta * DEG_TO_RAD));
		eZ = clean_float(cos(pntcrd[edgepts[i].end].theta * DEG_TO_RAD));

		//save data
		PRN << sX << ", " << sY << ", " << sZ << ", "
				 << eX << ", " << eY << ", " << eZ << '\n';
	}
	PRN.close();
}

//-------------------make an icosa face for sphere
void CGeodesic::icosa_sphere(long facenumber)
{
	double x_phi, x_theta, rot_factor;

	//Given the polyhedron face number calculate Icosa
	//face and store in array. Used in creating spheres
	if(classtype == 1){
		if(facenumber >=0 && facenumber <=4){
			//do a 72 degree "cap" rotation
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = pntcrd[i].phi + 72.0 * facenumber;
					sphere_pnt[i].theta = pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >= 5 && facenumber <=9){
			//do a downward  rotation then rotate appropriate angle
			k = facenumber - 5;		//phi rotation factor
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].theta = rotate_theta(36.0, pntcrd[i].phi, (180.0 - (atan(2.0) * RAD_TO_DEG)),
							pntcrd[i].theta);
					sphere_pnt[i].phi = (rotate_phi(36.0, pntcrd[i].phi, pntcrd[i].theta,
							sphere_pnt[i].theta) + 36.0) + 72.0 * k;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >= 10 && facenumber <=14){
			//do a downward rotation then rotate 72° x k phi and 180° theta
			k = facenumber - 10;		//phi rotation factor
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].theta = 180.0 - rotate_theta(36.0, pntcrd[i].phi, (180.0 - (atan(2.0) * RAD_TO_DEG)),
							pntcrd[i].theta);
					sphere_pnt[i].phi = ((rotate_phi(36.0, pntcrd[i].phi, pntcrd[i].theta,
							sphere_pnt[i].theta) + 36.0) + 72.0 * k) + 36.0;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);

				}
			}
		}
		else{
			//rotate bottom cap into place and the sphere is done!
			k = facenumber - 15;
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0) || !bucky_ball){
					sphere_pnt[i].phi = (pntcrd[i].phi + 72.0 * k) + 36.0;
					sphere_pnt[i].theta = 180.0 - pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
	}
	else if(classtype == 2){
		if(facenumber >=0 && facenumber <=4){
			//do a 72 degree "cap" rotation
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = pntcrd[i].phi + 72.0 * facenumber;
					sphere_pnt[i].theta = pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >= 5 && facenumber <= 29){
			if(facenumber >= 5 && facenumber <= 9){
				k = facenumber - 5;
				rot_factor = 0;
			}
			else if(facenumber >= 10 && facenumber <= 14){
				k = facenumber - 10;
				rot_factor = 1;
			}
			else if(facenumber >= 15 && facenumber <= 19){
				k = facenumber - 15;
				rot_factor = 2;
			}
			else if(facenumber >= 20 && facenumber <= 24){
				k = facenumber - 20;
				rot_factor = 3;
			}
			else if(facenumber >= 25 && facenumber <= 29){
				k = facenumber - 25;
				rot_factor = 4;
			}
			//do a 72 degree cap rotation then rotate appropriate angle
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					x_phi = pntcrd[i].phi + 72.0 * k;
					x_theta = pntcrd[i].theta;

					sphere_pnt[i].theta = rotate_theta(36.0 + (72.0 * rot_factor), x_phi, (atan(2.0) * RAD_TO_DEG),
							x_theta);
					sphere_pnt[i].phi = (rotate_phi(36.0 + (72.0 * rot_factor), x_phi, x_theta,
							sphere_pnt[i].theta) + 36.0) + 72.0 * k;

					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >=30 && facenumber <=34){
			//do a 72 degree "cap" rotation
			k = facenumber - 30;
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = 36.0 + pntcrd[i].phi + 72.0 * k;
					sphere_pnt[i].theta = 180.0 - pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >= 35 && facenumber <= 59){
			if(facenumber >= 35 && facenumber <= 39){
				k = facenumber - 35;
				rot_factor = 0;
			}
			else if(facenumber >= 40 && facenumber <= 44){
				k = facenumber - 40;
				rot_factor = 1;
			}
			else if(facenumber >= 45 && facenumber <= 49){
				k = facenumber - 45;
				rot_factor = 2;
			}
			else if(facenumber >= 50 && facenumber <= 54){
				k = facenumber - 50;
				rot_factor = 3;
			}
			else if(facenumber >= 55 && facenumber <= 59){
				k = facenumber - 55;
				rot_factor = 4;
			}
			//do a 72 degree cap rotation then rotate appropriate angle
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					x_phi = pntcrd[i].phi + 72.0 * k;
					x_theta = pntcrd[i].theta;

					sphere_pnt[i].theta = 180.0 - rotate_theta(36.0 + (72.0 * rot_factor), x_phi, (atan(2.0) * RAD_TO_DEG),
							x_theta);
					sphere_pnt[i].phi = 36.0 + (rotate_phi(36.0 + (72.0 * rot_factor), x_phi, x_theta,
							sphere_pnt[i].theta) + 36.0) + 72.0 * k;

					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
	}
}
//-------------------make an octa face for sphere
void CGeodesic::octa_sphere(long facenumber)
{

	double x_phi, x_theta, rot_factor;

	//Given the polyhedron face number calculate Icosa
	//face and store in array. Used in creating spheres.
	if(classtype == 1){
		if(facenumber >=0 && facenumber <=3){
			//Build a northern hemisphere face
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = pntcrd[i].phi + 90.0 * facenumber;
					sphere_pnt[i].theta = pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >= 4 && facenumber <=7){
			//Rotate and shift southern hemisphere face
			k = facenumber - 4;					//phi rotation factor
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = pntcrd[i].phi + 90.0 * k;
					sphere_pnt[i].theta = 180.0 - pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
					
				}
			}
		}
	}
	else if(classtype == 2){
		if(facenumber >=0 && facenumber <=3){
			//Build a northern hemisphere cap face
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = pntcrd[i].phi + 90.0 * facenumber;
					sphere_pnt[i].theta = pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
					
				}
			}
		}
		else if(facenumber >= 4 && facenumber <= 19){
			if(facenumber >= 4 && facenumber <= 7){
				k = facenumber - 4;
				rot_factor = 0;
			}
			else if(facenumber >= 8 && facenumber <= 11){
				k = facenumber - 8;
				rot_factor = 1;
			}
			else if(facenumber >= 12 && facenumber <= 15){
				k = facenumber - 12;
				rot_factor = 2;
			}
			else if(facenumber >= 16 && facenumber <= 19){
				k = facenumber - 16;
				rot_factor = 3;
			}
			//do a 90 degree cap rotation then rotate appropriate angle
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					x_phi = pntcrd[i].phi + 90.0 * k;
					x_theta = pntcrd[i].theta;

					sphere_pnt[i].theta = rotate_theta(45.0 + (90.0 * rot_factor), x_phi, 90.0,
							x_theta);
					sphere_pnt[i].phi = (rotate_phi(45.0 + (90.0 * rot_factor), x_phi, x_theta,
							sphere_pnt[i].theta) + 45.0) + 90.0 * k;

					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >=20 && facenumber <=23){
			//Build a northern hemisphere cap face
			k = facenumber - 20;
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = pntcrd[i].phi + 90.0 * k;
					sphere_pnt[i].theta = 180.0 - pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >= 24 && facenumber <= 39){
			if(facenumber >= 24 && facenumber <= 27){
				k = facenumber - 24;
				rot_factor = 0;
			}
			else if(facenumber >= 28 && facenumber <= 31){
				k = facenumber - 28;
				rot_factor = 1;
			}
			else if(facenumber >= 32 && facenumber <= 35){
				k = facenumber - 32;
				rot_factor = 2;
			}
			else if(facenumber >= 36 && facenumber <= 39){
				k = facenumber - 36;
				rot_factor = 3;
			}
			//do a 90 degree cap rotation then rotate appropriate angle
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					x_phi = pntcrd[i].phi + 90.0 * k;
					x_theta = pntcrd[i].theta;

					sphere_pnt[i].theta = 180.0 - rotate_theta(45.0 + (90.0 * rot_factor), x_phi, 90.0,
							x_theta);
					sphere_pnt[i].phi = (rotate_phi(45.0 + (90.0 * rot_factor), x_phi, x_theta,
							sphere_pnt[i].theta) + 45.0) + 90.0 * k;

					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}

	}
}
//-------------------make a tetra face for sphere
void CGeodesic::tetra_sphere(long facenumber)
{
	double rot_factor;

	//find A coordinate where angles turn
	long A_change = tetra_angle();

	//Given the polyhedron face number calculate Tetra
	//face and store in array. Used in creating spheres.
	if(classtype == 1){
		if(facenumber >=0 && facenumber <=2){
			//Build upper three faces
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = pntcrd[i].phi + 120.0 * facenumber;
					sphere_pnt[i].theta = pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber == 3){
			//Do bottom tetra face
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].theta = rotate_theta(240.0, pntcrd[i].phi, acos(-1.0/3.0) * RAD_TO_DEG,
							pntcrd[i].theta);
					sphere_pnt[i].phi = rotate_phi(240.0, pntcrd[i].phi, pntcrd[i].theta,
							sphere_pnt[i].theta);

					if(pnt_label[i].A >= A_change)
						sphere_pnt[i].phi = 180.0 - sphere_pnt[i].phi;
					if(sphere_pnt[i].phi < 0.0)
						sphere_pnt[i].phi += 360.0;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
	}
	else if(classtype == 2){
		if(facenumber >=0 && facenumber <=2 ){
			//Build top cap
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = pntcrd[i].phi + 120.0 * facenumber;
					sphere_pnt[i].theta = pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >= 3 && facenumber <= 5){
			rot_factor = facenumber - 3;
			//do a downward rotation then rotate into position.
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].theta = rotate_theta(60.0, pntcrd[i].phi, acos(-1.0/3.0) * RAD_TO_DEG,
							pntcrd[i].theta);
					sphere_pnt[i].phi = (rotate_phi(60.0, pntcrd[i].phi, pntcrd[i].theta,
							sphere_pnt[i].theta) + 60.0) + (120.0 * rot_factor);

					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >= 6 && facenumber <= 8){
			rot_factor = facenumber - 6;
			//do a downward rotation then rotate into position. Flip 180 deg
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].theta = rotate_theta(60.0, pntcrd[i].phi, acos(-1.0/3.0) * RAD_TO_DEG,
							pntcrd[i].theta);
					sphere_pnt[i].phi = (rotate_phi(60.0, pntcrd[i].phi, pntcrd[i].theta,
							sphere_pnt[i].theta) ) + (120.0 * rot_factor);

					sphere_pnt[i].theta = 180.0 - sphere_pnt[i].theta;

					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
		else if(facenumber >=9 && facenumber <=11 ){
			//Build bottom cap
			for(i=1; i<=vertex_calc; i++){
				if((bucky_ball && fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0) || !bucky_ball){
					sphere_pnt[i].phi = (pntcrd[i].phi + 120.0 * facenumber) + 60.0;
					sphere_pnt[i].theta = 180.0 - pntcrd[i].theta;
					if(!ellipse_flag && !ParabolaFlag)
						sphere_pnt[i].radius = 1.0;
					else if(ellipse_flag){
						sphere_pnt[i].theta = root_E(E, sphere_pnt[i].theta);
						sphere_pnt[i].radius = ellipse_radius(E, sphere_pnt[i].theta);
					}
					else
						sphere_pnt[i].radius = parabolic_radius(sphere_pnt[i].theta, ParabolaFocus);
				}
			}
		}
	}
}
//-------------------calculate root E theta correction------------
double CGeodesic::root_E(double E, double theta)
{
	//calculate the root E theta correction to even out elliptical edge length
	//distributions.

	double result, X, Y;

	if(theta == 90.0)
		result = 90.0;
	else if(theta == 270.0)
		result = 270.0;
	else{
		result = clean_float(atan(tan(theta * DEG_TO_RAD) / pow(E, 0.5)) * RAD_TO_DEG);
		if(result >= 0.0){
			X = result;
			Y = 180.0 + result;
			if(fabs(X - theta) <= fabs(Y - theta))
				result = X;
			else
				result = Y;
		}
		else{
			X = 180.0 + result;
			Y = 360.0 + result;
			if(fabs(X - theta) <= fabs(Y - theta))
				result = X;
			else
				result = Y;
		}
	}

	return(result);
}
//-------------------calculate radius of ellipse
double CGeodesic::ellipse_radius(double E, double theta)
{
	//calculate the radius for a given E and theta

	double result;

	result = pow(pow(E, 2.0) / (pow(E, 2.0) * pow(sin(theta * DEG_TO_RAD), 2.0)
					 + pow(cos(theta * DEG_TO_RAD), 2.0)), 0.5);

	return(result);
}
//-------------------Determine Axial Angle A
double CGeodesic::axial_angle_A(long index)
{
	//Calculates the first axial angle for chord "index"
	double sphi, stheta, ephi, etheta, length, result;

	sphi = pntcrd[edgepts[index].start].phi;
	stheta = pntcrd[edgepts[index].start].theta;
	ephi = pntcrd[edgepts[index].end].phi;
	etheta = pntcrd[edgepts[index].end].theta;
	length = chord_length(sphi, stheta, ephi, etheta);

	result = acos((length*length+pntcrd[edgepts[index].start].radius*
				pntcrd[edgepts[index].start].radius-pntcrd[edgepts[index].end].radius*
					pntcrd[edgepts[index].end].radius)/(2*length*
						pntcrd[edgepts[index].start].radius)) * RAD_TO_DEG;

	return(result);
}
//-------------------Determine Axial Angle B
double CGeodesic::axial_angle_B(long index, double angleA)
{
		//calculates second axial angle for chord "index" based on angle A

		double result;

		result = asin((pntcrd[edgepts[index].start].radius*
						sin(angleA * DEG_TO_RAD))/
							pntcrd[edgepts[index].end].radius) * RAD_TO_DEG;

		return(result);
}
//------------------Determine Face Angle A
double CGeodesic::face_angle_A(long index)
{
	//Calculates first face angle for face "index"
	double sphi, stheta, ephi, etheta, result;
	double l1, l2, l3;

		//get face lengths
		//First A to C
		sphi = pntcrd[polyface[index].cornerA].phi;
		stheta = pntcrd[polyface[index].cornerA].theta;
		ephi = pntcrd[polyface[index].cornerC].phi;
		etheta = pntcrd[polyface[index].cornerC].theta;
		l1 = chord_length(sphi, stheta, ephi, etheta);
		//Next A to B
		sphi = pntcrd[polyface[index].cornerA].phi;
		stheta = pntcrd[polyface[index].cornerA].theta;
		ephi = pntcrd[polyface[index].cornerB].phi;
		etheta = pntcrd[polyface[index].cornerB].theta;
		l2 = chord_length(sphi, stheta, ephi, etheta);
		//Next B to C
		sphi = pntcrd[polyface[index].cornerB].phi;
		stheta = pntcrd[polyface[index].cornerB].theta;
		ephi = pntcrd[polyface[index].cornerC].phi;
		etheta = pntcrd[polyface[index].cornerC].theta;
		l3 = chord_length(sphi, stheta, ephi, etheta);
		result = acos((l1*l1+l2*l2-l3*l3)/(2*l1*l2)) * RAD_TO_DEG;

		return(result);
}
//------------------Determine Face Angle B
double CGeodesic::face_angle_B(long index, double faceA)
{
	//Calculates face angle B for face "index" based on face angle A
	double sphi, stheta, ephi, etheta, result;
	double l1, l3;

		//get face lengths
		//First A to C
		sphi = pntcrd[polyface[index].cornerA].phi;
		stheta = pntcrd[polyface[index].cornerA].theta;
		ephi = pntcrd[polyface[index].cornerC].phi;
		etheta = pntcrd[polyface[index].cornerC].theta;
		l1 = chord_length(sphi, stheta, ephi, etheta);
		//...Then B to C
		sphi = pntcrd[polyface[index].cornerB].phi;
		stheta = pntcrd[polyface[index].cornerB].theta;
		ephi = pntcrd[polyface[index].cornerC].phi;
		etheta = pntcrd[polyface[index].cornerC].theta;
		l3 = chord_length(sphi, stheta, ephi, etheta);
		result = asin((l1*sin(faceA * DEG_TO_RAD)) / l3) * RAD_TO_DEG;

		return(result);

}
//------------------Determine Face Angle C
double CGeodesic::face_angle_C(double faceA, double faceB)
{
	//Calculate face angle C from face angles A & B
	double result;

	result = 180.0 - faceA - faceB;

	return(result);
}

//-------------------save data to ASCII file
void CGeodesic::save_ascii(char *filename)
{
	//function to save data to ASCII data file
	//save chord data for symmetry triangle to DAT file
	//filename should include .DAT extension on MS-DOS systems
	//Only data for the symmetry triangle is saved

	double sphi, stheta, ephi, etheta;

	ofstream ASCII(filename);
	//output ASCII report data
	ASCII << "Geodesic Dome Data" << '\n';
	ASCII << "----------------Dome Parameters-------------------------" << '\n';
	ASCII << "Polyhedron Type: ";
	if(polytype == 1)
		ASCII << "Icosahedron" << '\n';
	else if(polytype == 2)
		ASCII << "Octahedron" << '\n';
	else
		ASCII << "Tetrahedron" << '\n';
	ASCII << "Class ";
	if(classtype == 1)
		ASCII << "I" << '\n';
	else
		ASCII << "II" << '\n';
	ASCII << "Frequency: " << frequency << '\n';
	ASCII << "----------------Symmetry Triangle Data------------------" << '\n';
	ASCII << "Vertexia: " << vertex << '\n';
	ASCII << "Edges: " << edges << '\n';
	ASCII << "Faces: " << faces << '\n';
	ASCII << "-----------Symmetry Triangle Spherical Coordinates------" << '\n';
	ASCII << "Vertex #" << '\t' << "phi" << '\t' << "theta" << '\n';

	for(i=1; i<=vertex_calc; i++){
		ASCII << setiosflags(ios::right) << setw(3);
		ASCII << i << '\t';
		ASCII << setw(11) << setprecision(8) << setiosflags(ios::fixed);
		ASCII << pntcrd[i].phi << '\t';
		ASCII << pntcrd[i].theta << '\n';
	}
	ASCII << '\n';
	ASCII << "----------Symmetry Triangle Chord Data------------------" << '\n';
	ASCII << "Vertex" << '\n' << "start/end" << '\t' << "Chord Length" << '\n';
	for(i=1; i<=edges_calc; i++){
		ASCII << setiosflags(ios::right) << setw(3);
		ASCII << setw(12) << setprecision(8) << setiosflags(ios::fixed);
		ASCII << edgepts[i].start << "-" << edgepts[i].end << '\t';
		sphi = pntcrd[edgepts[i].start].phi;
		stheta = pntcrd[edgepts[i].start].theta;
		ephi = pntcrd[edgepts[i].end].phi;
		etheta = pntcrd[edgepts[i].end].theta;
		ASCII << chord_length(sphi, stheta, ephi, etheta) << '\n';
	}
	ASCII << "---------Symmetry Triangle Axial Angle Data-------------" << '\n';
	ASCII << "Vertex" << '\n' << "start/end" << '\t' << "Axial #1" << '\t' << "Axial #2" << '\n';
	for(i=1; i<=edges_calc; i++){
		ASCII << setiosflags(ios::right) << setw(3);
		ASCII << setw(12) << setprecision(8) << setiosflags(ios::fixed);
		ASCII << edgepts[i].start << "-" << edgepts[i].end << '\t';
		ASCII << axial_angle_A(i) << '\t' << axial_angle_B(i, axial_angle_A(i)) << '\n';
	}
	ASCII << "---------Symmetry Triangle Face Angle Data-------------" << '\n';
	ASCII << "Face" << '\t' << "Face Angle #1" << '\t' << "Face Angle #2" << '\t' << "Face Angle #3" << '\n';
	for(i=1; i<=face_calc; i++){
		ASCII << setiosflags(ios::right) << setw(3);
		ASCII << setw(12) << setprecision(8) << setiosflags(ios::fixed);
		ASCII << i << '\t';
		ASCII << face_angle_A(i) << '\t' << face_angle_B(i, face_angle_A(i)) << '\t';
		ASCII << face_angle_C(face_angle_A(i), face_angle_B(i, face_angle_A(i))) << '\n';
	}
	ASCII << "End Dome Data" << '\n';

	ASCII.close();
}

//-------------------calculate chord lengths
double CGeodesic::chord_length(double sphi, double stheta, double ephi, double etheta)
{
	//Function returns distance between two points given their
	//spherical coordinates. Radius is fixed at one.

	double result;

	result = pow((2-2 * (cos(stheta * DEG_TO_RAD) * cos(etheta * DEG_TO_RAD) +
		cos((sphi-ephi) *DEG_TO_RAD) * sin(stheta * DEG_TO_RAD) *
			sin(etheta * DEG_TO_RAD))),0.5);

	return result;

}
//------------------phi rotation function
double CGeodesic::rotate_phi(double phi, double phi1, double theta1, double theta2)
{
		//phi half of the rotation formula given in Kenner's "Geodesic Math &
		//How to Use it" Kenner credits this formula to Professor George Owen

		double result;

		result = (sin(theta1 * DEG_TO_RAD)* sin((phi-phi1) * DEG_TO_RAD)) /
						sin(theta2 * DEG_TO_RAD);

		//Apply correction for round-off errors
		if(result > 1.0)
			result = 1.0;
		else if(result < -1.0)
			result = -1.0;

		return asin(result) * RAD_TO_DEG;

}
//------------------theta rotation function
double CGeodesic::rotate_theta(double phi, double phi1, double theta, double theta1)
{
		//theta half of the rotation formula given in Kenner's "Geodesic Math &
		//How to Use it" Kenner credits this formula to Professor George Owen
		double result;

		result = cos(theta1 * DEG_TO_RAD) * cos(theta * DEG_TO_RAD) +
						sin(theta1 * DEG_TO_RAD) * sin(theta * DEG_TO_RAD) *
							cos((phi-phi1) * DEG_TO_RAD);

		//Apply correction for round-off errors
		if(result > 1.0)
			result = 1.0;
		else if(result < -1.0)
			result = -1.0;

		return acos(result) * RAD_TO_DEG;

}

//-------------Determine "A" Coordinate of bottom tetra face-----------
long CGeodesic::tetra_angle(void)
{

		//This function returns the "A" coordinate of the bottom tetrahedron face
		//where correction of the phi angle must begin. It is the result of quirks
		//in professor owens rotation formulae. I don't find this code very elegant
		//but it does provide the proper results.

		double sX, sphi, stheta;
		double eX, ephi, etheta;

		//Save the 0,0 coordinate as the start point and convert to an X value
		stheta = rotate_theta(240.0, pntcrd[1].phi, acos(-1.0/3.0) * RAD_TO_DEG,
						pntcrd[1].theta);
		sphi = rotate_phi(240.0, pntcrd[1].phi, pntcrd[1].theta, stheta);
		sX = clean_float(cos(sphi * DEG_TO_RAD) * sin(stheta * DEG_TO_RAD));

		//The goal of this routine is to find the "A" coordinate in which the
		//difference between the current X and the last X is positive. We start
		//at 0,0 then proceed down the tetra face edge until this is found
		//(i.e. 0,0 to 1,0; 1,0 to 2,0; 2,0 to 3,0;...f-1,0 to f,0)

		index = 1;
		for(i=2; i<=freq_calc; i++){
			index += (i-1);	//index given the array location of the tetra edge point

			//Save the next coordinate as the end point and convert to an X value
			etheta = rotate_theta(240.0, pntcrd[index].phi, acos(-1.0/3.0) * RAD_TO_DEG,
							pntcrd[index].theta);
			ephi = rotate_phi(240.0, pntcrd[index].phi, pntcrd[index].theta, etheta);

			eX = clean_float(cos(ephi * DEG_TO_RAD) * sin(etheta * DEG_TO_RAD));

			if((eX - sX) > 0.0)
				return pnt_label[index].A;
			else
				sX = eX;
		}

		//Handle special cases of frequency < 5
		if(freq_calc == 4)
			return 3;
		if(freq_calc == 1)
			return 1;
		else
			return 2;
}

//--------------------display data function
void CGeodesic::display_data(void)
{
		//Set the topology dependent on class type
		if(classtype == 1){
			vertex_calc = vertex;
			edges_calc = edges;
		}
		else if(classtype == 2){
			vertex_calc =vertexII;
			edges_calc = edgesII;
		}
		cout << "Geodesic Dome Data" << '\n';
		cout << "----------------Dome Parameters-------------------------" << '\n';
		cout << "Polyhedron Type: ";
		if(polytype == 1)
			cout << "Icosahedron" << '\n';
		else if(polytype == 2)
			cout << "Octahedron" << '\n';
		else
			cout << "Tetrahedron" << '\n';
		cout << "Class ";
		if(classtype == 1)
			cout << "I" << '\n';
		else
			cout << "II" << '\n';
		cout << "Frequency: " << frequency << '\n';

		if(bucky_ball)
			cout << "----------------Buckyball Data--------------------------" << '\n';
		else
			cout << "----------------Symmetry Triangle Data------------------" << '\n';

		if(bucky_ball){
			cout << "Vertexia: " << bucky_vertex << '\n';
			cout << "Edges: " << bucky_edges << '\n';
			cout << "Faces: " << bucky_face << '\n' << '\n';
		}
		else{
			cout << "Vertexia: " << vertex_calc << '\n';
			cout << "Edges: " << edges_calc << '\n';
			cout << "Faces: " << face_calc << '\n' << '\n';
		}

}

//--------------------The Buckyball chord function
void CGeodesic::bucky_factor(void)
{
	//Function cycles through each vertex and determines the chord connections for
	//buckyball constructs.
	index = 1;
	if(!show_status){
		cout << "Calculating Buckyball Chord Factors... ";
		status_count = 0;
	}

	//create array object for edges
	edgepts = new CChords[(edges + 1)];
	if (edgepts == NULL){
		cout << "Insufficient memory for chord array --- Execution Terminating" << '\n';
		exit(1);
	}

	for(i=1; i<vertex_calc; i++){
		if(!show_status){
			time_passage(status_count);
			status_count++;
			if(status_count > 3)
				status_count = 0;
		}
		if(fmod((pnt_label[i].A + pnt_label[i].B), 3.0) != 0.0){

			if(pnt_label[i].A == freq_calc && pnt_label[i].B == 0){
				//point is a left vertex; add 01
				//line #1 start point definition
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A && pnt_label[j].B == pnt_label[i].B + 1)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}
			}
			else if(pnt_label[i].A == pnt_label[i].B){
				//point is a right edge; add 10 & 11
				//line #1 start point definition
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}

				//Line #2 start point
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}
			}
			else if(pnt_label[i].A == freq_calc){
				//point is a bottom vertex; add 01
				//line #1 start point definition
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A && pnt_label[j].B == pnt_label[i].B + 1)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}
			}
			else if(pnt_label[i].B == 0){
				//point is a right edge; add 01, 10 & 11
				//line #1 start point definition
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A && pnt_label[j].B == pnt_label[i].B + 1)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}

				//Line #2 start point
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}
				//Line #3 start point
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}
			}
			else{
				//point is an interior vertex; add 01, 10 & 11
				//line #1 start point definition
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A && pnt_label[j].B == pnt_label[i].B + 1)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}

				//Line #2 start point
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}
				//Line #3 start point
				edgepts[index].start = i;

				for(j=i; j<=vertex_calc; j++){
					if(pnt_label[j].A == pnt_label[i].A + 1 && pnt_label[j].B == pnt_label[i].B + 1)
						break;
				}

				if(fmod((pnt_label[j].A + pnt_label[j].B), 3.0) != 0.0){
					//end point
					edgepts[index].end = j;
					index++;
				}
			}
		}
	}
	cout << '\r' << "                                        " << '\r'; //Clear status signal
}

//--------------------The zenith altitude function
double CGeodesic::zenith(double radius, double theta)
{
	//Calculate the zenith altitude for a given truncation plane specified by theta

	return radius - cos(theta * DEG_TO_RAD);
}

//--------------------The floor radius function
double CGeodesic::floor(double radius, double theta)
{
	//calculate the floor radius for a given truncation defined by theta

	return radius * sin(theta * DEG_TO_RAD);

}

double CGeodesic::parabolic_radius(double theta, double focus)
{
		//calculate the radius for a given focus location and theta

	double result;

	result = (2 * focus)/(1 + cos(theta * DEG_TO_RAD));

	return(result);
}


double CGeodesic::ParabolaTheta(double radius, double focus)
{
	//Calulate the theta angle for a given parabolic opening diameter

	double x;
	x = (radius * radius)/(4*focus);

	return atan(radius/(focus - x)) * RAD_TO_DEG;
}


double CGeodesic::FindTruncLimit(double Trunc)
{
	long end_face, indexinc;
	double peak = 0.0; // AR: set to 0.0 (?) otherwise compared before set


	//determine the number of faces required (first face = 0)
	if(sphere_flg)
		end_face = face_quantity(classtype, polytype);
	else
		end_face = 0;

	if(polytype == 1)
		indexinc = 5;
	else if(polytype == 2)
		indexinc = 4;
	else if(polytype == 3)
		indexinc = 3;
	
	//Allocate Peakdata array space
	Truncation = new CThetaTruncation[end_face + 1];
	if (Truncation == NULL){
		cout << "Insufficient memory for Truncation array --- Execution Terminating" << '\n';
		exit(1);
	}	
	
	for(j=0; j<=end_face; j+=indexinc){
		if(polytype == 1)
			icosa_sphere(j);
		else if(polytype == 2)
			octa_sphere(j);
		else if(polytype == 3)
			tetra_sphere(j);

		for(i=1; i<=vertex_calc; i++){
			//Find peak Theta
	
			if(sphere_pnt[i].theta <= Trunc){
				if(sphere_pnt[i].theta > Truncation[j].InitialPeakTheta){
					Truncation[j].InitialPeakTheta = sphere_pnt[i].theta;
					Truncation[j].InitialPeakA = pnt_label[i].A;
					Truncation[j].InitialPeakB = pnt_label[i].B;
				}	
			}
		}
		
		if(Truncation[j].InitialPeakTheta > peak){
			peak = Truncation[j].InitialPeakTheta;
			peakrow = Truncation[j].InitialPeakA;
		}


		Truncation[j].PeakTheta = Truncation[j].InitialPeakTheta;
		Truncation[j].PeakA = Truncation[j].InitialPeakA;
		Truncation[j].PeakB = Truncation[j].InitialPeakB;
		long LabelIndex = 0;


		LabelIndex++;
		for(i=1; i<=vertex_calc; i++){
			//Find peak Theta
			if((sphere_pnt[i].theta > Trunc) && (Truncation[j].InitialPeakA == pnt_label[i].A + LabelIndex || Truncation[j].InitialPeakA == pnt_label[i].A - LabelIndex)){
				if(sphere_pnt[i].theta <= 90.0){
						Truncation[j].PeakTheta = sphere_pnt[i].theta;
						Truncation[j].PeakA = pnt_label[i].A;
						Truncation[j].PeakB = pnt_label[i].B;
				}
			}	
			if(Truncation[j].PeakTheta > peak){
				peak = Truncation[j].PeakTheta;
				peakrow = Truncation[j].PeakA;
			}


		}
	}  


	return peak;

}



int CGeodesic::TestParabolicIcosa(CCartesian & A, CCartesian & B, double & TruncTheta, double & ThetaLimit)
{
	long crossA = 0; // AR: set to 0 to remove warning that may be used unset
   long crossB = 0; // AR: set to 0 to remove warning that may be used unset

	CCartesian T;
	int skip_value;

	skip_value = 0;
	if(sphere_pnt[edgepts[i].start].theta > ThetaLimit && sphere_pnt[edgepts[i].end].theta > ThetaLimit){
		skip_value = 1;
		goto ExitPoint;
	}
	else if(sphere_pnt[edgepts[i].start].theta <= TruncTheta && sphere_pnt[edgepts[i].end].theta <= TruncTheta){
		A.phi = sphere_pnt[edgepts[i].start].phi;
		B.phi = sphere_pnt[edgepts[i].end].phi;
		A.theta = sphere_pnt[edgepts[i].start].theta;
		B.theta = sphere_pnt[edgepts[i].end].theta;
		A.radius = sphere_pnt[edgepts[i].start].radius;
		B.radius = sphere_pnt[edgepts[i].end].radius;
	}

	else if(sphere_pnt[edgepts[i].start].theta <= TruncTheta && sphere_pnt[edgepts[i].end].theta > TruncTheta){
		A.phi = sphere_pnt[edgepts[i].start].phi;
		B.phi = sphere_pnt[edgepts[i].end].phi;
		A.theta = sphere_pnt[edgepts[i].start].theta;
		B.theta = TruncTheta;
		A.radius = sphere_pnt[edgepts[i].start].radius;
		B.radius = parabolic_radius(B.theta, ParabolaFocus);
	}
	else if(sphere_pnt[edgepts[i].start].theta > TruncTheta && sphere_pnt[edgepts[i].end].theta <= TruncTheta){
		A.phi = sphere_pnt[edgepts[i].start].phi;
		B.phi = sphere_pnt[edgepts[i].end].phi;
		A.theta = TruncTheta;
		B.theta = sphere_pnt[edgepts[i].end].theta;
		A.radius = parabolic_radius(A.theta, ParabolaFocus);
		B.radius = sphere_pnt[edgepts[i].end].radius;
	}
	else if(sphere_pnt[edgepts[i].start].theta > TruncTheta && sphere_pnt[edgepts[i].end].theta > TruncTheta){
		if(pnt_label[edgepts[i].start].A != pnt_label[edgepts[i].end].A){
			if(pnt_label[edgepts[i].start].A <= peakrow -1 || pnt_label[edgepts[i].end].A <= peakrow -1){
				skip_value = 1;
				goto ExitPoint;
			}

			//Set end point thetas to TruncTheta
			A.phi = sphere_pnt[edgepts[i].start].phi;
			B.phi = sphere_pnt[edgepts[i].end].phi;
			A.theta = TruncTheta;
			B.theta = TruncTheta;
			A.radius = parabolic_radius(A.theta, ParabolaFocus);
			B.radius = parabolic_radius(B.theta, ParabolaFocus);

							
			A.SphericalToCartesian();
			B.SphericalToCartesian();

		}
		else if( ((j>=5 && j<=9) && pnt_label[edgepts[i].end].A == peakrow + 1) ||
               ((j>=10 && j<=14) && pnt_label[edgepts[i].end].A == (frequency - peakrow) -1)) {
			//Set end point thetas to TruncTheta
			A.phi = sphere_pnt[edgepts[i].start].phi;
			B.phi = sphere_pnt[edgepts[i].end].phi;
			A.theta = TruncTheta;
			B.theta = TruncTheta;
			A.radius = parabolic_radius(A.theta, ParabolaFocus);
			B.radius = parabolic_radius(B.theta, ParabolaFocus);
		}
		else if(pnt_label[edgepts[i].start].A == pnt_label[edgepts[i].end].A){
			if(((j>=5 && j<=9) && pnt_label[edgepts[i].end].A == peakrow + 1) ||
            ((j>=10 && j<=14) && pnt_label[edgepts[i].end].A == (frequency - peakrow) -1)) {
				//Set end point thetas to TruncTheta
				A.phi = sphere_pnt[edgepts[i].start].phi;
				B.phi = sphere_pnt[edgepts[i].end].phi;
				A.theta = TruncTheta;
				B.theta = TruncTheta;
				A.radius = parabolic_radius(A.theta, ParabolaFocus);
				B.radius = parabolic_radius(B.theta, ParabolaFocus);
			}
			//Now check what type of triangle and verify that the points are on peakrow.				
			else if(((j>=5 && j<=9) && pnt_label[edgepts[i].end].A == peakrow) ||
                 ((j>=10 && j<=14) && pnt_label[edgepts[i].end].A == frequency - peakrow)) {
				A.phi = sphere_pnt[edgepts[i].start].phi;
				B.phi = sphere_pnt[edgepts[i].end].phi;
				A.theta = TruncTheta;
				B.theta = TruncTheta;
				A.radius = parabolic_radius(A.theta, ParabolaFocus);
				B.radius = parabolic_radius(B.theta, ParabolaFocus);					

				//convert endpoints to cartesian
				A.SphericalToCartesian();
				B.SphericalToCartesian();

				//find location of cross point
				if(j>=5 && j<=9){
					//Look for which endpoint has largest B label
					if(pnt_label[edgepts[i].start].B > pnt_label[edgepts[i].end].B)
						crossB = pnt_label[edgepts[i].start].B;
					else
						crossB = pnt_label[edgepts[i].end].B;

					//Add one to the start label
					crossA = pnt_label[edgepts[i].start].A + 1;
				}
				else if(j>=10 && j<=14){
					//Look for which endpoint has largest B label
					//Subtract one from it.
					if(pnt_label[edgepts[i].start].B > pnt_label[edgepts[i].end].B)
						crossB = pnt_label[edgepts[i].start].B - 1;
					else
						crossB = pnt_label[edgepts[i].end].B - 1;
					
					//Subtract one from the start label
					crossA = pnt_label[edgepts[i].start].A - 1;
				}

				//Search for cross point; exit loop when found. need the k value
				for(k=1; k<=vertex_calc; k++){
					if(pnt_label[k].A == crossA && pnt_label[k].B == crossB)
						break;
				}

				//Check that the check point is within truncation limit if not adjust
					if(sphere_pnt[k].theta > TruncTheta){
						T.theta = TruncTheta;
						T.radius = parabolic_radius(T.theta, ParabolaFocus);
						T.phi = sphere_pnt[k].phi;
					}
					else{
 						T.theta = sphere_pnt[k].theta;
						T.radius = sphere_pnt[k].radius;
						T.phi = sphere_pnt[k].phi;
					}

					T.SphericalToCartesian();

					if(TestChordCrossing(A, B, T) == 1){
						skip_value = 1;
						goto ExitPoint;
					}
				}
				else {
					skip_value = 1;
					goto ExitPoint;
				}
			}
			else{
				skip_value = 1;
				goto ExitPoint;
			}
		}
		else{
			skip_value = 1;
			goto ExitPoint;
		}
		
		//convert to cartesian

		A.SphericalToCartesian();
		B.SphericalToCartesian();

ExitPoint:
		return skip_value;

}


int CGeodesic::TestChordCrossing(CCartesian & A, CCartesian & B, CCartesian & T)
{
	//Used to test whether a truncated chord is crossing an inside vertex.
	double slope, TestY, DeltaY;

	int test_flag;

	test_flag = 0;

	slope = (B.Y - A.Y)/(B.X - A.X);					
	TestY = slope * (T.X - A.X) + A.Y;
		
	//Calculate the difference in Y
	DeltaY = TestY - T.Y;

	if(T.Y > 0.0){
		if(DeltaY <= 0.0){
			test_flag = 1;
		}
	}		
	else{
		if(DeltaY >= 0.0){
			test_flag = 1;
		}
	}

	return test_flag;
}
//--------------------End of class Geodesic
