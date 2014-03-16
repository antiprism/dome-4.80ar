//$Header:$

/*
DOME.CPP - A program for the calculation of geodesic dome properties
			  using geodesic class

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
	My wife & daughters
	Chris Fearnley
	Kirby Urner
		&
	R. Buckminster Fuller for changing the way I view Universe.


*/

#include "Command.h"
#include "Geodesic.h"

using std::cout;

CCommand	cmd;		//command object

//Function prototypes
void logo_display(void);
void help_display(void);
void get_cmd(CCommand &command, int param_count, char *param[]);
int main(int argc, char *argv[]);

void help_display(void)
{
	//Display usage and help then exit
	logo_display();
	cout << "Usage: dome [-fnnn] [-cn] [-px] [-s] [-sb] [-en] [-v] [-w] [-h] [filename.xxx]" << '\n';
	cout << "Where: -fnnn is geodesic frequency (default nnn=3)" << '\n';
	cout << "       -cn is class type (n=1 or 2; default n=1)" << '\n';
	cout << "       -en enables ellipse and specifies eccentricity (default = 1.0)" << '\n';
	cout << "           n > 0.0 and n < 2.0" << '\n';
	cout << "       -px sets the polyhedron type" << '\n';
	cout << "           where x is: i for icosahedron (default)" << '\n';
	cout << "                       o for octahedron" << '\n';
	cout << "                       t for tetrahedron" << '\n';
	cout << "       -s  generate full sphere data (default: symmetry triangle)" << '\n';
	cout << "       -sb generate buckyball. Sets Class I" << '\n';
	cout << "           Frequency must be a multiple of three" << '\n';
	cout << "       -w  Enable wire-frame DXF or VRML output (default: Face data)" << '\n';
	cout << "       -v  verbose data display at run-time" << '\n';
	cout << "       -dn enables parabolid and specifies focus location" << '\n';
	cout << "       -rn sets Outer radius of parabolid." << '\n';
	cout << "       -h  displays this help screen" << '\n';
	cout << "           filename.xxx is a standard DOS filename" << '\n';
	cout << "           where xxx is: DXF, WRL, DAT, POV, PRN or OFF" << '\n';
	exit(2);
}
void logo_display(void)
{
   /* AR: unused code
	char rev[]="$Revision: 17a $";
	int i, j;
	char level[10];

	//Get the revision level. I have been using GNU RCS to do revision control
	//on dome. This revision level will be used to automatically display the
	//correct rev level.
	j = 11;
	i = 0;
	while(rev[j] != ' '){
		level[i] = rev[j];
		j++;
		i++;
	}
	level[i] = '\0';
   */

	cout << "Dome 4.80, Copyright (C) 1995-2002 - Richard J. Bono" << '\n';
	cout << "Dome comes with ABSOLUTELY NO WARRANTY. This is free software," << '\n';
	cout << "and you are welcome to redistribute it under certain conditions." << '\n';
	cout << "See GNU General Public License for more details." << '\n' << '\n';
}

void get_cmd(CCommand &command, int param_count, char *param[])
{
	//Get and parse command line
	char cmd_parm[6];
	int t, j, k;

	for(t=1; t<param_count; ++t){
		if(param[t][0] == '-'){
			switch (tolower(param[t][1])){
				case ('p'):
					//Set polyhedron type
					if(tolower(param[t][2]) == 'i')
						//Set to Icosahedron
						command.polyt = 1;
					else if(tolower(param[t][2]) == 'o')
						//Set to octahedron
						command.polyt = 2;
					else if(tolower(param[t][2]) == 't')
						//Set to tetrahedron
						command.polyt = 3;
					else{
						cout << "Invalid Polyhedron Type --- Execution Terminating" << '\n';
						exit(1);
					}
					break;
				case ('h'):
					//Display help and exit. This overrides other parameters
					help_display();
					break;
				case ('v'):
					//Enable Parameter Display during execution
					command.verbose_flag = 1;
					break;
				case ('w'):
					//Enable wire frame output
					command.faceflag = 0;
					break;
				case ('s'):
					//generate buckyball or sphere
					if(tolower(param[t][2]) == 'b'){
						command.buckyball = 1;
						command.sphere_flag = 1;
					}
					else if(param[t][2] == '\0'){
						command.sphere_flag = 1;
						command.buckyball = 0;
					}
					else{
						cout << "Invalid command-line --- Execution Terminating" << '\n';
						exit(1);
					}
					break;
				case ('f'):
					//get frequency
					j=2;
					k=0;
					while(param[t][j]){
						cmd_parm[k] = param[t][j];
						j++;
						k++;
					}
					cmd_parm[k] = '\0';
					command.freq = atol(cmd_parm);
					if (command.freq <= 0){
						cout << "Invalid Frequency --- Execution Terminating" << '\n';
						exit(1);
					}
					break;
				case ('d'):
					//get focus location
					j=2;
					k=0;
					while(param[t][j]){
						cmd_parm[k] = param[t][j];
						j++;
						k++;
					}
					cmd_parm[k] = '\0';
					command.ParabolicFocus = atof(cmd_parm);
					if (command.ParabolicFocus <= 0.0){
						cout << "Invalid Focus - Execution Terminating" << '\n';
						exit(1);
					}
					command.ParabolicFlag = 1; //Enable paraboloid
					command.sphere_flag = 1;   //enable sphere	
					command.buckyball = 0;	   //disable buckyball

					break;
				case ('r'):
					//get parabolic radius
					j=2;
					k=0;
					while(param[t][j]){
						cmd_parm[k] = param[t][j];
						j++;
						k++;
					}
					cmd_parm[k] = '\0';
					command.ParabolicRadius = atof(cmd_parm);
					if (command.ParabolicRadius <= 0.0){
						cout << "Invalid Parabolic Diameter - Execution Terminating" << '\n';
						exit(1);
					}
					break;
				case ('e'):
					//get ellipse flag and eccentricity
					j=2;
					k=0;
					while(param[t][j]){
						cmd_parm[k] = param[t][j];
						j++;
						k++;
					}
					cmd_parm[k] = '\0';
					command.E = atof(cmd_parm);
					if (command.E <= 0.0 || command.E >= 2.0){
						cout << "Invalid Eccentricity --- Execution Terminating" << '\n';
						exit(1);
					}
					break;
				case ('c'):
					//get class type
					j=2;
					k=0;
					while(param[t][j]){
						cmd_parm[k] = param[t][j];
						j++;
						k++;
					}
					cmd_parm[k] = '\0';
					command.classt = atol(cmd_parm);
					if (command.classt < 1 || command.classt > 2){
						cout << "Invalid Class Type --- Execution Terminating" << '\n';
						exit(1);
					}
					break;
				default:
					cout << "Invalid command-line --- Execution Terminating" << '\n';
					exit(1);
					break;
			}
		}
		else{
			//Check to see if parameter is a valid filename
			j = 0;
			while(param[t][j]){
				if(param[t][j] != '.')
					j++;
				else{
					//Check for valid extension
					j++;
					k = 0;
					while(param[t][j]){
						cmd_parm[k] = tolower(param[t][j]);
						j++;
						k++;
					}
					cmd_parm[k] = '\0';
					if(!strcmp(cmd_parm,"dxf"))
						command.filet = 1;
					else if(!strcmp(cmd_parm, "dat"))
						command.filet = 2;
					else if(!strcmp(cmd_parm, "pov"))
						command.filet = 3;
					else if(!strcmp(cmd_parm, "prn"))
						command.filet = 4;
					else if(!strcmp(cmd_parm, "wrl"))
						command.filet = 5;
					else if(!strcmp(cmd_parm, "off"))
						command.filet = 6;
					else{
						cout << "Invalid command-line --- Execution terminating" << '\n';
						exit(1);
					}
					j = 0;
					k = 0;
					while(param[t][j]){
						command.filename[k] = tolower(param[t][j]);
						j++;
						k++;
					}
					command.filename[k] = '\0';
				}
			}
		}
	}

	//Check input for Buckyballs
	//Force class I. Stop if frequency is not multiple of three
	if(command.buckyball){
		if(fmod(command.freq, 3.0) != 0.0){
			cout << "Buckyball Frequency must be a Multiple of Three --- Execution Terminating" << '\n';
			exit(1);
		}
		//Override polyhedron and class values if needed...
		if(command.classt != 1){
			command.classt = 1;
			cout << "Class I set for Buckyball" << '\n';
		}
	}

	if(fmod(command.freq, 2.0) != 0.0 && command.classt == 2){
		cout << "Class II requires Even Frequency --- Execution Terminating" << '\n';
		exit(1);
	}
}

int main(int argc, char *argv[])
{
	//This main routines shows what can be done with this class. It implements
	//a simple program which displays dome parameters and optionally saves the
	//data in either DXF or ASCII file formats.

	//Some things that can be tried involve creating several instances of the
	//geodesic object. Geodesic shells can be created in this way given enough
	//memory.

	//Command line parameters are changing this routine into something more
	//permanent. Remember the class can be incorporated into your own programs.
	//This shell is just one possible implementation.

	//Get dome parameters
	if(argc == 1)
		//display usage and exit
		help_display();
	else
		//Get command line
		get_cmd(cmd, argc, argv);

	logo_display();

	
	//class instance
	CGeodesic geosys(cmd);
	//start calculations
	geosys.topology();			//calculate symmetry triangle topological abundance

	geosys.spherical(); 		//calculate spherical coordinates
	if(geosys.bucky_ball)
		geosys.bucky_factor();	//Determine buckyball edge factors
	else{
		geosys.chord_factor();	//determine chord factors
		geosys.face_factor();	//determine face factors
	}


	if(cmd.verbose_flag)
		geosys.display_data();

	if(cmd.filet == 1){
		//output file in DXF format
		if(cmd.buckyball){
			geosys.save_buckydxf(cmd.filename);
		}
		else if(cmd.faceflag)
			geosys.save_dxf(cmd.filename);
		else
			geosys.save_dxf_wire(cmd.filename);
	}
	else if(cmd.filet == 2){
		//output all data in ASCII report format
		if(cmd.buckyball){
			cout << "ASCII Report not Valid for Buckyball Formation" << '\n';
			exit(1);
		}
		else if(cmd.E != 1.0){
			cout << "ASCII Report not Valid for Elliptical Formation" << '\n';
			exit(1);
		}
		else
			geosys.save_ascii(cmd.filename);
	}
	else if(cmd.filet == 3){
		//output data in POV-Ray script format
		if(cmd.buckyball)
			geosys.save_buckypov(cmd.filename);
		else
			geosys.save_POV(cmd.filename);
	}
	else if(cmd.filet == 4){
		//output data in PRN Text format
		if(cmd.buckyball){
			cout << "PRN Data File not Valid for Buckyball Formation" << '\n';
			exit(1);
		}
		else if(cmd.E != 1.0){
			cout << "PRN Data File not Valid for Elliptical Formation" << '\n';
			exit(1);
		}
		else
			geosys.save_PRN(cmd.filename);
	}
	else if(cmd.filet == 5){
		//output data in VRML format
		if(cmd.buckyball){
			geosys.save_buckywrl_2(cmd.filename);
		}
		else if(cmd.faceflag)
			geosys.save_WRL_2(cmd.filename);
		else
			geosys.save_WRL_wire_2(cmd.filename);
	}
	if(cmd.filet == 6){
		//output file in OFF format
		if(cmd.buckyball){
         geosys.save_buckyoff(cmd.filename);
      }
      else
         geosys.save_OFF(cmd.filename);
	}

	if(cmd.filet == 0)
		cout << "Execution Complete" << '\n';
	else{
		cout << '\r' << "                      " << '\r';
		cout << "Execution Complete -- Output File: " << cmd.filename << '\n';
	}

	return(0);

}
