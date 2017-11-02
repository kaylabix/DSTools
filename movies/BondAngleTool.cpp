#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdio.h>

using namespace std;

#define SHELLSCRIPT "\
#/bin/bash \n\
gnuplot GNUSCHELLSCRIPT \n\
rm GNUSCHELLSCRIPT\n\
"

void createGNUSS()
{
    ofstream outFile("GNUSCHELLSCRIPT");
    outFile <<
    "set terminal png \n\
    set datafile separator\",\" \n\
    set output \"BondAngleGraph.png\" \n\
    plot 'BondAngle.csv' with lines\
    ";
}

void usage()
{
    cout << "usage: BondAngleTool [XYZ_FILE] [ATOM_1_INDEX] [MIDDLE_ATOM_INDEX] [ATOM_3_INDEX]" << endl;
}

vector<float> getAtomCoord(ifstream& xyz_file)
{
    vector<float> atom_coord;
    string str_atom_coord;
    getline(xyz_file, str_atom_coord);
    istringstream iss(str_atom_coord);
    string atomic_symbol;
    iss >> atomic_symbol;

    for(int i = 0; i < 3 ; i++)
    {
        float coordinate;
        iss >> coordinate;
        atom_coord.push_back(coordinate);
    }

    return atom_coord;
}

map<int, vector<float> > getAtomCoordMap(int atom_num, ifstream& xyz_file, int first_wanted_atom, int second_wanted_atom, int third_wanted_atom)
{
    map<int, vector<float> > atoms_coords;

    for(int i = 0; i <= atom_num; i++)
    {
        if(i == first_wanted_atom)
        {
            atoms_coords[first_wanted_atom] = getAtomCoord(xyz_file);
        }

        else if(i == second_wanted_atom)
        {
            atoms_coords[second_wanted_atom] = getAtomCoord(xyz_file);
        }

        else if(i == third_wanted_atom)
        {
            atoms_coords[third_wanted_atom] = getAtomCoord(xyz_file);
        }

        else
        {
            string ignore;
            getline(xyz_file, ignore);
        }
    }

    return atoms_coords;
}

float calcBondLength(map<int, vector<float> > atoms_coords, int first_wanted_atom, int second_wanted_atom)
{
    float sum = 0;

    for (int i = 0; i < atoms_coords[first_wanted_atom].size(); i++)
    {
        float first_atom_coord = atoms_coords[first_wanted_atom][i];
        //cout << "First: " << first_atom_coord << endl;
        float second_atom_coord = atoms_coords[second_wanted_atom][i];
        //cout << "Second: " << second_atom_coord << endl;
        float coord_diff = second_atom_coord - first_atom_coord;
        //cout << "Difference: " << coord_diff << endl;
        float coord_squared = pow(coord_diff, 2);
        //cout << "Square: " << coord_squared << endl;
        sum = coord_squared + sum;
        //cout << "Sum" << sum << endl;
    }
    float bond_length = sqrt(sum);
    //cout << "Square Root: " << coord_root << endl;
    return bond_length;
}

float calcBondAngle(map<int, vector<float> > atoms_coords, int first_wanted_atom, int second_wanted_atom, int third_wanted_atom)
{
    float a = calcBondLength(atoms_coords, first_wanted_atom, second_wanted_atom);
    float b = calcBondLength(atoms_coords, second_wanted_atom, third_wanted_atom);
    float c = calcBondLength(atoms_coords, first_wanted_atom, third_wanted_atom);
    float PI = 3.14159265;
    float to_degrees = 180/PI;

    float numerator = pow(c, 2) - pow(a, 2) - pow(b,2);
    float denominator = (-2) * (a) * (b);
    float fraction = numerator/denominator;
    float bond_angle = acos(fraction) * to_degrees;
    return bond_angle;
}

void exportBondAngleFile(vector<float> bond_angles)
{
    ofstream outFile("BondAngle.csv");
    for (int i = 0; i < bond_angles.size(); i++)
    {
        cout << i << "," << bond_angles[i] << endl;
    }
}

int main(int argc, char* argv[])
{
    if(argc <= 4)
    {
        usage();
        return 1;
    }
    ifstream xyz_file;
    xyz_file.open(argv[1]);
    string str_first_wanted_atom = argv[2];
    string str_second_wanted_atom = argv[3];
    string str_third_wanted_atom = argv[4];
    int first_wanted_atom = atoi(str_first_wanted_atom.c_str());
    int second_wanted_atom = atoi(str_second_wanted_atom.c_str());
    int third_wanted_atom = atoi(str_third_wanted_atom.c_str());
    //check to make sure none of these are the same
    vector<float> bond_angles;
    while(xyz_file.peek() != EOF)
    {
        string atom_num_string;
        getline(xyz_file, atom_num_string);
        int atom_num = atoi(atom_num_string.c_str());
        if(first_wanted_atom > atom_num || second_wanted_atom > atom_num || third_wanted_atom > atom_num)
        {
            cout << "ERROR: The atoms you specified are out of bounds.  Please enter two different atoms" << endl;
            return 1;
        }
        else if (first_wanted_atom < 1 || second_wanted_atom < 1 || third_wanted_atom < 1)
        {
            cout << "Error the atoms you sepcified are out of bounds. Please enter two different atoms" << endl;
            return 1;
        }

        map<int, vector<float> > atom_coord = getAtomCoordMap(atom_num, xyz_file, first_wanted_atom, second_wanted_atom, third_wanted_atom);
        float bond_angle = calcBondAngle(atom_coord, first_wanted_atom, second_wanted_atom, third_wanted_atom);
        bond_angles.push_back(bond_angle);
    }

    exportBondAngleFile(bond_angles);
    //createGNUSS();
    //system(SHELLSCRIPT);
    return 0;
}



























