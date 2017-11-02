#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

/*
    First param: Name of .xyz file you want to open
    Second param: Atom number #1 (corresponds to GaussView)
    Third param: Atom numer #2
    Fourth param: Not included yet, but possible multiply times timestep, if you wanted seconds instead???
*/

#define SHELLSCRIPT "\
#/bin/bash \n\
gnuplot GNUSHELLSCRIPT \n\
rm GNUSHELLSCRIPT\n\
"

void createGNUSS()
{
    ofstream outFile("GNUSHELLSCRIPT");
    outFile <<
"set terminal png \n\
set datafile separator\",\" \n\
set output \"BondLengthGraph.png\" \n\
set ylabel \"Bond Length (Å)\" \n\
set xlabel \"Timesteps\" \n\
set title \"Bond Length\" \n\
plot 'BondLength.csv' with lines\
";
    outFile.close();
}

void usage()
{
    cout << "usage: BondLengthTool XYZ_FILE ATOM_1_INDEX ATOM_2_INDEX" << endl;
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


map<int, vector<float> > getAtomCoordMap(int atom_num, ifstream& xyz_file, int first_wanted_atom, int second_wanted_atom)
{
    map<int, vector<float> > atoms_coords;

    for (int i = 0; i <= atom_num; i++)
    {
        if (i == first_wanted_atom)
        {
            atoms_coords[first_wanted_atom] = getAtomCoord(xyz_file);
        }
        else if (i == second_wanted_atom)
        {
            atoms_coords[second_wanted_atom] = getAtomCoord(xyz_file);
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
    float coord_root = sqrt(sum);
    //cout << "Square Root: " << coord_root << endl;
    return coord_root;
}

void exportBondLengthFile(vector<float> bond_lengths)
{
    ofstream outFile("BondLength.csv");
    for (int i = 0; i < bond_lengths.size(); i++)
    {
        cout << bond_lengths[i] << endl;
        //cout << bond_lengths[i] << endl;
    }
    outFile.close();
}

int main(int argc, char* argv[])
{
    if (argc <= 3)
    {
        usage();
        return 1;
    }

    ifstream xyz_file;
    xyz_file.open(argv[1]);
    string str_first_wanted_atom = argv[2];
    string str_second_wanted_atom = argv[3];
    int first_wanted_atom = atoi(str_first_wanted_atom.c_str());
    int second_wanted_atom = atoi(str_second_wanted_atom.c_str());
    if (first_wanted_atom == second_wanted_atom)
    {
        cout << "ERROR: The two atoms you specified are the same. Please enter two different atoms." << endl;
    }

    vector<float> bond_lengths;
    while(xyz_file.peek() != EOF)
    {
        string atom_num_str;
        getline(xyz_file, atom_num_str);
        int atom_num = atoi(atom_num_str.c_str());
        if (first_wanted_atom > atom_num || second_wanted_atom > atom_num)
        {
            cout << "ERROR: The atoms you specified are out of bounds. Please enter two different atoms." << endl;
            return 1;
        }
        else if (first_wanted_atom < 1 || second_wanted_atom < 1)
        {
            cout << "ERROR: The atoms you specified are out of bounds. Please enter two different atoms." << endl;
            return 1;
        }

        map<int, vector<float> > atom_coord = getAtomCoordMap(atom_num, xyz_file, first_wanted_atom, second_wanted_atom);
        float bond_length = calcBondLength(atom_coord, first_wanted_atom, second_wanted_atom);
        bond_lengths.push_back(bond_length);
    }

    exportBondLengthFile(bond_lengths);
    //createGNUSS();
    //system(SHELLSCRIPT);
    return 0;
}
