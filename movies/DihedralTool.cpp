#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

#define SHELLSCRIPT "\
#/bin/bash \n\
gnuplot GNUSHELLSCRIPT \n\
rm GNUSHELLSCRIPT \n\
"

void createGNUSS()
{
    ofstream outFile("GNUSHELLSCRIPT");
    outFile <<
"set terminal png \n\
set datafile separator\",\" \n\
set output \"DihedralAngleGraph.png\" \n\
set ylabel \"Dihedral Angle (°)\" \n\
set xlabel \"Timesteps\" \n\
set title \"Dihedral Angles\" \n\
plot 'DihedralAngles.csv' with lines\
";
    outFile.close();
}

void usage()
{
    cout << "usage: DihedralAngleTool XYZ_FILE ATOM_1_INDEX ATOM_2_INDEX ATOM_3_INDEX ATOM_4_INDEX" << endl;
}

vector<float> getAtomCoord(ifstream& xyz_file)
{
    vector<float> atom_coord;
    string str_atom_coord;
    getline(xyz_file, str_atom_coord);
    istringstream iss(str_atom_coord);
    string atomic_symbol;
    iss >> atomic_symbol;

    for(int i = 0; i < 3; i++)
    {
        float coordinate;
        iss >> coordinate;
        atom_coord.push_back(coordinate);
    }

    return atom_coord;
}

map<int, vector<float> > getAtomCoordMap(int atom_num, ifstream& xyz_file, int first_wanted_atom, int second_wanted_atom, int third_wanted_atom, int fourth_wanted_atom)
{
    map<int, vector<float> > atoms_coords;

    for (int i = 0; i <= atom_num; i++)
    {
        if(i == first_wanted_atom)
        {
            atoms_coords[first_wanted_atom] = getAtomCoord(xyz_file);
        }

        else if(i == second_wanted_atom)
        {
            atoms_coords[second_wanted_atom] = getAtomCoord(xyz_file);
        }

        else if (i == third_wanted_atom)
        {
            atoms_coords[third_wanted_atom] = getAtomCoord(xyz_file);
        }

        else if (i == fourth_wanted_atom)
        {
            atoms_coords[fourth_wanted_atom] = getAtomCoord(xyz_file);
        }

        else
        {
            string ignore;
            getline(xyz_file, ignore);
        }
    }

    return atoms_coords;
}

void exportDihedralFile(vector<float> dihedral_angles)
{
    ofstream outFile("DihedralAngles.csv");
    for(int i = 0; i < dihedral_angles.size(); i++)
    {
        cout << dihedral_angles[i] << endl;
    }

    outFile.close();
}

vector<float> subtractVectors(vector<float> point_1, vector<float> point_2)
{
    vector<float> resulting_vector;
    float x_coord = point_2[0] - point_1[0];
    float y_coord = point_2[1] - point_1[1];
    float z_coord = point_2[2] - point_1[2];
    resulting_vector.push_back(x_coord);
    resulting_vector.push_back(y_coord);
    resulting_vector.push_back(z_coord);
    return resulting_vector;
}

vector<float> crossVectors(vector<float> point_1, vector<float> point_2)
{
    vector <float> resulting_vector;
    float x_coord = (point_1[1] * point_2[2]) - (point_1[2] * point_2[1]);
    float y_coord = (point_1[2] * point_2[0]) - (point_1[0] * point_2[2]);
    float z_coord = (point_1[0] * point_2[1]) - (point_1[1] * point_2[0]);
    resulting_vector.push_back(x_coord);
    resulting_vector.push_back(y_coord);
    resulting_vector.push_back(z_coord);
    return resulting_vector;
}

float dotVectors(vector<float> point_1, vector<float> point_2)
{
    float x_coord = point_1[0] * point_2[0];
    float y_coord = point_1[1] * point_2[1];
    float z_coord = point_1[2] * point_2[2];
    return x_coord + y_coord + z_coord;
}

float normsProdVectors(vector<float> point_1, vector<float> point_2)
{
    float sq_pt_1 = pow(point_1[0],2) + pow(point_1[1],2) + pow(point_1[2],2);
    float norm_pt_1 = sqrt(sq_pt_1);
    float sq_pt_2 = pow(point_2[0],2) + pow(point_2[1],2) + pow(point_2[2],2);
    float norm_pt_2 = sqrt(sq_pt_2);
    return norm_pt_1 * norm_pt_2;
}

void printVec(vector<float> vector)
{
    for (auto k:vector)
    {
        cout << k << " ";
    }

    cout << endl;
}

float calcDihedralAngle(map<int, vector<float> > atoms_coords, int first_wanted_atom, int second_wanted_atom, int third_wanted_atom, int fourth_wanted_atom)
{
    float deg_conversion = 180/3.141592654;
    vector <float> atom_1_coords = atoms_coords[first_wanted_atom];
    vector <float> atom_2_coords = atoms_coords[second_wanted_atom];
    vector <float> atom_3_coords = atoms_coords[third_wanted_atom];
    vector <float> atom_4_coords = atoms_coords[fourth_wanted_atom];
    vector <float> sub_res_1 = subtractVectors(atom_1_coords, atom_2_coords);
    vector <float> sub_res_2 = subtractVectors(atom_2_coords, atom_3_coords);
    vector <float> sub_res_3 = subtractVectors(atom_3_coords, atom_4_coords);
    vector <float> cross_res_1 = crossVectors(sub_res_1, sub_res_2);
    vector <float> cross_res_2 = crossVectors(sub_res_2, sub_res_3);
    float dot_res = dotVectors(cross_res_1, cross_res_2);
    float norms_prod = normsProdVectors(cross_res_1, cross_res_2);
    float dihedral_angle = acos(dot_res/norms_prod) * deg_conversion;
    //cout << dihedral_angle << endl;
    return dihedral_angle;
}

int main(int argc, char* argv[])
{
    if (argc <= 5)
    {
        usage();
        return 1;
    }

    ifstream xyz_file;
    xyz_file.open(argv[1]);
    string str_first_wanted_atom = argv[2];
    string str_second_wanted_atom = argv[3];
    string str_third_wanted_atom = argv[4];
    string str_fourth_wanted_atom = argv[5];
    int first_wanted_atom = atoi(str_first_wanted_atom.c_str());
    int second_wanted_atom = atoi(str_second_wanted_atom.c_str());
    int third_wanted_atom = atoi(str_third_wanted_atom.c_str());
    int fourth_wanted_atom = atoi(str_fourth_wanted_atom.c_str());
    //make sure none of them equal eachother
    vector<float> dihedral_angles;
    while(xyz_file.peek() != EOF)
    {
        string atom_num_str;
        getline(xyz_file, atom_num_str);
        int atom_num = atoi(atom_num_str.c_str());
        if(first_wanted_atom > atom_num || second_wanted_atom > atom_num || third_wanted_atom > atom_num || fourth_wanted_atom > atom_num)
        {
            cout << "ERROR: The atoms you specified are out of bounds. Please enter different atoms." << endl;
            return 1;
        }
        map<int, vector<float> > atom_coord = getAtomCoordMap(atom_num, xyz_file, first_wanted_atom, second_wanted_atom, third_wanted_atom, fourth_wanted_atom);
        float dihedral_angle = calcDihedralAngle(atom_coord, first_wanted_atom, second_wanted_atom, third_wanted_atom, fourth_wanted_atom);
        dihedral_angles.push_back(dihedral_angle);
    }

     exportDihedralFile(dihedral_angles);
    // createGNUSS();
    // system(SHELLSCRIPT);
    return 0;
}











