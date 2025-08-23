#include <cstdlib>
#include <map>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "JUNO_PMTs.h"

using namespace std;

void Read_file (string FileName, std::map <int,tuple<int,double,double,double,double,double>>& out_map, bool is_Wp) {

    ifstream PMT_File(FileName);

    if (!PMT_File.is_open()) {
        cerr << "Error: Could not open file.\n";
        exit(1);
    }

    string line,value;

    for (int i=0;i<4;i++) getline (PMT_File,line); //Four header lines in file

    while (getline(PMT_File,line)) {

        stringstream ss(line);
        vector <string> separated_values;

        while (getline(ss,value,',')) {
            separated_values.push_back(value);
        }

        tuple <int,double,double,double,double,double> PMT_Pos;

        if (is_Wp) PMT_Pos = make_tuple(stoi(separated_values[6]),stod(separated_values[10]),stod(separated_values[11]),stod(separated_values[12]),stod(separated_values[13]),stod(separated_values[14]));
        else PMT_Pos = make_tuple(stoi(separated_values[6]),stod(separated_values[9]),stod(separated_values[10]),stod(separated_values[11]),stod(separated_values[12]),stod(separated_values[13]));

        out_map[stoi(separated_values[0])] = PMT_Pos;

        separated_values.clear();

    }

    PMT_File.close();


};

JUNO_PMTs::JUNO_PMTs () {};

JUNO_PMTs::JUNO_PMTs (std::string FileName) {

    CdPmts_Path = FileName;

    map <int,tuple<int,double,double,double,double,double>> Data_map;

    Read_file(FileName,Data_map,false);

    Map_PMTs = Data_map;
    
};

JUNO_PMTs::JUNO_PMTs (string CdFileName, string WpFileName, string BottomFilename) {

    CdPmts_Path = CdFileName;
    WpPmts_Path = WpFileName;
    BottomPmts_Path = BottomFilename;

    map <int,tuple<int,double,double,double,double,double>> Data_map;
    
    Read_file(CdFileName,Data_map,false);
    Read_file(WpFileName,Data_map,true);
    Read_file(BottomFilename,Data_map,true);

    Map_PMTs = Data_map;
    
}

JUNO_PMTs::~JUNO_PMTs () {
    Map_PMTs.clear();
}

bool JUNO_PMTs::isHama(int PmtNo) {

    if (get<0>(Map_PMTs[PmtNo]) == 1) return true;
    else return false;

}

bool JUNO_PMTs::isNNVT(int PmtNo) {

    if (get<0>(Map_PMTs[PmtNo]) != 1) return true;
    else return false;
    
}

void JUNO_PMTs::SetCdPmts(string CdFileName) {

    if (!CdPmts_Path.empty()) {
        cout << "WARNING: you are resetting an already present Cd PMTs file " << endl;
    } 

    Read_file(CdFileName,Map_PMTs,false);

    return;
}

void JUNO_PMTs::SetWpPmts(string WpFileName) {

    if (!WpPmts_Path.empty()) {
        cout << "WARNING: you are resetting an already present Wp PMTs file " << endl;
    } 

    Read_file(WpFileName,Map_PMTs,true);

    return;
}

void JUNO_PMTs::SetBottomPmts(string BottomFileName) {

    if (!BottomPmts_Path.empty()) {
        cout << "WARNING: you are resetting an already present Bottom PMTs file " << endl;
    } 

    Read_file(BottomFileName,Map_PMTs,true);

    return;
}
