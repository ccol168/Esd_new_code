#include <cstdlib>
#include <map>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

#ifndef JUNO_PMTS_H
#define JUNO_PMTS_H

class JUNO_PMTs {

    private :
        map <int,tuple<int,double,double,double,double,double>> Map_PMTs;
        string CdPmts_Path;
        string WpPmts_Path;
        string BottomPmts_Path;

    public : 

        JUNO_PMTs();
        JUNO_PMTs(string);
        JUNO_PMTs(string, string, string);

        ~JUNO_PMTs();

        void SetCdPmts(string);
        void SetWpPmts(string);
        void SetBottomPmts(string);

        string GetCdFileName() {return CdPmts_Path;}
        string GetWpFileName() {return WpPmts_Path;}
        string GetBottomFileName() {return BottomPmts_Path;}

        int GetType(int PmtNo) {return get<0>(Map_PMTs[PmtNo]);}
        double GetX(int PmtNo) {return get<1>(Map_PMTs[PmtNo]);}
        double GetY(int PmtNo) {return get<2>(Map_PMTs[PmtNo]);}
        double GetZ(int PmtNo) {return get<3>(Map_PMTs[PmtNo]);}
        double GetTheta(int PmtNo) {return get<4>(Map_PMTs[PmtNo]);}
        double GetPhi(int PmtNo) {return get<5>(Map_PMTs[PmtNo]);}

        tuple <double, double, double> GetCartesian (int PmtNo) {return make_tuple(get<1>(Map_PMTs[PmtNo]),get<2>(Map_PMTs[PmtNo]),get<3>(Map_PMTs[PmtNo]));}
        tuple <double, double> GetSpherical (int PmtNo) {return make_tuple(get<4>(Map_PMTs[PmtNo]),get<5>(Map_PMTs[PmtNo]));}
        tuple <int,double,double,double,double,double> GetAll (int PmtNo) {return Map_PMTs[PmtNo];}

        bool isHama(int PmtNo);
        bool isNNVT(int PmtNo);

};

#endif
