
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "MolarMass.h"
#include "ListOfMass.h"

using namespace std;
typedef map<string,double> element;
typedef element::iterator IT;
typedef pair<string,string> PAIR;
typedef vector<PAIR> VECT;


int main()
{
    //string compound="CH3(CH2(CH2)2)^2+CH3";
    string compound="[CH3(CH2(CH2)2)^2+CH3 C12H24O12 NH4^+ CaHCO3]";

    MolMass prova;
    ListOfMass provaL;
    provaL.import(compound);
    ListOfMass copia=provaL;
    std::vector<double> v;
    std::vector<string> s;

    for(int i=0;i<copia.length();i++)
    {
        cout<<"massa molare di "<<copia[i].name<<" = "<<copia[i].MMol<<endl;
    }

    v=copia.get_MM();
    s=copia.get_namelist();


    // prova.setval(compound);
    // cout<<"massa molare di "<<prova.name<<" = "<<prova.MMol<<endl;

    // for(int i=0;i<mmlist.size();i++)
    // {
    //     cout<<"MM="<<mmlist[i]<<endl;
    // }
    cout<<"\n";

    for(int i=0;i<v.size();i++)
    {
        cout<<"massa molare di "<<s[i]<<" = "<<v[i]<<endl;
    }

}