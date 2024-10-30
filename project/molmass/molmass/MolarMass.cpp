#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include "MolarMass.h"

using namespace std;
typedef std::map<std::string,double> element;
typedef std::pair<std::string,std::string> PAIR;
typedef std::vector<PAIR> VECT;



MolMass::MolMass()
{
    MMol=0;
    name="";
}

MolMass::MolMass(const MolMass & MM)
{
    MMol=MM.MMol;
    name=MM.name;
}

MolMass MolMass::operator=(const MolMass & MM)
{
    MMol=MM.MMol;
    name=MM.name;
    return *this;
}

void MolMass::setval(const string & compound)
{
        // molar masses of the elements:
    double H = 1.00794, He = 4.002602, Li = 6.941, Be = 9.012182, B = 10.811, C = 12.011, N = 14.00674,
        O = 15.9994, F = 18.9984032, Ne = 20.1797, Na = 22.989768, Mg = 24.305, Al = 26.981539, Si = 28.0855,
        P = 30.973762, S = 32.066, Cl = 35.4527, Ar = 39.948, K = 39.0983, Ca = 40.078, Sc = 44.95591,
        Ti = 47.88, V = 50.9415, Cr = 51.9961, Mn = 54.93805, Fe = 55.847, Co = 58.9332, Ni = 58.69,
        Cu = 63.546, Zn = 65.39, Ga = 69.723, Ge = 72.61, As = 74.92159, Se = 78.96, Br = 79.904, Kr = 83.8,
        Rb = 85.4678, Sr = 87.62, Y = 88.90585, Zr = 91.224, Nb = 92.90638, Mo = 95.94, Tc = 98.9063,
        Ru = 101.07, Rh = 102.9055, Pd = 106.42, Ag = 107.8682, Cd = 112.411, In = 114.82, Sn = 118.71,
        Sb = 121.75, Te = 127.6, I = 126.90447, Xe = 131.29, Cs = 132.90543, Ba = 137.327, La = 138.9055,
        Ce = 140.115, Pr = 140.90765, Nd = 144.24, Pm = 146.9151, Sm = 150.36, Eu = 151.965, Gd = 157.25,
        Tb = 158.92534, Dy = 162.5, Ho = 164.93032, Er = 167.26, Tm = 168.93421, Yb = 173.04, Lu = 174.967,
        Hf = 178.49, Ta = 180.9479, W = 183.85, Re = 186.207, Os = 190.2, Ir = 192.22, Pt = 195.08,
        Au = 196.96654, Hg = 200.59, Tl = 204.3833, Pb = 207.2, Bi = 208.98037, Po = 208.9824, At = 209.9871,
        Rn = 222.0176, Ac = 223.0197, Th = 226.0254, Pa = 227.0278, U = 232.0381, Np = 231.0359, Pu = 238.0289,
        Am = 237.0482, Cm = 244.0642, Bk = 243.0614, Cf = 247.0703, Es = 247.0703, Fm = 251.0796, Md = 252.0829,
        No = 257.0951, Lr = 258.0986, Rf = 259.1009, Db = 260.1053, Sg = 261.1087, Bh = 262.1138, Hs = 263.1182,
        Mt = 262.1229,
        Nn	= 1;	//Nn = Not named - for not named substances

    element loe; //loe=list of elements;
    loe["H"]=H;
    loe["He"]=He;
    loe["Li"]=Li;
    loe["Be"]=Be;
    loe["B"]=B;
    loe["C"]=C;
    loe["N"]=N;
    loe["O"]=O;
    loe["F"]=F;
    loe["Ne"]=Ne;
    loe["Na"]=Na;
    loe["Mg"]=Mg;
    loe["Al"]=Al;
    loe["Si"]=Si;
    loe["P"]=P;
    loe["S"]=S;
    loe["Cl"]=Cl;

    loe["Ca"]=Ca;

    name=compound;
    int elem_ind=-1,br_level=0;
    vector<int> fact;
    //vector<double> mass;
    vector<string> chartype(compound.size(),"");
    PAIR empt("","");
    VECT name_vect;

    string nametemp;

    vector<double> MM;    //molar masses of each element
    MolMass MMsub;
    string br_sub;

    for(int ind=0;ind<compound.size();ind++)
    {
        //cout<<"ind="<<ind<<endl;
        if ((compound[ind]>= 'A') && (compound[ind]<= 'Z'))
        {
            elem_ind++;
            fact.push_back(1);

            chartype[ind]="L";        //L stands for "upper case letter"

            string character(1,compound[ind]); //convert char to string

            empt.first=character;
            name_vect.push_back(empt);

        }

        else if ((compound[ind]>= 'a') && (compound[ind]<= 'z'))
        {
            string character(1,compound[ind]);
            name_vect[elem_ind].second=character;
            chartype[ind]="l";          //l stands for "lower case letter"
        }

        else if ((compound[ind]>= '0') && (compound[ind]<= '9'))
        {
            if(compound[ind-1]!= '^')
            {
                string character(1,compound[ind]);
                stringstream ss(character);
                int nu;
                ss>>nu;
                chartype[ind]="num";
                if (chartype[ind-1].compare("num")==0)
                {
                    fact[elem_ind]=fact[elem_ind]*10+nu;
                }
                else
                {
                    fact[elem_ind]=nu;
                }
            }
        }

        else if(compound[ind]== '^')
        {
            chartype[ind]="apex";
            while((compound[ind]!='+') && (compound[ind]!='-'))
            {
                ind++;
                chartype[ind]="charge";
            }
        }

        else if (compound[ind]=='(')
        {
            chartype[ind]="br_op";
            elem_ind++;
            fact.push_back(1);
            empt.first="(";
            empt.second="";
            name_vect.push_back(empt);

            br_level=1;
            int br_ind=-1;           //index of the substring inside the bracket
            while(br_level>0)
            {
                br_ind++;
                ind++;
                if(compound[ind] == '(')
                {
                    br_level++;
                }
                else if	(compound[ind] == ')')
                {
                    br_level--;
                }

                if (br_level > 0)
                {
                    br_sub.push_back(compound[ind]);
                    chartype[ind]="sub";
                }
                else
                {
                    chartype[ind]="br_cl";
                }
            }
            //cout<<"substring: "<<br_sub<<endl;
            MMsub.setval(br_sub);
        }
    }
    // cout<<"factor size: "<<fact.size()<<endl;
    // cout<<"massa sub: "<<MMsub<<endl;

    MM.reserve(fact.size());

    for(int k=0;k<fact.size();k++)
    {
        double mass_temp=0;
        if(name_vect[k].first.compare("(")!=0)
        {
            nametemp=name_vect[k].first+name_vect[k].second;
            //cout<<nametemp<<"-"<<fact[k]<<endl;
            mass_temp=loe[nametemp];
        }
        else
        {
            mass_temp=MMsub.MMol;
        }
        MM.push_back(mass_temp);
        //cout<<MM[k]<<endl;
    }

    // molar mass
    double MMol_temp=0;
    double prod=0;
    for(int k=0;k<fact.size();k++)
    {
        prod=MM[k]*fact[k];
        MMol_temp=MMol_temp+prod;
    }
    MMol=MMol_temp;

}

double MolMass::mass()
{
    return MMol;
}

std::string MolMass::subst()
{
    return name;
}

