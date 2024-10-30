#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include "ListOfMass.h"

using namespace std;
typedef std::map<std::string,double> element;
typedef std::pair<std::string,std::string> PAIR;
typedef std::vector<PAIR> VECT;


ListOfMass::ListOfMass(){}

ListOfMass::ListOfMass(const ListOfMass & LL)
{
    CompList=LL.CompList;
}

ListOfMass ListOfMass::operator=(const ListOfMass & LL)
{
    CompList=LL.CompList;
    return *this;
}


void ListOfMass::import(const string & list)
{
    MolMass el;
    bool fl=0;

    //check syntax
    string permitted_char="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[];, +-.^";
    for(int i=0;i<list.size();i++)
    {
        if (permitted_char.find(list[i])==string::npos)
        {
            cout<<"UNEXPECTED CHARACTER"<<endl;
            el.MMol=NAN;
            CompList.push_back(el);
            fl=1;
            break;
        }
    }

    if (!fl)
    {
        if(list[0]=='[')
        {
            int count=0;
            string list_clean;

            for(int i=0;i<list.size();i++)
            {
                if((list[i]!='[') && (list[i]!=']'))
                {
                    list_clean.push_back(list[i]);
                }
            }

            stringstream ss(list_clean);
            while (ss.good())
            {
                string substr;
                getline(ss, substr, ' ');
                el.setval(substr);
                CompList.push_back(el);
                //el.MMol=0;
            }
    //        for (size_t i = 0; i < names.size(); i++)
    //            cout << names[i] << "----"; cout<<endl;
        }
        else
        {
            el.setval(list);
            CompList.push_back(el);
        }
    }
}


MolMass ListOfMass::operator[](const int ind)
{
    return CompList[ind];
}

size_t ListOfMass::length()
{
    return CompList.size();
}

vector<double> ListOfMass::get_MM()
{
    vector<double> v;
    for(int i=0;i<this->length();i++)
    {
        v.push_back(CompList[i].MMol);
    }
    return v;
}

vector<string> ListOfMass::get_namelist()
{
    vector<string> s;
    for(int i=0;i<this->length();i++)
    {
        s.push_back(CompList[i].name);
    }
    return s;
}