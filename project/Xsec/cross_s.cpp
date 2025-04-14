#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>
#include "cross_s.h"

cross_sect::cross_sect(const string& path)
{
    MAP int_map;
    int_map["EFFECTIVE"]=EFFECTIVE;
    int_map["EXCITATION"]=EXCITATION;
    int_map["IONIZATION"]=IONIZATION;
    int_map["ATTACHMENT"]=ATTACHMENT;

    //string filename="N2.txt";     //nome del file
    ifstream file;                //stream object
    string line;                  //riga da acquisire
    string energia,formula;
    bool fl=0;
    unsigned n=0;
    unsigned counter=0;

    file.open(path,ifstream::in);  //apri il file
    if(file.is_open())
    {
        while(getline(file,line))
        {
            if (line.compare("EFFECTIVE")==0 | line.compare("EXCITATION")==0 | line.compare("IONIZATION")==0 | line.compare("ATTACHMENT")==0)
            {
                n++;
            }
        }
    }
    file.close();

    //riapri il file
    file.open(path,ifstream::in);

    //VARIABILE di estrazione
    len=n;
    tab.reserve(n);

    if(file.is_open())
    {
        int i=0;
        while(getline(file,line))
        {
            if (line.compare("EFFECTIVE")==0 | line.compare("EXCITATION")==0 | line.compare("IONIZATION")==0 | line.compare("ATTACHMENT")==0)
            {
                fl=1;
                MAP_IT it=int_map.find(line);
                cout<<"ok"<<endl;
                tab[i].interact=it->second;
                cout<<"ok"<<endl;
                getline(file,formula);
                tab[i].react=formula;
                cout<<"ok"<<endl;
                getline(file,energia);
                stringstream ss_en(energia);
                cout<<"ok"<<endl;
                ss_en>>tab[i].en_avg;
                cout<<endl<<"tipo di interazione "<<"n."<<i+1<<": "<<line<<endl;
                cout<<"formula "<<"n."<<i+1<<": "<<tab[i].react<<endl;
                cout<<"energia "<<"n."<<i+1<<": "<<tab[i].en_avg<<endl<<endl;
            }

            if(fl==1)
            {
                if(line[0]=='-')
                {
                    counter++;
                }
                stringstream ss(line);
                double x,y;

                if(counter<2)
                {
                    if(ss>>x>>y)
                    {
                        ss>>x>>y;
                        cout<<x<<" "<<y<<endl;
                        tab[i].energy.push_back(x);
                        tab[i].sect.push_back(y);
                        cout<<"ok"<<endl;
                    }

                }
                else
                {
                    fl=0;
                    counter=0;
                    i++;
                }
            }
        }
    }
}

//cross_sect::~cross_sect()
//{
//    delete[](tab);
//}


//string cross_sect::formula(const int & ind)
//{
//    string form;
//    form=tab[ind].react;
//    return form;
//}


//INTER cross_sect::interaction(const int & ind)
//{
//    INTER inter;
//    inter=tab[ind].interact;
//    return inter;
//}


//VECT cross_sect::energy(const int & ind)
//{
//    VECT en_v;
//    en_v=tab[ind].energy;
//    return en_v;
//}

//double cross_sect::energy(const int & ind,const int & pos)
//{
//    double val;
//    val=tab[ind].energy[pos];
//    return val;
//}

//VECT cross_sect::section(const int & ind)
//{
//    VECT en_v;

//    en_v=tab[ind].sect;

//    return en_v;
//}

//double cross_sect::section(const int & ind,const int & pos)
//{
//    double val;
//    val=tab[ind].sect[pos];
//    return val;
//}

//void cross_sect::remove(const int & ind,const int & off)
//{
//    tab[ind].energy.erase(tab[ind].energy.begin()+off);
//    tab[ind].sect.erase(tab[ind].sect.begin()+off);
//}

void cross_sect::print()
{
    map<INTER,string> int_map;
    int_map[EFFECTIVE]="EFFECTIVE";
    int_map[EXCITATION]="EXCITATION";
    int_map[IONIZATION]="IONIZATION";
    int_map[ATTACHMENT]="ATTACHMENT";

    map<INTER,string>::iterator it;

    for(int ind=0;ind<len;ind++)
    {
        cout<<"INTERAZIONE n."<<ind+1<<"\n";
        it=int_map.find(tab[ind].interact);
        cout<<"Tipo interazione: "<<it->second<<"\n";
        cout<<"Formula: "<<tab[ind].react<<"\n";
        cout<<"Energia media: "<<tab[ind].en_avg<<" eV"<<"\n\n";
        cout<<"E[eV]"<<"\t"<<"Sez[m2]"<<"\n";
        for(int k=0;k<tab[ind].energy.size();k++)
        {
            cout<<tab[ind].energy[k]<<"\t"<<tab[ind].sect[k]<<"\n";
        }
        cout<<"-----------"<<"\n\n";

    }
}
