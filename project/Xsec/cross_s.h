#ifndef CROSS_S_H
#define CROSS_S_H

#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

enum INTER
{
    EFFECTIVE,
    EXCITATION,
    IONIZATION,
    ATTACHMENT
};

typedef vector<double> VECT;
typedef map<string,INTER> MAP;
typedef MAP::iterator MAP_IT;

class tabella
{
public:
    VECT energy,sect;
    double en_avg;
    string react;
    INTER interact;
};


class cross_sect
{
public:
    vector<tabella> tab;
    unsigned len;

public:
    cross_sect(const string & path);
//    ~cross_sect();
    //accesso ai membri
//    VECT energy(const int & ind);
//    double energy(const int & ind,const int & pos);
//    VECT section(const int & ind);
//    double section(const int & ind,const int & pos);
//    string formula(const int & ind);
//    INTER interaction(const int & ind);
//    void remove(const int & ind, const int & off);
//    //stampa
    void print();
};

#endif // CROSS_S_H
