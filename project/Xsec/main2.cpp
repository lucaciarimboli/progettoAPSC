#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>
#include "cross_s.h"

using namespace std;


int main(int argc, char *argv[])
{
    string path("C:/Disco D/Fasani/interazioni/H2O.txt");

    cross_sect h2o(path);

    int k=h2o.len-1;
    //    VECT en=h2o.energy(k);
    //    VECT sec=h2o.section(k);
//    for(int i=0;i<h2o.tab[k].energy.size();i++)
//    {
//        cout<<h2o.energy(k,i)<<"\t"<<h2o.section(k,i)<<endl;
//    }

//    int r=0;
//    h2o.remove(k,r);
//    for(int i=0;i<h2o.tab[k].energy.size();i++)
//    {
//        cout<<h2o.energy(k,i)<<"\t"<<h2o.section(k,i)<<endl;
//    }

    h2o.print();

}
