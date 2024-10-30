#ifndef LISTOFMASS_H
#define LISTOFMASS_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include "MolarMass.h"


class ListOfMass
{
  public:
    typedef std::map<std::string,double> element;
    typedef std::pair<std::string,std::string> PAIR;
    typedef std::vector<PAIR> VECT;

  public:
    std::vector<MolMass> CompList;

  public:
      ListOfMass();
      ListOfMass(const ListOfMass & LL);
      ListOfMass operator=(const ListOfMass & LL);
  public:
      //void print();
      void import(const std::string & list);

  public:
      MolMass operator[](const int ind);
      std::vector<double> get_MM();
      std::vector<std::string> get_namelist();

      size_t length();

};

#endif
