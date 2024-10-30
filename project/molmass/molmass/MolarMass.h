#ifndef MOLARMASS_H
#define MOLARMASS_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>

//using namespace std;


class MolMass
{
  public:
    typedef std::map<std::string,double> element;
    typedef std::pair<std::string,std::string> PAIR;
    typedef std::vector<PAIR> VECT;

  public:
    double MMol;
    std::string name;

  public:
      MolMass();
      MolMass(const MolMass & MM);
      //MolMass(const std::string& compound);
      MolMass operator=(const MolMass & MM);
  public:
    //void print();
    void setval(const std::string & compound);

  public:
    double mass();
    std::string subst();

};

#endif

