// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Entry.hpp
/// \author Stephane Breuils
/// \brief define an entry for a convertor from Garamon to Geogebra files


#ifndef GA_ENTRY_HPP__
#define GA_ENTRY_HPP__

#include <iostream>
#include <vector>
#include <string>


class Entry
{
public:

    // object equation or geogebra formulation
    std::string _equation;

    // object name
    std::string _objectName;

    // color (from 0 to 1.)
    std::vector<double> _color;


public:

    Entry();

    //Entry(const std::string &equation, const std::string &objectName, const std::vector<int> &color);

    Entry(const std::string &equation, const std::string &objectName, const double &red = -1, const double &green = -1, const double &blue = -1);

    ~Entry();

    void display() const;
};


#endif // GA_ENTRY_HPP__
