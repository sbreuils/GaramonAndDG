// Copyright (c) 2020 by National Institute of Informatics
// Entry.cpp
// Author: Stephane Breuils
// Contact: breuils@nii.ac.jp
// Licence MIT
// A a copy of the MIT License is given along with this program


#include "Entry.hpp"


// Entry::Entry(const std::string &equation, const std::string &objectName, const std::vector<unsigned int> &color)
//    : _equation(equation), _objectName(objectName), _color(color)
//    {
//    }


Entry::Entry(const std::string &equation, const std::string &objectName, const double &red, const double &green, const double &blue)
   : _equation(equation), _objectName(objectName)
   {
      _color = std::vector<double>(3);

      double redComponent = red;
      double greenComponent = green;
      double blueComponent = blue;

      if(red>1.0 && red <256.0)
        redComponent = (red)/255.0;
      if(green>1.0 && green <256.0)
        greenComponent = (green)/255.0;
      if(blue>1.0 && blue<256.0)
        blueComponent = (blue)/255.0;

      _color[0] = redComponent;
      _color[1] = greenComponent;
      _color[2] = blueComponent;
   }



Entry::~Entry()
{}

void Entry::display() const
{
    std::cout << "Entry " << std::endl;
    std::cout << "  object  : " << _objectName << std::endl;
    std::cout << "  equation: " << _equation << std::endl;
    std::cout << "  color   : (";
    for(int i=0; i<(int)_color.size()-1; ++i)
        std::cout << std::to_string(_color[i]) << ", ";
    std::cout << std::to_string(_color[_color.size()-1]) << ")" << std::endl;
    std::cout << std::endl;
}
