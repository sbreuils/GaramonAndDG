// Copyright (c) 2019 by University Paris-Est Marne-la-Vallee
// Geogebra_c3ga.hpp
// Authors: Vincent Nozick and Stephane Breuils
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Geogebra_c3ga.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief define an entry for a convertor from Garamon to Geogebra files


#ifndef GA_GEOGEBRA_C3GA_HPP__
#define GA_GEOGEBRA_C3GA_HPP__


#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <limits>

#include <c2ga/Mvec.hpp>

#include "c2gaTools.hpp"
#include "Entry.hpp"
#include "Directory.hpp"

const double epsi = 1e-7;


class Viewer_c2ga
{
public:

    // list of entries
    std::list<Entry> entries;


public:

    Viewer_c2ga();

    ~Viewer_c2ga();

    void display() const;

    void removeNameDoublons();

    void render(const std::string &filename);

    template<typename T>
    inline void pushPolygon(std::list<c2ga::Mvec<T> > &mvList, const unsigned int &red = -1, const unsigned int &green = -1, const unsigned int &blue = -1 ){

        // init equation
        std::string equation = "Polygon( ";

        // add all points for polygon
        for(auto& e : mvList){

            // homogeneous coordinates
            if(fabs(e[c2ga::E0]) > std::numeric_limits<T>::epsilon())
                e = e / e[c2ga::E0];

            // add point to equation
            equation += "Point({" + std::to_string(e[c2ga::E1]) + "," + std::to_string(e[c2ga::E2]) +  "})";

            // last element ?
            if(&e != &mvList.back())
                equation += ",";
        }

        // close equation
        equation += ")";

        // push the entry
        entries.push_back(Entry(equation, "poly", red, green, blue));


        // add all points (just points)
        for(auto& e : mvList){
            // add point to equation
            push(e, red*0.5, green*0.5, blue*0.5);
        }
    }


    template<typename T>
    inline int push(c2ga::Mvec<T> mv, const unsigned int &red = -1, const unsigned int &green = -1, const unsigned int &blue = -1 ){
        return push(mv,"",red,green,blue);
    }

    template<typename T>
    inline int push(c2ga::Mvec<T> mv, std::string objectName = "", const unsigned int &red = -1, const unsigned int &green = -1, const unsigned int &blue = -1 ){

        // exrtact multivector type (point / sphere / ....)
        mv.roundZero(epsi);
        std::string mvType = c2ga::whoAmI(mv);
        std::cout << "type = " <<  mvType << std::endl;

        // remove space in the name
        objectName.erase(std::remove(objectName.begin(), objectName.end(), ' '), objectName.end());

        // final equation
        std::string equation;

        // grade 1 ///////////////////////////////////

        // Euclidean vector
        if(mvType == "Euclidean vector") {
            equation = "Vector((" + std::to_string(mv[c2ga::E1]) + "," + std::to_string(mv[c2ga::E2]) +  "))";

            // put a default name
            if(objectName == "")
                objectName = "vec";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // points
        if(mvType == "point (dual tangent trivector)"){

            // homogeneous coordinate to 1
            mv /= mv[c2ga::E0];
            equation = " Point({" + std::to_string(mv[c2ga::E1]) + "," + std::to_string(mv[c2ga::E2]) +  "})";

            // put a default name
            if(objectName == "")
                objectName = "pt";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;        }

        // dual sphere
        if(mvType == "dual sphere" || mvType == "imaginary dual sphere"){
            // back to correct scale
            mv /= mv[c2ga::E0];

            // extract signedradius
            double squaredRadius = mv | mv;
            double radius = sqrt(squaredRadius);

            // extract center
            c2ga::Mvec<double> center = mv;
            center /= center[c2ga::E0];

            // equation
            equation = " Sphere((" + std::to_string(mv[c2ga::E1]) + "," + std::to_string(mv[c2ga::E2]) + std::to_string(fabs(radius)) + ")";

            // put a default name
            if(objectName == "")
                objectName = "sphere";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // dual plane
        if(mvType == "dual plane"){
            equation = std::to_string(mv[c2ga::E1]) + " x ";
            if(mv[c2ga::E2] >= 0) equation += " + ";
            equation += std::to_string(mv[c2ga::E2]) + " y ";
            if(-mv[c2ga::Ei] >= 0) equation += " + ";
            equation += std::to_string(-mv[c2ga::Ei]) + " = 0";

            // put a default name
            if(objectName == "")
                objectName = "plane";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }



        // grade 2 ///////////////////////////////////

        // tangent vector (dual tangent bivector)
        if(mvType == "tangent vector (dual tangent bivector)") {

            // position and orientation
            c2ga::Mvec<T> pos,dir;
            c2ga::extractTangentVector(mv, pos, dir);
            equation = " Vector( Point({" + std::to_string(pos[c2ga::E1]) + "," + std::to_string(pos[c2ga::E2]) + "}), Point({"
                       + std::to_string(pos[c2ga::E1] + dir[c2ga::E1]) + "," + std::to_string(pos[c2ga::E2] + dir[c2ga::E2]) +  "}) )";

            // put a default name
            if(objectName == "")
                objectName = "tangentVec";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // add tangent bivector
            c2ga::Mvec<T> pt1 = pos - 1.0e-5 * dir;
            c2ga::Mvec<T> pt2 = pos + 1.0e-5 * dir;

            equation = "Cylinder(Point({" + std::to_string(pt1[c2ga::E1]) + "," + std::to_string(pt1[c2ga::E2]) + "}), Point({"
                       + std::to_string(pt2[c2ga::E1]) + "," + std::to_string(pt2[c2ga::E2]) + "}),0.3)";
            objectName = "tanBiv";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // pair point (imaginary dual circle)
        if(mvType == "pair point (imaginary dual circle)" || mvType == "imaginary pair point (dual circle)") {

            // extract the 2 points
            c2ga::Mvec<T> pt1,pt2;
            c2ga::extractPairPoint(mv, pt1, pt2);

            // equation of point 1
            equation = " Point({" + std::to_string(pt1[c2ga::E1]) + "," + std::to_string(pt1[c2ga::E2]) +  "})";

            // put a default name
            if(objectName == "")
                objectName = "pp";

            // push the first point
            entries.push_back(Entry(equation, objectName + std::to_string(1), red, green, blue));

            // equation of point 2
            equation = "Point({" + std::to_string(pt2[c2ga::E1]) + "," + std::to_string(pt2[c2ga::E2]) + "})";

            // push the entry
            entries.push_back(Entry(equation, objectName + std::to_string(2), red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // flat point
        if(mvType == "flat point") {

            // extract the point
            c2ga::Mvec<T> pt;
            c2ga::extractFlatPoint(mv, pt);

            // equation of point
            equation = " Point({" + std::to_string(pt[c2ga::E1]) + "," + std::to_string(pt[c2ga::E2]) + "})";

            // put a default name
            if(objectName == "")
                objectName = "flatPt";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }







        // grade 3 ///////////////////////////////////

        // tangent vector (dual tangent bivector)
        if(mvType == "tangent bivector (dual tangent vector)") {

            // position and orientation
            c2ga::Mvec<T> pos,dir;
            c2ga::extractTangentBivector(mv, pos, dir);
            equation = " Vector( Point({" + std::to_string(pos[c2ga::E1]) + "," + std::to_string(pos[c2ga::E2]) + "}), Point({"
                       + std::to_string(pos[c2ga::E1] + dir[c2ga::E1]) + "," + std::to_string(pos[c2ga::E2] + dir[c2ga::E2]) + "}) )";

            // put a default name
            if(objectName == "")
                objectName = "tangentVec";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // add tangent bivector
            c2ga::Mvec<T> pt1 = pos - 1.0e-5 * dir;
            c2ga::Mvec<T> pt2 = pos + 1.0e-5 * dir;

            equation = "Cylinder(Point({" + std::to_string(pt1[c2ga::E1]) + "," + std::to_string(pt1[c2ga::E2]) + "}), Point({"
                       + std::to_string(pt2[c2ga::E1]) + "," + std::to_string(pt2[c2ga::E2]) +  "}),0.3)";
            objectName = "tanBiv";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // circle
        if(mvType == "circle (imaginary dual pair point)" || mvType == "imaginary circle (dual pair point)") {

            T radius;
            c2ga::Mvec<T> center, direction;
//            extractDualCircle(mv.dual(), radius,  center, direction);
            equation = "Circle( Point({" + std::to_string(center[c2ga::E1]) + "," + std::to_string(center[c2ga::E2]) +"}), "
                       + std::to_string(fabs(radius))
                       + ", Vector((" + std::to_string(direction[c2ga::E1]) + "," + std::to_string(direction[c2ga::E2]) + ")))";

            // put a default name
            if(objectName == "")
                objectName = "circle";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }




        // dual flat point
        if(mvType == "dual flat point") {

            // extract the point
            c2ga::Mvec<T> pt;
            c2ga::extractFlatPoint(mv.dual(), pt);

            // equation of point
            equation = " Point({" + std::to_string(pt[c2ga::E1]) + "," + std::to_string(pt[c2ga::E2]) +  "})";

            // put a default name
            if(objectName == "")
                objectName = "dualFlatPt";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;        }




        return EXIT_FAILURE;
    }

};




#endif // GA_GEOGEBRA_C3GA_HPP__
