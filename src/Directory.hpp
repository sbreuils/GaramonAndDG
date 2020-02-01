// Author: Stephane Breuils and Vincent Nozick
// Contact: breuils@nii.ac.jp
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Directory.hpp
/// \author Stephane Breuils
/// \brief Basic tools to manipulate files and directories.


#ifndef GA_DIRECTORY_HPP__
#define GA_DIRECTORY_HPP__

#include <string>
#include <vector>

void makeDirectory(const std::string &dirName);

bool directoryExists(const std::string &dirName);

bool directoryOrFileExists(const std::string &dirName);

bool directoryOrFileExists_ifstream(const std::string& name);

std::string readFile(const std::string& fileName);

bool writeFile(const std::string& data, const std::string& fileName);

void substitute(std::string &data, const std::string &pattern, const std::string &replaceBy);

bool copyBin(const std::string &src, const std::string &dest);

bool copyText(const std::string &src, const std::string &dest);


#endif  // GA_DIRECTORY_HPP__