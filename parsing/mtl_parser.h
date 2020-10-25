#pragma once

#include <map>
#include <string>
#include <fstream>
#include "../entities/material.h"
#include "commons.h"

class MtlParser {
public:
    static std::map<std::string, Material> Parse(const std::string& filename) {
        std::map<std::string, Material> result;
        std::string current_name;
        std::ifstream file(filename);
        if (file.is_open()) {
            std::string line;
            while (getline(file, line)) {
                std::string type = GetFirstStringBeforeSpace(line);
                auto split_line = SplitBySpacesAndTabs(line);
                if (type == "illum") {
                    double illum = atof(split_line.back().c_str());
                    result[current_name].illum = illum;
                } else if (type == "newmtl") {
                    current_name = line.substr(7);
                    result[current_name] = {};
                } else if (type == "Ka") {
                    Vector rgb(atof(split_line[1].c_str()), atof(split_line[2].c_str()),
                               atof(split_line[3].c_str()));
                    result[current_name].Ka = rgb;
                } else if (type == "Kd") {
                    Vector rgb(atof(split_line[1].c_str()), atof(split_line[2].c_str()),
                               atof(split_line[3].c_str()));
                    result[current_name].Kd = rgb;
                } else if (type == "Ks") {
                    Vector rgb(atof(split_line[1].c_str()), atof(split_line[2].c_str()),
                               atof(split_line[3].c_str()));
                    result[current_name].Ks = rgb;
                } else if (type == "Ke") {
                    Vector rgb(atof(split_line[1].c_str()), atof(split_line[2].c_str()),
                               atof(split_line[3].c_str()));
                    result[current_name].Ke = rgb;
                } else if (type == "Ns") {
                    double value = atof(split_line[1].c_str());
                    result[current_name].Ns = value;
                } else if (type == "Ni") {
                    double value = atof(split_line[1].c_str());
                    result[current_name].Ni = value;
                } else if (type == "d") {
                    double value = atof(split_line[1].c_str());
                    result[current_name].d = value;
                    result[current_name].Tr = 1 - value;
                } else if (type == "Tr") {
                    double value = atof(split_line[1].c_str());
                    result[current_name].Tr = value;
                    result[current_name].d = 1 - value;
                }
            }
        }
        return result;
    }
};
