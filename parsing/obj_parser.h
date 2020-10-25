#pragma once

#include <vector>
#include <string>
#include "mtl_parser.h"
#include "../entities/sphere.h"
#include "../entities/light.h"
#include "../entities/figure.h"
#include "../entities/vertex.h"
#include "../geometry/vector.h"

class ObjParser {
public:
    void Parse(const std::string& filename, std::vector<Figure>* figures,
               std::vector<Sphere>* spheres, std::vector<Light>* lights) {
        std::vector<Vector> vertices;
        std::vector<std::pair<double, double>> vt;
        std::vector<Vector> vn;
        std::string current_mtl;
        std::map<std::string, Material> mtls;
        std::ifstream file(filename);
        if (file.is_open()) {
            std::string line;
            while (getline(file, line)) {
                std::string type = GetFirstStringBeforeSpace(line);
                auto split_line = SplitBySpacesAndTabs(line);
                if (type == "mtllib") {
                    std::string mtl_filename = line.substr(7);
                    mtls = MtlParser::Parse(GetMtlFilename(filename, mtl_filename));
                } else if (type == "usemtl") {
                    current_mtl = line.substr(7);
                } else if (type == "v") {
                    vertices.push_back({atof(split_line[1].c_str()), atof(split_line[2].c_str()),
                                        atof(split_line[3].c_str())});
                } else if (type == "vt") {
                    vt.push_back({atof(split_line[1].c_str()), atof(split_line[2].c_str())});
                } else if (type == "vn") {
                    vn.push_back({atof(split_line[1].c_str()), atof(split_line[2].c_str()),
                                  atof(split_line[3].c_str())});
                } else if (type == "P") {
                    Vector point = {atof(split_line[1].c_str()), atof(split_line[2].c_str()),
                                    atof(split_line[3].c_str())};
                    Vector rgb = {atof(split_line[4].c_str()), atof(split_line[5].c_str()),
                                  atof(split_line[6].c_str())};
                    lights->push_back({point, rgb});
                } else if (type == "S") {
                    Vector center = {atof(split_line[1].c_str()), atof(split_line[2].c_str()),
                                     atof(split_line[3].c_str())};
                    double radius = atof(split_line[4].c_str());
                    spheres->emplace_back(mtls[current_mtl], center, radius, spheres->size());
                } else if (type == "f") {
                    auto parsed_vertices = ParseFigure(split_line, vertices, vt, vn);
                    figures->emplace_back(mtls[current_mtl], parsed_vertices, figures->size());
                }
            }
        }
    }

    // TODO refactoring
    std::vector<Vertex> ParseFigure(const std::vector<std::string>& split_line,
                                    const std::vector<Vector>& vertices,
                                    const std::vector<std::pair<double, double>>& vt,
                                    const std::vector<Vector>& vn) {
        const char delimiter = '/';
        int v_n = vertices.size(), vt_n = vt.size(), vn_n = vn.size();
        int count = split_line.size();
        std::vector<Vertex> result;
        for (int i = 1; i < count; ++i) {
            auto current = split_line[i];
            int delimiter_index = current.find_first_of(delimiter);
            int v_index = atoi(current.substr(0, delimiter_index).c_str());
            if (v_index > 0) {
                --v_index;
            } else {
                v_index += v_n;
            }
            Vertex current_v(vertices[v_index]);
            if (delimiter_index != std::string::npos) {
                int delimiter_index2 = current.find_first_of(delimiter, delimiter_index + 1);
                if (delimiter_index2 != delimiter_index + 1) {
                    int vt_index = atoi(
                        current.substr(delimiter_index + 1, delimiter_index2 - delimiter_index - 1)
                            .c_str());
                    if (vt_index > 0) {
                        --vt_index;
                    } else {
                        vt_index += vt_n;
                    }
                    current_v.vt = vt[vt_index];
                }
                if (delimiter_index2 != std::string::npos) {
                    int vn_index = atoi(current.substr(delimiter_index2 + 1).c_str());
                    if (vn_index > 0) {
                        --vn_index;
                    } else {
                        vn_index += vn_n;
                    }
                    current_v.vn = vn[vn_index];
                    current_v.is_normal_set = true;
                }
            }
            result.push_back(current_v);
        }
        return result;
    }

private:
    std::string GetMtlFilename(const std::string& obj_filename, const std::string& mtl_filename) {
        int last_delimiter = obj_filename.find_last_of("/");
        return obj_filename.substr(0, last_delimiter + 1) + mtl_filename;
    }
};
