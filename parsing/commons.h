#pragma once

#include <string>
#include <vector>

std::string GetFirstStringBeforeSpace(const std::string& str) {
    int i = 0;
    int n = str.size();
    while (i < n && (str[i] == ' ' || str[i] == '\t')) {
        ++i;
    }
    std::string result;
    while (i < n && !(str[i] == ' ' || str[i] == '\t')) {
        result += str[i++];
    }
    return result;
}

std::vector<std::string> SplitBySpacesAndTabs(const std::string& str) {
    int i = 0;
    int n = str.size();
    std::vector<std::string> result;
    while (i < n) {
        while (i < n && (str[i] == ' ' || str[i] == '\t')) {
            ++i;
        }
        std::string current;
        while (i < n && !(str[i] == ' ' || str[i] == '\t')) {
            current += str[i++];
        }
        if (!current.empty()) {
            result.push_back(current);
        }
    }
    return result;
}