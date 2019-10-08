#ifndef TOOLS___H_
#define TOOLS___H_

#include <string>
#include <vector>

bool does_file_exist(const std::string& name);

// checks if a string is a proper number which is useful before converting it as such
bool is_string_number(const std::string& s);

// these next two somehow split a string into a vector of strings based on delim
std::vector<std::string> split(const std::string &s, char delim);

std::string get_working_path();

unsigned dna_to_bin(const std::string& dna);

unsigned dna_to_revcompl_bin(const std::string& dna);

void make_uppercase(std::string& str);

#endif