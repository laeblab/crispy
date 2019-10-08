#include <sstream>
#include <sys/stat.h>
#include <sys/param.h>
#include <unistd.h>
#include <algorithm>

#include "crispy++.h"


bool does_file_exist(const std::string& name) {
	struct stat buffer;

	return stat(name.c_str(), &buffer) == 0;
}

bool is_string_number(const std::string& s) {
    std::string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::string get_working_path()
{
   char temp[MAXPATHLEN];
   return ( getcwd(temp, MAXPATHLEN) ? std::string( temp ) : std::string("") );
}

unsigned dna_to_bin(const std::string& dna) {
    unsigned bin_int = 0;

    //iterating forward to make PAM most nuc least significant bit
    for (const auto nuc : dna) {
        if (nuc == 'A') bin_int = (bin_int << 2);
        else if (nuc == 'T') bin_int = (bin_int << 2) | 1;
        else if (nuc == 'C') bin_int = (bin_int << 2) | 2;
        else if (nuc == 'G') bin_int = (bin_int << 2) | 3;
    }

    return bin_int;
}

unsigned dna_to_revcompl_bin(const std::string& dna) {
    unsigned bin_int = 0;

    // NOW SINCE THIS IS THE OPPOSITE STRAND AND WE WANT 5'-3' gRNA's WE NEED TO TRANSLATE IN OPPOSITE
    // DIRECTION AND WITH COMPLEMENTARY BASES
    //		fwd	rev
    // A = 	00	01
    // T = 	01	00
    // C =	10	11
    // G =	11	10
    for (size_t i = 0; i < dna.size(); ++i) {
        const auto nuc = dna[i];
        const auto shift = 2 * i;

        if (nuc == 'A') bin_int |= 1 << shift;
        else if (nuc == 'T') bin_int |= 0 << shift;
        else if (nuc == 'C') bin_int |= 3 << shift;
        else if (nuc == 'G') bin_int |= 2 << shift;
    }

    return bin_int;
}

void make_uppercase(std::string& str) {
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}
