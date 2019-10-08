/*============================================================================
 * Name        : crispy++
 * Author      : Lasse Ebdrup Pedersen
 * Version     : 2.0.0
 * Copyright   :
 * Description : crispy++ has 2 functions
 *
 * 1:	count:	(indexing) find and count all PAMs in genome.
 * 2:	score:	Use index to score list of targets
 *
 *============================================================================
*/
#include <iostream>

#include "crispy++.h"
int debug = 0; // set to 1 to see count debug

void crispy_cmd_line(std::vector<std::string> args) {
	std::string prg_name = args[0].substr(args[0].find_last_of('/')+1);
	std::cerr << "Usage: " << prg_name << " <count|score|version>" << std::endl
			<< "Additional help is available by typing e.g." << std::endl
			<< "\t" << prg_name << " count" << std::endl
			<< "\tor" << std::endl
			<< "\t" << prg_name << " score" << std::endl;
}

int main(int argc, char **argv) {
	//start the program
	std::cout << "<<< !CRISPy TIME! >>>" << std::endl;
	auto start=std::chrono::system_clock::now();	//start a program timer

	//convert command line to vector string
	std::vector<std::string> args(argv, argv + argc);
	if (argc < 2) {
		crispy_cmd_line(args);
	} else {
		if (args[1] == "count") {
			count(args);
		} else if ((args[1] == "score")) {
			score(args);
		} else if ((args[1] == "version")) {
			std::cout << std::endl << "Version 2.0.0, by Lasse Ebdrup Pedersen, lasse.ebdrup@gmail.com" << std::endl;
		} else {
			crispy_cmd_line(args);
		}
	}
	std::cout << std::endl << "Total runtime: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() << " secs" << std::endl;
}
