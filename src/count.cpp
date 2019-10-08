#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include "crispy++.h"
#include "tools.h"


void count_forward_gRNA(std::unordered_map<unsigned, unsigned>& target_store, const std::string& gRNA) {
	// check if there N's or god forbid other characters in our target sequence
	if (gRNA.find_first_not_of("ATCG") == std::string::npos) {
		//is nice DNA =)
		//then we need to do a binary encoding
		const unsigned bin_int = dna_to_bin(gRNA);

		if (!target_store.emplace(bin_int, 1).second) target_store[bin_int]++;
	}
}


void count_reverse_gRNA(std::unordered_map<unsigned, unsigned>& target_store, const std::string& gRNA) {
	// check if there N's or god forbid other characters in our target sequence
	if (gRNA.find_first_not_of("ATCG") == std::string::npos) {
		//is nice DNA =)
		//then we need to do a reverse binary encoding
		const unsigned bin_int = dna_to_revcompl_bin(gRNA);

		if (!target_store.emplace(bin_int, 1).second) target_store[bin_int]++;
	}
}


void find_mad7_targets(std::unordered_map<unsigned, unsigned>& target_store, const std::string& fasta_string) {
	// Checks for guideRNA targets in a (genomic) DNA fasta string and saves an encoded version to target_store
	unsigned size_to_save = 13;			//how many bp's should be used in a target
	// MAD7 PAM is YTTN (Y=T/C)
	// I'm gonna search for both TTT and CTT
	const std::string PAM_fwd_1 = "TTT";
	const std::string PAM_fwd_2 = "CTT";
	const std::string PAM_rev_1 = "AAA";
	const std::string PAM_rev_2 = "AAG";

	if (debug==1) std::cout << "reading" << std::endl << fasta_string << std::endl;

	std::string::size_type pos = -1;
	while ((pos=fasta_string.find(PAM_fwd_1, pos+1)) != std::string::npos) {
		//unlike Cas9 gRNA, the target is downstream of PAM
		if (pos<=fasta_string.length()-size_to_save-4) {
			const std::string gRNA = fasta_string.substr(pos+4,size_to_save);
			if (debug==1) std::cout << "found:" << gRNA << std::endl;

			count_forward_gRNA(target_store, gRNA);
		}
	}

	//now second fwd PAM
	pos = -1;
	while ((pos=fasta_string.find(PAM_fwd_2, pos+1)) != std::string::npos) {
		//unlike Cas9 gRNA, the target is downstream of PAM
		if (pos<=fasta_string.length()-size_to_save-4) {
			const std::string gRNA = fasta_string.substr(pos+4,size_to_save);
			if (debug==1) std::cout << "found:" << gRNA << std::endl;

			count_forward_gRNA(target_store, gRNA);
		}
	}

	//now opposite strand
	pos = -1;
	while ((pos=fasta_string.find(PAM_rev_1, pos+1)) != std::string::npos) {
		//found a PAM at pos
		if (pos>=size_to_save+1) {
			const std::string gRNA = fasta_string.substr(pos-size_to_save-1,size_to_save);
			if (debug==1) std::cout << "found_Rev:" << gRNA << std::endl;

			count_reverse_gRNA(target_store, gRNA);
		}
	}

	//now opposite strand second PAM
	pos = -1;
	while ((pos=fasta_string.find(PAM_rev_2, pos+1)) != std::string::npos) {
		//found a PAM at pos
		if (pos>=size_to_save+1) {
			const std::string gRNA = fasta_string.substr(pos-size_to_save-1,size_to_save);
			if (debug==1) std::cout << "found_rev:" << gRNA << std::endl;

			count_reverse_gRNA(target_store, gRNA);
		}
	}
}

void find_cas9_targets(std::unordered_map<unsigned, unsigned>& target_store, const std::string& fasta_string) {
	// Checks for guideRNA targets in a (genomic) DNA fasta string and saves an encoded version to target_store
	unsigned size_to_save = 13;			//how many bp's should be used in a target
	const std::string PAM_fwd = "GG";	//PAM is NGG
	const std::string PAM_rev = "CC";	//Reverse...should make a function

	if (debug==1) std::cout << "reading" << std::endl << fasta_string << std::endl;
	std::string::size_type pos = -1;
	while ((pos=fasta_string.find(PAM_fwd, pos+1)) != std::string::npos) {
		//found a PAM at pos
		if (pos>=size_to_save+1) {
			const std::string gRNA = fasta_string.substr(pos-size_to_save-1,size_to_save);
			if (debug==1) std::cout << "found:" << gRNA << std::endl;

			count_forward_gRNA(target_store, gRNA);
		}
	}

	//now one more time with the opposite strand
	pos = -1;
	while ((pos=fasta_string.find(PAM_rev, pos+1)) != std::string::npos) {
		//found a PAM at pos
		if (pos<=fasta_string.length()-size_to_save-3) {
			const std::string gRNA = fasta_string.substr(pos+3,size_to_save);
			if (debug==1) std::cout << "found-rev:" << gRNA << std::endl;

			count_reverse_gRNA(target_store, gRNA);
		}
	}
}

void count_cmd_line(std::vector<std::string> args) {
	std::string prg_name = args[0].substr(args[0].find_last_of('/')+1);
	std::cerr << "Usage: " << prg_name << " count <enzyme> <fasta>" << std::endl
			  << "    enzyme: Options are: Cas9 or Mad7 (case sensitive)" << std::endl
			  << "    fasta: File path for fasta file in which to find targets." << std::endl
			  << "Example:" << std::endl
			  << prg_name << " count Mad7 /my/path/genome.fa" << std::endl;
}

int count(const std::vector<std::string> &args) {
	std::cout << std::endl << "<<< welcome to PAM counting >>>" << std::endl;
	//command line checking
	if (args.size()!= 4) {
		std::cout << std::endl << "wrong number of arguments!" << std::endl;
		count_cmd_line(args);
		return 1;
	}

	auto enzyme_arg = args[2];
	make_uppercase(enzyme_arg);
	const auto& fasta_arg = args[3];

	std::unordered_set<std::string> implemented_enzymes = {"CAS9", "MAD7"};

	if (implemented_enzymes.count(enzyme_arg)==0) {
		std::cout << std::endl << std::endl << "Error: Unrecognized enzyme.\nEnzyme must be either Cas9 or Mad7 (case sensitive!).\nYou wrote: " << enzyme_arg << std::endl << std::endl;
		return 1;
	}
	//does file exist? 
	//no quality check is done on the fasta...we will attempt to parse it no matter what
	//could be unicode for all we know...that will crash a few things I think..
	//see https://stackoverflow.com/questions/313970/how-to-convert-stdstring-to-lower-case
	if(!does_file_exist(fasta_arg)){
		std::cout << std::endl << std::endl << "Error: Can't find fasta file.\nLast parameter must be path to fasta file.\nLast parameter was: " << fasta_arg << std::endl << std::endl;
		return 1;
	}

	std::cout << std::endl << std::endl << "!!!Great, Lets GO!!! 1" << std::endl;

	std::unordered_map<unsigned, unsigned> genome_targets;	//>target, exact copies>
	std::ifstream genome_fasta (fasta_arg);	//open the genome fasta

	if (genome_fasta.is_open()) {
		//BEGIN DB GENERATION
		std::cout << "\nfinding and counting all target sites." << std::endl;
		unsigned line_no = 0;
		std::string line;
		bool first = true;	//first line needs to be read different than the rest
		std::string fasta_string = "";

		while (std::getline(genome_fasta, line)) {
			if (line_no%1000000==0) std::cout << "Line no: " << line_no << "\r" << std::flush;
			//if (line_no>10000) break; //for debugging

			if (line[0]=='>') { //we have encountered the beginning of a sequence
				if (first) first = false; //is it the first time? then just ignore because we haven't collected any sequence yet
				else if (enzyme_arg == "CAS9") find_cas9_targets(genome_targets, fasta_string);
				else if (enzyme_arg == "MAD7") find_mad7_targets(genome_targets, fasta_string);

				fasta_string = ""; 	//reset fasta string
			} else {
				make_uppercase(line);
				fasta_string.append(line);
			}

			++line_no;
		}

		//the last fasta doesn't end with a > and so must be processed here after the last line of the file
		if (enzyme_arg == "CAS9") find_cas9_targets(genome_targets, fasta_string);
		else if (enzyme_arg == "MAD7") find_mad7_targets(genome_targets, fasta_string);
		std::cout << "Line no: " << line_no << std::endl;

		//END db generation
	} else {
		std::cout << std::endl << "Couldn't open fasta file...weird...it seems to be there..." << std::endl;
		return 1;
	}

	std::cout << "unique targets found: " << genome_targets.size() << std::endl;
	//close down the fasta file
	genome_fasta.close();

	//write file
	std::string filename;
	if (enzyme_arg == "CAS9") filename = "cas9_counts.laeb";
	else if (enzyme_arg == "MAD7") filename = "mad7_counts.laeb";

	std::ofstream counts (filename, std::ofstream::binary);

	// Record enzyme for which PAMs were collected
	counts.write(enzyme_arg.c_str(), enzyme_arg.size() + 1);

	// Add zero-padding to ensure that counts are aligned,
	// and for backwards compatibility
	const size_t blocksize = 2 * sizeof(unsigned);
	for (size_t i = enzyme_arg.size() + 1; i < blocksize; ++i) {
		counts.put(0);
	}

	for (auto &tal: genome_targets) {
		counts.write(reinterpret_cast<const char *>(&tal.first), sizeof(tal.first));
		counts.write(reinterpret_cast<const char *>(&tal.second), sizeof(tal.second));
	}
	counts.close();

	return 0;
}
