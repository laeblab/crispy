#include <atomic>
#include <fstream>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <algorithm>

#include "crispy++.h"
#include "tools.h"

class offtarg {
public:
	unsigned target;	//the sequence (converted to unsigned)
	unsigned char seed;		//number of muts in seed (first 5 bp)
	unsigned char rest;		//number of muts in rest of seq
	offtarg (unsigned itarget, unsigned char iseed, unsigned char irest) {
		target = itarget;
		seed = iseed;
		rest = irest;
	}
};

void looper(std::vector<offtarg>& targets, const offtarg& target,unsigned depth, const unsigned& str_len, unsigned pos=0) {
	// the purpose of this function is to generate all combinations of mutations with "depth or less" simultaneous mutations
	// It needs to return a list of mutated sequences as well as info regarding how many of the muts are in seed (first 5 bp) and how many are in the rest for each mut seq

	const unsigned seed_size = 5; 		//the size of the seed region is 5 bp
	const unsigned max_seed_muts = 2;	//throw away anything with more than 2 seed muts
	const unsigned max_total = 4;		//throw away anything with more than 4 total mutations
					//HMM max_total and depth sort of does the same thing...

	offtarg mut_target(0,0,0);	//init a resuable instance of class offtarg

	for (;pos < str_len; pos++) {
		//reset and update which pos we are mutating
		unsigned char seed_muts = target.seed;	//reset seed
		unsigned char rest_muts = target.rest;	//reset rest
		if (pos<seed_size) {++seed_muts;} else {++rest_muts;}	//update
		if (seed_muts>max_seed_muts) continue;
		if (seed_muts+rest_muts>max_total) continue;

		//identify base at current pos
		switch((target.target >> 2*pos) & 3) {
		case 0:
			//its A
			//mut to T
			mut_target = offtarg(target.target + (1<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			//mut to C
			mut_target = offtarg(target.target + (2<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			//mut to G
			mut_target = offtarg(target.target + (3<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			break;
		case 1:
			//its T
			//mut to A
			mut_target = offtarg(target.target - (1<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			//mut to C
			mut_target = offtarg(target.target + (1<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			//mut to G
			mut_target = offtarg(target.target + (2<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			break;
		case 2:
			//its C
			//mut to A
			mut_target = offtarg(target.target - (2<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			//mut to T
			mut_target = offtarg(target.target - (1<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			//mut to G
			mut_target = offtarg(target.target + (1<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			break;
		case 3:
			//its G
			//mut to A
			mut_target = offtarg(target.target - (3<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			//mut to T
			mut_target = offtarg(target.target - (2<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			//mut to C
			mut_target = offtarg(target.target - (1<<2*pos),seed_muts,rest_muts);
			targets.emplace_back(mut_target);
			if (depth!=0) looper(targets, mut_target, depth-1,str_len, pos+1);
			break;
		}
	}
}

void scoring_thread(std::vector<std::pair<std::string, unsigned long long>>& target_data_ref, std::unordered_map<unsigned, unsigned>& count_data_ref, const std::vector<std::vector<unsigned>>& scoring_matrix_ref, const unsigned& size_to_save_ref, std::vector<unsigned>& progress_ref, const unsigned pots, const unsigned ti, const std::string& enzyme) {
	//step 1: determine target range for this thread
	const unsigned start_pos = pots * ti;
	const unsigned end_pos = std::min<unsigned>(start_pos + pots, target_data_ref.size());

	//step 2: iterate over that part of target_data and score it
	for (unsigned current_pos = start_pos, i=1; current_pos<end_pos; ++current_pos, ++i) {
		//start with a clean slate
		unsigned long long score = 0;
		//convert to uppercase
		std::string target = split(target_data_ref[current_pos].first,'\t')[0]; //save string
		std::transform(target.begin(), target.end(), target.begin(), ::toupper);
		//check if target is nice DNA
		if (target.find_first_not_of("ATCG") == std::string::npos) {
			//it was, so assume target
			//convert target to binary
			std::string part_target;
			if (enzyme == "CAS9") {
				part_target = target.substr(20-size_to_save_ref,size_to_save_ref);
			} else if (enzyme == "MAD7") {
				part_target = target.substr(4, size_to_save_ref);
			} else {
				continue;  // todo: shouldn't it fail here?
			}

			const unsigned targ_seq = dna_to_bin(part_target);

			//try to find the target in the count data
			std::unordered_map<unsigned, unsigned>::const_iterator target_in_count = count_data_ref.find(targ_seq);
			if (target_in_count != count_data_ref.end()) {
				//add wt count * penalty for wt
				score += target_in_count->second * scoring_matrix_ref[0][0];
			}
			// create an offtarget template based on current target
			offtarg test_this(targ_seq,0,0);
			// create offtarget list to save generated offtargets in
			std::vector<offtarg> off_targets;
			//generate list of all offtargets
			looper(off_targets, test_this,2, size_to_save_ref);
			//try to find each offtarget in count data and if found, add their counts to the score * penalty for specific mutations
			for (auto& offtargs : off_targets) {
				target_in_count = count_data_ref.find (offtargs.target);
				if (target_in_count != count_data_ref.end())
					score += target_in_count->second * scoring_matrix_ref[offtargs.seed][offtargs.rest];
			}
			//save score
			target_data_ref[current_pos].second = score;
		}
		//update progress
		progress_ref[ti]=i;
	}
}


void thread_watch(const std::vector<unsigned>& progress, const std::atomic<bool>& notdone, const unsigned pots) {
	//print current score_thread status, update each sec
	const char esc = 27;

	std::cout << "....~~~~~≃≃≃≃≃========≃≃≃≃≃~~~~~...." << std::endl;
	for (bool first_loop = true; notdone; first_loop = false) {
		if (!first_loop) {
			std::this_thread::sleep_for(std::chrono::seconds(1));
			for (unsigned i = 0; i < progress.size() + 2; ++i) {
				std::cout << esc << "[1A\r";
			}
		}

		unsigned sum = 0;
		for (size_t i = 0; i < progress.size(); ++i) {
			const auto s = progress.at(i);
			std::cout << "Thread " << i << ": " << s << "/" << pots << "\n";

			sum += s;
		}

		std::cout << "....~~~~~≃≃≃≃≃========≃≃≃≃≃~~~~~....\nSum: " << sum << std::endl;
	}
}


void score_cmd_line(std::vector<std::string> args) {
	std::string prg_name = args[0].substr(args[0].find_last_of('/')+1);
	std::cerr << "Usage: " << prg_name << " score <number of threads> <file with targets to score> [count.laeb]" << std::endl
			<< "Target file is parsed line by line. Each line is split on tab and the first element is then checked for characters other than ATCG and if found, the line is ignored." << std::endl
			<< prg_name << " score 8 /my/path/targetList" << std::endl;
}

int score(const std::vector<std::string> &args) {
	std::cout << std::endl << "<<< welcome to target scoring >>>" << std::endl;
	//variables
	unsigned num_threads = 2;
	unsigned size_to_save = 13;
	std::string count_file_name = get_working_path() + "/counts.laeb";
	std::string target_file_name;
	std::unordered_map<unsigned, unsigned> count_data; //target, count
	std::vector<std::pair<std::string, unsigned long long>> target_data; //target, score

	//deal with cmd line
	if (args.size()<4 or args.size()>5) {
		std::cout << std::endl << "wrong number of arguments! There should be 4 or 5. You had " << args.size() << std::endl;
		score_cmd_line(args);
		return 1;
	}

	//check if num_threads is a number
	if (is_string_number(args[2])) {
		num_threads = stoul(args[2]);
		if (num_threads<2) num_threads = 2;
	} else {
		std::cout << std::endl << "parameter 2 is not a number! It should be the number of threads. E.g. 4. It was: '" << args[2] << "'" << std::endl;
		score_cmd_line(args);
		return 1;
	}

	//check if targetFile is a file
	if(does_file_exist(args[3])){
		target_file_name = args[3];
	} else {
		std::cout << std::endl << "parameter 3 is not a file. It was: '" << args[3] << "'" << std::endl;
		score_cmd_line(args);
		return 1;
	}

	//check count file
	if (args.size()==5) { //cmd line count filepath
		if(does_file_exist(args[4])){
			count_file_name = args[4];
		} else {
			std::cout << std::endl << "parameter 4 is not a file. It was: '" << args[4] << "'" << std::endl;
			score_cmd_line(args);
			return 1;
		}
	} else { //std count filepath
		if(!does_file_exist(count_file_name)){
			std::cout << std::endl << "No count index specified and couldn't find default count index, counts.laeb, in current directory! Please specify count index as last parameter." << std::endl;
			score_cmd_line(args);
			return 1;
		}
	}

	//good all parameters were fine
	//read count data
	std::cout << "\tcommand line good. Loading count data into memory..." << std::endl;
	std::ifstream count_file (count_file_name, std::ios::in|std::ios::binary|std::ios::ate);
	std::streampos size = count_file.tellg();
	char * memblock;
	memblock = new char [size];
	count_file.seekg(0, std::ios::beg);
	count_file.read(memblock, size);
	count_file.close();

	if ((unsigned)size < 2 * sizeof(unsigned)) {
		std::cerr << "ERROR: Invalid or truncated counts data file "
				  << "'" << count_file_name << "'" << std::endl;
		return 1;
	}

	// Ensure that header is NUL terminated
	memblock[2 * sizeof(unsigned) - 1] = '\0';
	const std::string enzyme = memblock;
	if (enzyme != "MAD7" && enzyme != "CAS9") {
		std::cerr << "ERROR: Invalid count data header in "
				  << "'" << count_file_name << "'" << std::endl;
		return 1;
	};
	std::cout << "\tCount data for " << enzyme << "..." << std::endl;

	count_data.reserve(size / (2 * sizeof(unsigned)));
	// The first block is skipped, as this contains the header
	for (unsigned i = 2*sizeof(unsigned); i < size; i+=2*sizeof(unsigned)) {
		int n = *(reinterpret_cast<int *>(memblock+i));
		int m = *(reinterpret_cast<int *>(memblock+i+sizeof(unsigned)));
		count_data.emplace(n,m);
	}

	//read target data
	std::cout << "\tCount data memorized. Loading target data into memory..." << std::endl;
	std::string line;
	std::ifstream target_file (target_file_name);
	while (getline(target_file, line)) {
		target_data.emplace_back(make_pair(line,0)); //a new target is born perfect
	}
	target_file.close();

	//split the_data into num_threads vectors
	std::cout << "\tTarget data memorized. Splitting it unto threads and scoring..." << std::endl;

	std::vector<std::vector<unsigned>> scoring_matrix(3, std::vector<unsigned>(5));
	// mismatch offtarget scoring matrix
	// seed mismatch on column and not-seed mismatch on row
	// 		0										1						2					3				4
	// 0	zero seed mm and zero non-seed mm		0 seed and 1 non-seed	0 seed 2 non-seed	0 seed 3 non	0 seed 4 non
	// 1	1 seed mm and 0 non-seed muts			1 seed and 1 non-seed	1 seed 2 non-seed	1 seed 3 non	not-a-target
	// 2	2 seed mm and 0 non-seed muts			2 seed and 1 non-seed	2 seed 2 non-seed	not-a-target	not-a-target

	// the values are going to be summed up.
	// E.g. for a specific gRNA there might just be 1 potential offtarget but it this example it has 0 seed mismatches and 0 additional mm's
	// so that gRNA gets an offtarget score of 500
	// another gRNA has 499 potential offtargets, but they are all with 2 seed and 2 additional mm's each giving a score of 1 for a total of
	// 499 offtarget score.
	// So in this example the second gRNA is deemed better than the first.
	// Depending on your use-case you might want to change this.
	// E.g. you could imagine you wanted to perform a test of gRNA edits at all offtargets and for that reason you wanted as few potential offtargets as possible.
	// Then you might set the score of all situations to 1 and then just pick the gRNA with the lowest score.
	scoring_matrix[0][0] = 500;
	scoring_matrix[0][1] = 100;
	scoring_matrix[0][2] = 50;
	scoring_matrix[0][3] = 20;
	scoring_matrix[0][4] = 3;
	scoring_matrix[1][0] = 80;
	scoring_matrix[1][1] = 30;
	scoring_matrix[1][2] = 15;
	scoring_matrix[1][3] = 2;
	scoring_matrix[2][0] = 20;
	scoring_matrix[2][1] = 5;
	scoring_matrix[2][2] = 1;

	//split data among threads
	const unsigned pots = (target_data.size() + num_threads - 1) / num_threads;

	//vector of threads
	std::vector<std::thread> mythreads;
	//thread progress counters
	std::vector<unsigned> thread_prog;
	//start each thread
	for (unsigned i=0;i<num_threads;i++) {
		thread_prog.emplace_back(0);
		mythreads.emplace_back(std::thread(scoring_thread,std::ref(target_data), std::ref(count_data), std::ref(scoring_matrix), std::ref(size_to_save), std::ref(thread_prog), pots, i, std::ref(enzyme)));
	}

	//start the thread watcher
	std::atomic<bool> notdone = true;
	std::thread pro(thread_watch, std::ref(thread_prog), std::ref(notdone), pots);
	std::cout << "threads running" << std::endl;

	//hold here until they all return
	for (auto& a_thread : mythreads) a_thread.join();
	//scoring complete
	notdone=false;
	//close thread watcher
	pro.join();

	//write scored targets to file
	std::cout << "Writing scores to file scored.laeb" << std::endl;
	std::ofstream counts ("scored.laeb");
	for (auto &score: target_data)
		if (score.first[0] == '>') {
			counts << score.first << std::endl;
		} else {
			counts << score.first << "\t" << score.second << std::endl;
		}
	counts.close();
	return 0;
}
