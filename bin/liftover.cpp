#define PROGRAM "liftover"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <utility>
#include <cstring>
#include <string.h>
#include <getopt.h>
#include <cassert>
#include <htslib/hts.h> 
#include <htslib/sam.h>

#define PROGRAM "liftover"

static const char VERSION_MESSAGE[] =
	PROGRAM " Version 0.0.0 \n";

static const char USAGE_MESSAGE[]=
PROGRAM " Version 0.0.0 \n"
"Usage: liftover [OPTIONS] \n"
"\n"
"Requires: BED file containing breaktigs and BAM alignment file in SAM or BAM format.\n"
"Outputs a SAM file of liftover coordinates in alignment file. Can pipe into samtools to conver into BAM file.\n"
"\n"
" Options:\n"
"	-t		threads (default is 1); Not yet implemented\n"
"	-b		breaktig bed file\n"
"	-a		alignment file in SAM or BAM format\n";

using namespace std;

namespace opt {
	unsigned threads=1;
	string breaktigfile;
	string alignmentfile;
	string outputfile="liftover-output";
	string r="rb";
	string w="wb";
}

static const char shortopts[] = "t:b:a:";

enum {OPT_HELP = 1, OPT_VERSION};

static const struct option longopts[]= {
	{"threads",required_argument, NULL, 't'},
	{"tigmint-output", required_argument, NULL, 'b'},
	{"alignment-file", required_argument, NULL, 'a'},
	{"help", no_argument, NULL, OPT_HELP},
	{"version", no_argument, NULL, OPT_VERSION},
	{NULL, 0, NULL, 0}
};

struct newtig {
	string newname; 
	int startpos; 
	int endpos; 
	int refid;
};

struct tiginfo {
	std::map<string, std::vector<newtig>> tig;
	std::map<string, std::pair<int, int>> unchangedtig; 
	std::unordered_set<string> changed;
	int numbrokentigs;
};


tiginfo parseTig(string breaktigfile) {
	ifstream tigFile; 
	tigFile.open(breaktigfile.c_str()); 
	
	tiginfo Tigs;
	Tigs.numbrokentigs = 0;

	string line;
	while (getline(tigFile, line)) {
		std::istringstream iss(line); 
		string part; 
		std::vector<string> pl;
		while (getline(iss, part, '\t')){
			pl.push_back(part);
		}
		if (pl[0] != pl[3]) {
			newtig data;
			data.newname = pl[3];
			data.startpos = stoi(pl[1]);
			data.endpos = stoi(pl[2]);
			data.refid = -1;
			if (Tigs.changed.find(pl[0]) == Tigs.changed.end()) {
				Tigs.changed.insert(pl[0]);
			}
			Tigs.tig[pl[0]].push_back(data); 
			Tigs.numbrokentigs++; 
		} else {
			std::pair<int, int> data; 
			data = std::make_pair(stoi(pl[1]), stoi(pl[2])); 
			Tigs.unchangedtig[pl[0]] = data; 
		}
	}
	 
	return Tigs;
}
	
void liftoverReads(string bwa_output, string outfile_name, string r, string  w, tiginfo Tigs) {

	samFile *bwafile = hts_open(bwa_output.c_str(), "r");
	htsFile *output = hts_open("-", "w"); //writes to stdout
	
	// Manipulate the header to remove the old and 
	bam_hdr_t *bwaHdr = sam_hdr_read(bwafile);

	std::unordered_map<string, int> newHdrDict; //to store all the new refids	

	// CHANGING THE HEADER
	bam_hdr_t *newHdr;	
	newHdr = bam_hdr_init();
	newHdr->sdict = NULL;

       	char *origtxt = bwaHdr->text;	
	char *HD = strtok(origtxt, "\n"); 
	std::cout << HD << endl; 
	int headerlen = bwaHdr->n_targets - Tigs.changed.size() + Tigs.numbrokentigs; 

 	newHdr->n_targets = headerlen; 
	newHdr->target_len = (uint32_t *)malloc(sizeof(uint32_t) * newHdr->n_targets); 
	newHdr->target_name = (char**) malloc(sizeof(char*) * newHdr->n_targets); 
	
	std::unordered_set<string> changed_unbroken; 

	int index = 0;
	for (int i=0; i<bwaHdr->n_targets; i++) {
		char *origname = bwaHdr->target_name[i]; 
		char  *SQ = strtok(NULL, "\n");
		if (Tigs.changed.find(origname) == Tigs.changed.end()) {
			newHdr->target_name[index] = origname;
			int bedLength = Tigs.unchangedtig[origname].second - Tigs.unchangedtig[origname].first;	
			if (bwaHdr->target_len[i] == bedLength) {
				newHdr->target_len[index] = bwaHdr->target_len[i];
				std::cout << SQ << endl;
			} else {
				newHdr->target_len[index] = bedLength; 
				std::cout << "@SQ\tSN:" << origname << "\tLN:" << bedLength << endl; 
				changed_unbroken.insert(origname); 
			}
			newHdrDict[origname] = index;  
		        index++; 
		}
	}
	
	for (auto &broken: Tigs.changed) {
		std::vector<newtig> pieces = Tigs.tig[broken];
	       	int numpieces = 0;	
		for (auto &piece: pieces) {
			newHdr->target_name[index] = strdup(piece.newname.c_str());
			newHdr->target_len[index] = piece.endpos - piece.startpos; 
			Tigs.tig[broken][numpieces].refid = index; 
			std::cout << "@SQ\tSN:" << piece.newname << "\tLN:" << (piece.endpos-piece.startpos) << endl; 
			newHdrDict[piece.newname] = index; 
			index++; 
		}
	}
	
	std::cout << "@PG\tID:cliftover\tPN:cliftover\tCL: see /projects/btl_scratch/jzhang/liftover/src/liftover --help" << endl; 
	


	bam1_t *aln = bam_init1();
	while(sam_read1(bwafile, bwaHdr, aln) >= 0 ) {

		char *contig = bwaHdr->target_name[aln->core.tid];	
		int32_t startpos = aln->core.pos; 
		if (Tigs.changed.find(contig) != Tigs.changed.end()) {
			std::vector<newtig> segments = Tigs.tig[contig]; 
			for (auto &segment : segments) {
				if (startpos >= segment.startpos && startpos < segment.endpos) {
					aln->core.tid = newHdrDict[segment.newname]; 
					aln->core.pos = startpos-segment.startpos; 
					break; 
				}
			}
		} else if (changed_unbroken.find(contig) != changed_unbroken.end()){
			aln->core.tid = newHdrDict[contig]; 
			aln->core.pos = startpos-Tigs.unchangedtig[contig].first;
		} else {
			aln->core.tid = newHdrDict[contig]; 
		}

		char *matecontig = bwaHdr->target_name[aln->core.mtid];
		int32_t mstartpos = aln->core.mpos;

		if (Tigs.changed.find(matecontig) != Tigs.changed.end()) {
			std::vector<newtig> segments = Tigs.tig[matecontig]; 
			for (auto &segment: segments) {
				if (mstartpos >= segment.startpos && mstartpos < segment.endpos) {
					aln->core.mtid = newHdrDict[segment.newname]; 
					aln->core.mpos = mstartpos-segment.startpos;
					break;
				}
			}
		} else if (changed_unbroken.find(matecontig) != changed_unbroken.end()) {
			aln->core.mtid = newHdrDict[matecontig];
			aln->core.mpos = mstartpos-Tigs.unchangedtig[matecontig].first;
		} else {
			aln->core.mtid = newHdrDict[matecontig]; 
		}
			        	 
		if (sam_write1(output, newHdr, aln) < 0) {
			std::cerr << "Error writing reads" << endl; 
			exit(EXIT_FAILURE); 
		} 
	}
	bam_destroy1(aln); 
	sam_close(bwafile); 
	sam_close(output); 

}
	 

int main(int argc, char * const argv[]) {
	bool die = false; 
	for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
		std::istringstream arg(optarg != NULL ? optarg: ""); 
		switch(c) {
			case '?': 
				die = true; 
				break;
			case 't':
				arg >> opt::threads;
				break;
			case 'b':
				arg >> opt::breaktigfile;				
				break;
			case 'a':
				arg >> opt::alignmentfile;
				break;
			case OPT_HELP: 
				std::cerr << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cerr << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			std::cerr << PROGRAM ": invalid option: " << (char) c << optarg << "\n"; 
			exit (EXIT_FAILURE);
		}
	}

	std::ifstream bedfile(opt::breaktigfile);
	if (!bedfile.good()) {
	std::cerr << "Cannot open breaktig file" << opt::breaktigfile << "\n";
		die = true; 
	}

	std::ifstream alignfile(opt::alignmentfile);
	if (!alignfile.good()) {
		std::cerr << "Cannot open alignment file" << opt::alignmentfile << "\n"; 
		die = true; 
	}

	std::ostringstream outfilename;
	if (opt::outputfile == "") {
		outfilename << opt::alignmentfile << "_liftover";
	}
	opt::outputfile = outfilename.str(); //only the base of the name not the appendix

	if (die) {
		std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
		exit (EXIT_FAILURE);
	}

	//Extract all of the broken tigs from the breaktigfile provided by Tigment
	tiginfo brokenTigs = parseTig(opt::breaktigfile);
	
	//Liftover read contig names and coordinates 
	liftoverReads(opt::alignmentfile, opt::outputfile, opt::r, opt::w, brokenTigs);
	return 0;	
}
