/*
 * R package developed by:
 * Developer : David Mosen-Ansorena
 * All code (c)2012 David Mosen-Ansorena all rights reserved
 * all Credit goes to David Mosen-Ansorena, this is an adaptation as I couldn't get the R package to work.
 */

//libraries
#include <sstream>
#include <string>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
  //main variables
  int slide = 50;
  string inFile;
  string outFile;
  //get values from argv
  int opt;
  while((opt = getopt(argc, argv, "s:i:o:")) != -1) {
    switch(opt) {
      case 's':
        slide = atoi(optarg);
        break;
      case 'i':
        inFile = optarg;
        break;
      case 'o':
        outFile = optarg;
        break;
      default:
        abort();
      }
  }

  slide = slide * 1000;

  printf("\n sliding= %d \n in file= %s \n out file= %s \n", slide, inFile.c_str(), outFile.c_str());
  //chromosome definitions
  const int chr_num = 24;
  std::string chr_str[chr_num] = {"chr1","chr2","chr3","chr4","chr5","chr6",
								"chr7","chr8","chr9","chr10","chr11","chr12",
								"chr13","chr14","chr15","chr16","chr17","chr18",
								"chr19","chr20","chr21","chr22","chrX","chrY"};
  std::string chrnum_str[chr_num] = {"1","2","3","4","5","6",
                	              "7","8","9","10","11","12",
                								"13","14","15","16","17","18",
                								"19","20","21","22","X","Y"};
  double chr_max_len[chr_num] = {	275e+06,268e+06,218e+06,211e+06,200e+06,189e+06,176e+06,162e+06,156e+06,150e+06,
                              		149e+06,148e+06,127e+06,119e+06,113e+06,100e+06,90e+06,86e+06,66e+06,70e+06,53e+06,57e+06,
                              		171e+06,66e+06};
  //read type definitions
  const int types = 7;
	int flags[types][12] = { {0,16},{99,147,83,163},{81,161,97,145},{67,131,115,179},{65,129,113,177},{73,133,89,121,165,181,101,117,153,185,69,137},{-1} };
	int lengths[types] = {2,4,4,4,4,8,1};
  //vector initializations
	std::vector< std::vector< std::vector<unsigned int> > > counts;
	std::vector< std::vector<float> > gc_mean, mapq_mean;
	counts.resize(chr_num);
	gc_mean.resize(chr_num);
	mapq_mean.resize(chr_num);
	int chr_size;
	for (unsigned int chr=0; chr<counts.size(); chr++) {
		chr_size = (int) floor(chr_max_len[chr]/slide);
		gc_mean[chr].reserve(chr_size);
		mapq_mean[chr].reserve(chr_size);
		for (int j=0; j<chr_size; j++) {
			gc_mean[chr][j] = 0;
			mapq_mean[chr][j] = 0;
	}
	counts[chr].resize(types);
	for (unsigned int i=0; i<types-1; i++)
		counts[chr][i].reserve(chr_size);
	  }
  for (int chr=0; chr<chr_num; chr++){
  	for (unsigned int i=0; i<types-1; i++) {
  		for (unsigned int j=0; j<counts[chr][i].capacity(); j++) {
  			counts[chr][i][j] = 0;
      }
    }
  }

  char line[2048]; // should be greater than the expected longest line
	long linenum = 0;
	std::string sym = ".";
	int flag, chr, mapq;
	int is_paired = 0; // considered PEM if at least one improper read
	float gc, win_count;
	std::string seg;
	unsigned int index;
	unsigned int maxindex[chr_num];
	for (int i=0; i<chr_num; i++) {
		maxindex[i] = 0;
  }
  //setbuf(stdout, NULL);
  // set cmd
  std::string cmd = "samtools view ";
  cmd.append(inFile);
  FILE* pipe = popen(cmd.c_str(), "r");
  //loop over file and analyze
  next:
  while (fgets(line,sizeof(line),pipe)) {
    linenum++;
		if (linenum%200000 == 0) {
			printf(sym.c_str());
			if (linenum%10000000 == 0) {
				printf(" %dM reads processed.\r", linenum/1000000);
				if (sym == ".") sym = "|"; else sym = ".";
			}
		}
		std::stringstream ss(line);
		std::vector<std::string> elems;
		for (std::string each; std::getline(ss, each, '\t'); elems.push_back(each));
		if (elems.size() < 10) goto next; // if so, not a SAM read line
		chr = -1;
		for (int i=0;i<chr_num;i++)
			if (elems[2].compare(chr_str[i])==0 || elems[2].compare(chrnum_str[i])==0) {
				chr = i;
				goto nochr;
			}
		nochr:
		if (chr<0 || chr>=chr_num) goto next;	// skip if not in a chromosome of interest
		for (int i=0;i<types-1;i++) {
			flag = atoi(elems[1].c_str());
			for (int j=0;j<lengths[i];j++)
				if (flags[i][j]==flag) {
					index = (int)floor(atoi(elems[3].c_str())/slide);
					if (index<0) goto next;
					if (index > maxindex[chr])
						maxindex[chr] = index;
					if (i==1 || i==0) { // if the read is proper, update window GC and quality means
						seg = elems[9];
						gc = count(seg.begin(), seg.end(), 'G') + count(seg.begin(), seg.end(), 'C');
						gc = gc / seg.length();
						mapq = atoi(elems[4].c_str());
						win_count = (float)counts[chr][i][index];
						if (win_count==0) {
							gc_mean[chr][index] = gc;
							mapq_mean[chr][index] = mapq;
						} else {
							gc_mean[chr][index] *= win_count/(win_count+1);
							gc_mean[chr][index] += gc * 1/win_count;
							mapq_mean[chr][index] *= win_count/(win_count+1);
							mapq_mean[chr][index] += mapq * 1/win_count;
						}
					} else { // if the read is improper, reads are paired-end
						if (is_paired == 0)
							is_paired = 1;
					}
					counts[chr][i][index]++;
					goto next;
				}
		}
  }
  //save results
  printf("\nWriting results to disk...\n");
  std::ofstream myfile;
  myfile.open(outFile.c_str());
  myfile << "chrom" << "\t" << "win.start" << "\t" << "reads.gc" << "\t" << "reads.mapq";
  if (is_paired == 1) {
		for (int i=1;i<types-1;i++) {
			myfile << "\t" << "type" << i;
    }
  }
	else {
		myfile << "\t" << "counts";
  }
	myfile << "\n";
	for (int ch=0;ch<chr_num;ch++) {
		if ( ! (maxindex[ch] == 0 && gc_mean[ch][0] == 0 && mapq_mean[ch][0] == 0)) {
  		for (unsigned int w=0;w<=maxindex[ch];w++) {
  			myfile << chr_str[ch] << "\t" << w*slide << "\t" << gc_mean[ch][w] << "\t" << mapq_mean[ch][w];
  			if (is_paired == 1) {
  				for (int i=1;i<types-1;i++) {
  					myfile << "\t" << counts[ch][i][w];
          }
        }
  			else {
  				myfile << "\t" << counts[ch][0][w];
        }
  			myfile << "\n";
  		}
    }
	}
	myfile.close();
	printf("\nDone!\n");

  return 0;
}
