#ifndef UTILS_CPP_
#define UTILS_CPP_

#include <iostream>
#include <vector>
#include <string> 
#include <map>
#include "math.h"
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include "GenomicRange.h"
#include "kfunc.h"
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BFC.h"

#define TWO_BIT_MASK (3)
#define BITS_PER_BYTE (8)
#define BIG_ENOUGH (1024)

using namespace std;
using namespace SeqLib;

float poisson_pmf(int, double);

//std::map<std::string, GenomicRange> chomosomeDict( std::string bamFile) {
//	std::map<std::string, GenomicRange> mapOfChromosomes;
	


double inline getAlleleFrequency(const std::string& bamFile, const std::string& chr, int start, int end, int readSupport) {
    std::string search_area = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
    int regionLength = end - start + 1; // Adjusted for potential off-by-one error

    BamReader br;
    if (!br.Open(bamFile)) {
        // Handle error
        return -1; // Or throw an exception
    }
    BamRecord r;

    int count = 0;
    GenomicRegion gr(search_area, br.Header());
	br.SetRegion(gr);

    while (br.GetNextRecord(r)) {
        count++;
    }

    double alleleFrequency = 0.00;
    if (count > 0) {
        alleleFrequency = count / regionLength; // Cast to double for floating-point division
    }

	cout << bamFile << " " << count << " " << regionLength << " "<< alleleFrequency << endl;
    return alleleFrequency;
}
	

//######################################### 
double inline getKmerDiversity(std::vector<std::string>& reads) {
    std::vector<double> div_vec;
    const int defaultKSize = 5;

    for (const auto& read : reads) {
        std::set<std::string> hash;
        int kSize = (read.length() < defaultKSize) ? 3 : defaultKSize;

        for (size_t j = 0; j <= read.length() - kSize; ++j) {
            std::string kmer = read.substr(j, kSize);
            hash.insert(kmer);
        }

        int num_kmers = hash.size();
        int possible_kmers = (read.length() >= kSize) ? read.length() - kSize + 1 : 1;
        double kdiv = static_cast<double>(num_kmers) / possible_kmers;     
        div_vec.push_back(kdiv);
    }
    return std::accumulate(div_vec.begin(), div_vec.end(), 0.0) / div_vec.size();
}

//######################################### Calculate median
float inline computeGC ( std::string chr, int start, int end, RefGenome& reference ) {

	std::string sequence = reference.QueryRegion(chr, start, end);
	int count = 0;
	for (auto& ntd : sequence) {
		if (ntd == 'C' || ntd == 'c' || ntd == 'G' || ntd == 'g' ) {
			count++;
		}
	}
	float gc_content = (count/ (float) sequence.length()) * (float) 100;
	return gc_content;
}


//######################################### Calculate median
void inline printProgBar( int percent ) {
  string bar;

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  cout<< "\r " "[" << bar << "] ";
  cout.width( 3 );
  cout<< percent << "%     " << std::flush;
}

//######################################### Calculate median
 std::tuple<double, double> inline reciprocalOverlap ( int& start_A, int& end_A, int& start_B, int& end_B ) {

	double olap_a;
	double olap_b;

	int length_a = end_A - start_A;
	int length_b = end_B - start_B;

	// case1: exact overlap
	if ( start_A == start_B && end_A == end_B ) {
		olap_a = 1.00;
		olap_b = 1.00;
	}
	// case2: B is inside A
	if (  (start_B >= start_B && end_A > end_B)  || (start_B > start_B && end_A >= end_B )) {
		olap_a = 1.00 - (((start_B-start_A) + (end_A-end_B) ) / (double) length_a) ;
		olap_b = 1.00;
	}
	// case3: B is adjacent to A
	if (start_A <= start_B && end_B > end_A) {
		olap_a = 1.00 - ( (start_B - start_A )/ (double) length_a);
		olap_b = 1.00 - ( (end_B - end_A)/ (double) length_b);
		
		if (olap_a < 0) { olap_a = 0; }
		if (olap_b < 0) { olap_b = 0; }
	}
	//cout << start_A << "-" << end_A << "\t" << start_B << "-" << end_B << "\tlength" << olap_a << "," << olap_b << endl;

	return std::make_tuple(olap_a, olap_b);
 }

//######################################### Calculate median
string inline fillQualString ( std::string& seq ) {

	string outQual;

	for (int i = 0; i < seq.length(); i++) {
		outQual+="I";
	}
	return outQual;
}

//######################################### Calculate median
int inline calculateMedian ( vector<int>& Array ) {

   int size = Array.size();
   int median;
   if (size  % 2 == 0)
   {
      median = (Array[size / 2 - 1] + Array[size / 2]) / 2;
   }
   else 
   {
      median = Array[size / 2];
   }

   return median;
}

//######################################### return 1 is chr is somatic; 0 if no
int inline isChrSomatic ( string chrom ) {

	if (chrom == "23") { 
		return 0;
	}
	if (chrom == "chr23") {
		return 0;
	}
	if (chrom == "chrX") {
		return 0;
	}
	if (chrom == "X") {
		return 0;
	}
	return 1;
}

//######################################### return 1 is chr is somatic; 0 if no
int inline isChrX ( string chrom ) {

	if (chrom == "23") { 
		return 1;
	}
	if (chrom == "chr23") {
		return 1;
	}
	if (chrom == "chrX") {
		return 1;
	}
	if (chrom == "X") {
		return 1;
	}
	return 0;
}

//######################################### return chromosme lexicographical format
string inline returnChromLexicoGraphical ( string chrom ) {

	if (chrom == "23") { 
		chrom = "chrX";
	}
	if (chrom == "24") { 
		chrom = "chrY";
	}
	if (chrom == "chr23") {
		chrom = "chrX";
	}
	if (chrom == "chr24") {
		chrom = "chrY";
	}
	if (chrom == "chrM") {
		chrom = "chrM";
	}
	if (chrom == "M") {
		chrom = "chrM";
	}
	if (chrom == "MT") {
		chrom = "chrM";
	}
	return chrom;
}

//######################################### Binary Search
 bool inline binarySearch2 ( vector<GenomicRange>& vec, int n, string chr, int start, int end ) {

	int high   = n-1;
	int low    = 0;
	int result = -1;

	int iteration = 0;

	while (low <= high) {
		iteration++;
		int mid = (low+high)/2;
		
		if ( (start >= vec[mid].chromosome_start && end <= vec[mid].chromosome_end) ) {
			high = mid -1;
			result = mid;
		}
		else if (vec[mid].chromosome_start > start) {
			high = mid -1;
		}
		else if (vec[mid].chromosome_start < start) {
        	low = mid + 1;
		}
		
		if (iteration > 10000) {
			break;
		}
	}
	
	if (result == -1) {
		return 0;
	}
	else {
		return 1;
	}
}

//######################################### Binary Search
 vector<GenomicRange> inline binarySearch ( vector<GenomicRange>& vec, int n, string chr, int start, int end ) {

	int high   = n-1;
	int low    = 0;
	int result = -1;

	int iteration = 0;
    	std::vector<GenomicRange> intersected;
	std::map<std::string, GenomicRange> unique;

	while (low <= high) {
		iteration++;
		int mid = (low+high)/2;
		
		if ( (vec[mid].start >= start && vec[mid].end <= end) ) {

			high = mid -1;
			result = mid;
		}
		else if (vec[mid].start > start) {
			high = mid -1;
		}
		else if (vec[mid].start < start) {
        	low = mid + 1;
		}
		
		if (iteration > 10000) {
			break;
		}
	}
	
	if (result == -1) {
		return intersected;
	}
	int flag = 0;
	for (int i = result; i < n; i++) {
		if ( vec[i].start >= start && vec[i].end <= end ) {
			
			std::string key = vec[i].chromosomeA + std::to_string(vec[i].start) + std::to_string(vec[i].start) + vec[i].softclip_type;
			GenomicRange grange;
			grange.chromosomeA = vec[i].chromosomeA;
			grange.chromosomeB = vec[i].chromosomeB;
			grange.start = vec[i].start;
			grange.end   = vec[i].end;
			grange.softclip_type  = vec[i].softclip_type;
			unique.insert ( std::pair<std::string,GenomicRange>(key, grange) );
		}
		else {
			break;
		}
	}
	for (auto& r: unique) {
		GenomicRange grange;
		grange.chromosomeA = r.second.chromosomeA;
		grange.chromosomeB = r.second.chromosomeB;
		grange.start = r.second.start;
		grange.end   = r.second.end;
		grange.softclip_type  = r.second.softclip_type;
		intersected.push_back(grange);
	}

	return intersected;
}

//#########################################
bool inline liesOnBlackRegion( std::vector<GenomicRange>& badRegions, std::string& readChrom, int& readPos, int& readLen ) {

	bool doesOverlap = binarySearch2( badRegions, badRegions.size(), readChrom, readPos, readLen );

	if (!doesOverlap) {
		return 0;
	}
	else {
		return 1;
	}
}

// double getDP(const std::string& bamFile, const std::string& chr, int start, int end) {
//     BamReader br;
//     if (!br.Open(bamFile)) {
//         std::cerr << "Failed to open BAM file." << std::endl;
//         return 0.0;
//     }

//     // Define regions
//     int st_A = std::max(0, start - 10);
//     int st_B = start + 10;
//     int ed_A = end - 10;
//     int ed_B = end + 10;

//     std::vector<std::string> regions = {
//         chr + ":" + std::to_string(st_A) + "-" + std::to_string(st_B),
//         chr + ":" + std::to_string(ed_A) + "-" + std::to_string(ed_B)
//     };

//     int totalDepth = 0;
//     int regionLength = (st_B - st_A) + (ed_B - ed_A);
//     for (const auto& region : regions) {
//         GenomicRegion gr(region, br.Header());
//         br.SetRegion(gr);
//         BamRecord r;
//  		while (br.GetNextRecord(r)) {
// 	      	totalDepth++;
// 		}
//     }
//     br.Close();

//     // Calculate the average coverage depth
//     double averageDepth = regionLength > 0 ? static_cast<double>(totalDepth) / regionLength : 0.0;
//     return averageDepth;
// }

//#########################################
double inline getCoverage( std::string& bamFile, string chr, int start, int end) {

	int st_A = start-1;
	int st_B = start;
	int ed_A = end;
	int ed_B = end+1;
	if (st_A < 0) {
		st_A = ed_A;
		st_B = ed_B;
	}

	std::string upstream_area, downstream_area;

	upstream_area   = chr + ":" + std::to_string(st_A) + "-" + std::to_string(st_B);
	downstream_area = chr + ":" + std::to_string(ed_A) + "-" + std::to_string(ed_B);

	// cout << upstream_area << endl;
	// cout << downstream_area << endl;

	BamReader br;
	br.Open(bamFile);
	BamRecord r;

	std::vector<std::string> regions;
	regions.push_back(upstream_area);
	regions.push_back(downstream_area);
	int count = 0;
	for (auto& i : regions ) {
		GenomicRegion gr(i, br.Header());
		br.SetRegion(gr);
		while (br.GetNextRecord(r)) {
	      	count++;
		}
	}
	double cov;
	if (count == 0) {
		cov = 0.00;
	}
	else {	
		cov = count/2;
	}
	// cout << cov << endl;
	return cov;		
}
//#########################################
// Split function (works faster than boost implementation)
vector<string> inline split( std::string const& original, char separator ) {
    std::vector<std::string> results;
    std::string::const_iterator start = original.begin();
    std::string::const_iterator end = original.end();
    string::const_iterator next = std::find( start, end, separator );
    while ( next != end ) {
        results.push_back( string( start, next ) );
        start = next + 1;
        next = std::find( start, end, separator );
    }
    results.push_back( std::string( start, next ) );
    return results;
}

//#########################################
vector<string> inline splitSoftClip( string read ) {

 vector <string> tmp;
 tmp = split (read, '\t');
 string cigar    = tmp[7];
 string sequence = tmp[11];

 int softclipped, matched;
 string clipped_seq, matched_seq;
 vector<string> s;
 vector<string> m;
 if (tmp[0] == "LEFT") {
	s = split(cigar, 'S');
	m = split(s[1], 'M');

	softclipped = atoi(s[0].c_str());
	matched     = atoi(m[0].c_str());

	clipped_seq = sequence.substr(0,softclipped);
	matched_seq = sequence.substr(softclipped);
 }
 if (tmp[0] == "RIGHT") {
	m = split(cigar, 'M');
	s = split(m[1], 'S');

	matched     = atoi(m[0].c_str());
	softclipped = atoi(s[0].c_str());

	clipped_seq = sequence.substr(matched);
	matched_seq = sequence.substr(0, matched);
 }

 vector<string> seqs;//this vector will be returned 

 seqs.push_back(matched_seq);
 seqs.push_back(clipped_seq);

 return seqs;
}

//#########################################
int inline computeMedian (std::vector<int>& vec)  {

	std::sort(vec.begin(), vec.end());
	int size = vec.size();
	int median;
	if (size  % 2 == 0)
	{
	      median = (vec[size / 2 - 1] + vec[size / 2]) / 2;
	}
	else 
	{
	      median = vec[size / 2];
	}

	return median;
}


//#########################################
string inline revComp(string seq) {

	string rev_seq(seq);
	reverse(rev_seq.begin(),rev_seq.end());
	for (int i = 0; i < rev_seq.length(); i++) {
		switch (rev_seq[i]) 
		{
			case 'A': rev_seq[i] = 'T'; break;
			case 'T': rev_seq[i] = 'A'; break;
			case 'C': rev_seq[i] = 'G'; break;
			case 'G': rev_seq[i] = 'C'; break;
			case 'a': rev_seq[i] = 't'; break;
			case 't': rev_seq[i] = 'a'; break;
			case 'c': rev_seq[i] = 'g'; break;
			case 'g': rev_seq[i] = 'c'; break;
		}
	}
	return rev_seq; 
}

//#########################################
int inline mostFrequentPosition ( vector<int> positions ) {

	map<int, int> m;
	int maxCount = 0;
	int currentMax;
	int mostCommon = 0;
	for(int i=0;i < positions.size(); i++)
	{
	    int updated = m[positions[i]]++;  //Increment the value of key for counting occurances        
	    updated++; // due to post increment 
	    if (maxCount < updated) {
		 maxCount = updated;
		 currentMax = i;
		 mostCommon = positions[i];
	    }
	}
	return mostCommon;
}

//#########################################

double inline shannonEntropy ( std::string teststring ) {

   std::map<char , int> frequencies ;
   for ( char c : teststring )
     frequencies[ c ] ++ ;
   int numlen = teststring.length( ) ;
   double infocontent = 0 ;
   for ( std::pair<char , int> p : frequencies ) {
      double freq = static_cast<double>( p.second ) / numlen ;
      infocontent += freq * log2( freq ) ;
   }
   infocontent *= -1 ;

   return infocontent;
}
//#########################################

double inline computeDiscordantClusterSignificance ( int& cumulativeSizes, long& genomeSize, int& numInserts, int& nDiscordants  ) {
	// Based on BreakDancer's confidence interval
	double lambda = static_cast<double>(cumulativeSizes*numInserts)/genomeSize;
	double pvalue   = poisson_pmf(nDiscordants, lambda);

	return pvalue;
}

//#########################################
 std::tuple<double, double, double, int, int, int> inline calculatePvalueCoverage( std::string chr, int start , int end, std::string svtype, std::string bamFile)  {

	// This function calculates the significance (poisson dist) between inner-breakpoint counts and outter-breakpoint counts
	//         <- outer       inner ->   <- inner      outer ->
	//eg.    #############|_________________________|##############
	//       ######### break1 ################### break2 ##########

	BamReader br;
	br.Open(bamFile);
	BamRecord r;

	int size = end-start;

	if (size > 1000) {
		size = 1000;
	}

	int start_upstream = start - size;
	int end_upstream   = start;
	std::string upstream_area = chr + ":" + std::to_string(start_upstream) + "-" + std::to_string(end_upstream);

	int inner_start = start;
	int inner_end   = start + size;
	std::string inner_area = chr + ":" + std::to_string(inner_start) + "-" + std::to_string(inner_end);

	int start_downstream = end;
	int end_downstream   = end + size;
	std::string downstream_area = chr + ":" + std::to_string(start_downstream) + "-" + std::to_string(end_downstream);

	std::vector<std::string> regions;
	regions.push_back(upstream_area);
	regions.push_back(inner_area);
	regions.push_back(downstream_area);

    std::vector<int> v_counts;

	for (auto& i : regions ) {
	
		GenomicRegion gr(i, br.Header());
		br.SetRegion(gr);
		int count = 0;
		while (br.GetNextRecord(r)) {
		      if (r.MapQuality() < 20) {
				continue;
		      }
		      count++;
		}
		v_counts.push_back(count);
	}


   	double fisher_left_p, fisher_right_p, fisher_twosided_p;
   	kt_fisher_exact(v_counts[0], v_counts[1], v_counts[1], v_counts[2], &fisher_left_p, &fisher_right_p, &fisher_twosided_p);

	return std::make_tuple(fisher_left_p, fisher_right_p, fisher_twosided_p, v_counts[0], v_counts[1], v_counts[2] );
 }




#endif /* UTILS_CPP_ */
