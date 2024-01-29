#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

// grapes headers
#include "Assembler.h"
#include "Aligner.h"
#include "alignContig.h"
#include "dna2bit.h"
#include "GenomicRange.h"
#include "sam_t.h"
#include "SV.h"
#include "Utils.cpp"
#include "vcf_t.h"
#include "varCall.h"
#include "smallSV.h"
#include "Trimmer.h"
#include "ssw_cpp.h"

// Seqlib headers
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BFC.h"
#include "SeqLib/FermiAssembler.h"
#include <boost/progress.hpp>

using namespace std;
using namespace SeqLib;


int getSupportingContigReads(std::vector<std::string>&reads, std::string& contig) {

	int supportingReads = 0;
	for(auto& read_seq: reads) {

		// Declares a default Aligner
		StripedSmithWaterman::Aligner ssw_aligner;

		// Declares a default filter
		StripedSmithWaterman::Filter filter;

		// Declares an alignment that stores the result
		StripedSmithWaterman::Alignment alignment;

		ssw_aligner.Align(read_seq.c_str(), contig.c_str(), contig.length(), filter, &alignment);
		int mismatches = alignment.mismatches;
		int q_begin = alignment.query_begin;
		int q_end = alignment.query_end;
		int q_size = q_end-q_begin;

		float minQuerySize = read_seq.length()*0.9;
		if (q_size >= minQuerySize && mismatches < 3) {
			supportingReads++;
		}
		else {

			// Declares a default Aligner
			StripedSmithWaterman::Aligner ssw_aligner;

			// Declares a default filter
			StripedSmithWaterman::Filter filter;

			// Declares an alignment that stores the result
			StripedSmithWaterman::Alignment alignment;

			std::string revcomp_seq = revComp(read_seq);

			ssw_aligner.Align(revcomp_seq.c_str(), contig.c_str(), contig.length(), filter, &alignment);
			int mismatches = alignment.mismatches;
			int q_begin = alignment.query_begin;
			int q_end = alignment.query_end;
			int q_size = q_end-q_begin;

			float minQuerySize = read_seq.length()*0.9;
			if (q_size >= minQuerySize && mismatches < 3) {
				supportingReads++;
			}

		}
	}
	return supportingReads;
}


void Analyze ( string& bamFile, vector<string>& contigs, string& chr, int& start, string& sctype, 
	int& totalReads, int& totalAssembled, string& refSequence, std::vector<vcf_t>& vcf_v, std::string& wes, double&kmerDiversity, double& meanMapQual ) {

	std::vector<sam_t> v_sam;

	for (int i = 0; i < contigs.size(); i++) {

		BamRecord r;
		BamRecordVector brv;
		bool hardclip = false;
		float secondary_cutoff = 0.90; // secondary alignments must have score >= 0.9*top_score
		int secondary_cap = 10; // max number of secondary alignments to return

		// Make an in-memory BWA-MEM index of region
		BWAWrapper bwa;
		UnalignedSequenceVector usv = {{chr, refSequence}};
		bwa.ConstructIndex(usv);
		string num = std::to_string(i);
		string Qname = chr + "_" + num;
		BamRecordVector results; // alignment results (can have multiple alignments)
		
		bwa.SetAScore(2);
		bwa.SetGapOpen(30);
		bwa.SetGapExtension(1);			
		bwa.SetBandwidth(1000);	
		bwa.SetMismatchPenalty(4);		
		bwa.SetReseedTrigger(1.5);
		bwa.Set3primeClippingPenalty(5);
		bwa.Set5primeClippingPenalty(5);
			
		BamHeader myHeader;
		myHeader = bwa.HeaderFromIndex();

		bwa.AlignSequence(contigs[i], Qname, results, hardclip, secondary_cutoff, secondary_cap);
		vector<string> bwa_results;

		for (auto& i : results) {

			sam_t sam_read;
			sam_read.read_name = i.Qname();
			sam_read.chr       = myHeader.IDtoName(i.ChrID());
			sam_read.pos       = i.Position();
			sam_read.align_pos = i.AlignmentPosition();
			sam_read.align_end = i.AlignmentEndPosition();
			sam_read.cigar     = i.CigarString();
			sam_read.strand    = i.ReverseFlag() == true ? '-' : '+';
			sam_read.order     = i.FirstFlag()   == true ? '1' : '2';
			int mapq_ascii     = i.MapQuality();
			sam_read.mapq      = mapq_ascii;
			sam_read.seq       = i.Sequence();
			sam_read.clipped   = i.NumClip();
	
			v_sam.push_back(sam_read);
		 }
	 }

	 SV call(v_sam, chr, chr, start, start, sctype, totalReads, totalAssembled, bamFile, refSequence, "small", 0, 0.00, kmerDiversity, meanMapQual);
	 std::string undetermined = "undetermined";
	 if (call.classify(vcf_v, undetermined) == 1) {

	 }
}

void smallSV::callSmallSV (map<string, vector<sam_t>>& breakReadClusters, std::map<std::string, std::vector<GenomicRange>>& breakReadLocations, std::string& fastq, std::vector<vcf_t>& vcf_v, std::string& wes, std::string& outdir) {

   // Goal here is to jointly consider all brek-reads that fall inside +/-500 bp from each soft-clipping position
   // as assembly candidates.

   int N = breakReadClusters.size();
   int count = 0;
	
   string dummy = "INIT";
   dna2bit DNA(dummy);

   RefGenome ref;
   ref.LoadIndex(reference);

   // opening bam file	
   BamReader reader;
   reader.Open(bamFile);

   BamHeader myHeader;
   myHeader = reader.Header();

   std::multimap<std::string, std::string> visitedBreak;

   std::string logF = outdir + "/" + "grapes.log";
   std::ofstream logFile;

   logFile.open(logF, std::ofstream::app);
   int minTrimSize = readLength*0.20;

	std::map<std::string, std::vector<std::string>> breaksInWindows;

	for (auto& breakC : breakReadClusters) {
		vector<string> tmpCoord = split(breakC.first, '\t');
		string chr = tmpCoord[0];

		int pos = std::stoi(tmpCoord[1]);		
		int start = pos - 250;
		int end = pos + 250;

		vector<GenomicRange> breakClusters = binarySearch(breakReadLocations[chr], 
			breakReadLocations[chr].size(), chr, start, end);

		for (auto& c : breakClusters) {
			std::string mkey = chr + "\t" + std::to_string(c.start) + "\t" + c.softclip_type;
			visitedBreak.insert(std::pair<std::string, std::string>(mkey, mkey));
			if (visitedBreak.count(mkey) <= 1) {
				breaksInWindows[breakC.first].push_back(mkey);
			}
		}
	}


  	boost::progress_display show_progress( breaksInWindows.size() );

	// Using dynamic schedule for assigning chunks on runtime. This works better when we have unbalanced tasks between threads.
	// http://jakascorner.com/blog/2016/06/omp-for-scheduling.html
	#pragma omp parallel 
	{
	    std::vector<vcf_t> vcf_v_private;
	
	    int tid = omp_get_thread_num();
	
		#pragma omp for schedule(dynamic, 1)

   		for (int i = 0; i < breaksInWindows.size(); i++) {

			auto breakC = breaksInWindows.begin();
       		advance(breakC, i);

			vector<string> tmpCoord = split ( breakC->first, '\t');
			string chr = tmpCoord[0];
			string coordinates = tmpCoord[0]+"_"+tmpCoord[1]+"_"+tmpCoord[2];

			int pos = std::stoi(tmpCoord[1]);		
			int start = pos - 500;
			int end = pos + 500;

			int totalReads;
			int totalReadsAssembled;
			int meanMapQual;
			int meanBaseQual;

			std::vector<std::string> breakReadVec;			
			std::vector<int> mapqQualities;
			std::vector<int> baseQualities;
			std::vector<string> contigsJoint;
			std::map<std::string, std::string> seenReadsMap;

			for (auto& subclusterKey: breakC->second) {

				vector<string> tmpClust = split( subclusterKey, '\t');			
				std::string mkey = tmpClust[0] + "\t" + tmpClust[1] + "\t" + tmpClust[2];
				auto it1 = breakReadClusters.find(mkey);

				if (it1 != breakReadClusters.end()) {
					for (auto& read: it1->second) {

						if (seenReadsMap.find(read.read_name) != seenReadsMap.end()) {
							continue;
						}

						string seq  = DNA.decompressDNA(read.bitSeq, read.seqLen);
						mapqQualities.push_back(read.mapq);
						baseQualities.push_back(read.mean_qual);

						// Perform a naive 5-end trimming
						Trimmer trim1(seq, read.qual);
						string trimmedseq = trim1.trim();

						if (trimmedseq.length() >= minTrimSize) {
							breakReadVec.push_back(trimmedseq);
						}
						seenReadsMap[read.read_name] = read.read_name;
					}
				}
			}
			if (breakReadVec.size() > 0) {

				meanMapQual = accumulate( mapqQualities.begin(), 
					mapqQualities.end(), 0.0)/mapqQualities.size();

				meanBaseQual = accumulate(baseQualities.begin(),
					baseQualities.end(), 0.0)/baseQualities.size();
				// cout << "Mean Base Qual " << " " << meanBaseQual << endl;

				int numAssembledA = 0;
				int cidx = 0;
				UnalignedSequenceVector unalignVect;
				for (auto&candidateRead : breakReadVec) {
					std::string Qual= "";
					for (auto&ntd : candidateRead) {
						Qual+= "E";
					}
					UnalignedSequence unSeq(std::to_string(cidx), candidateRead, Qual);
					unalignVect.push_back(unSeq);
					cidx++;
				}
				FermiAssembler f;

				// add the reads and error correct them  
				f.AddReads(unalignVect);
				f.CorrectReads();

				// // peform the assembly
				f.PerformAssembly();
				std::vector<std::string> contigs = f.GetContigs();
				contigsJoint = contigs;

				f.ClearReads();
				f.ClearContigs();

				totalReads = breakReadVec.size();
			}
			std::vector<sam_t> samAlignments;
			std::vector<vcf_t> smallSV;

			if (contigsJoint.size() > 0) {
				
				RefGenome ref;
 				ref.LoadIndex(reference);

        		int lengthChrom = myHeader.GetSequenceLength(chr) ;

				start = start >= 0 ? start : 0;
				end = end <= lengthChrom ? end : lengthChrom;
				std::string localReference = ref.QueryRegion(chr, start, end);

				BamRecord r;
				BamRecordVector results; // alignment results (can have multiple alignments)
				bool hardclip = false;
				float secondary_cutoff = 0.95; // secondary alignments must have score >= 0.9*top_score
				int secondary_cap = 5; // max number of secondary alignments to return

				// Make an in-memory BWA-MEM index of region
				BWAWrapper bwamem;
				UnalignedSequenceVector usv = {{chr, localReference}};
				bwamem.ConstructIndex(usv);
				bwamem.SetAScore(2);
				bwamem.SetGapOpen(8);
				bwamem.SetGapExtension(0.5);			
				bwamem.SetBandwidth(2000);	
				bwamem.SetMismatchPenalty(5);		
				bwamem.SetReseedTrigger(1.5);
				bwamem.Set3primeClippingPenalty(5);
				bwamem.Set5primeClippingPenalty(5);

				int idxContig = 0;

				double kmerDiversity = 0.0;
				if (contigsJoint.size() > 0) {
					kmerDiversity = getKmerDiversity(contigsJoint);
				}
				for (auto& contig : contigsJoint) {
					int totalReadsAssembled = getSupportingContigReads(breakReadVec, contig);
					string Qname = coordinates + "_" + std::to_string(idxContig);
					bwamem.AlignSequence(contig, Qname, results, hardclip, secondary_cutoff, secondary_cap);
					for (auto& i : results) {
						sam_t sam_aln;
						sam_aln.read_name = i.Qname();
						sam_aln.chr       = i.ChrName();
						sam_aln.pos       = i.Position();
						sam_aln.align_pos = i.AlignmentPosition();
						sam_aln.align_end = i.AlignmentEndPosition();
						sam_aln.cigar     = i.CigarString();
						sam_aln.strand    = i.ReverseFlag() == true ? '-' : '+';
						sam_aln.order     = i.FirstFlag() == true ? '1' : '2';
						int mapq_ascii    = i.MapQuality();
						sam_aln.mapq      = mapq_ascii;
						sam_aln.seq       = i.Sequence();
						sam_aln.mean_qual = meanBaseQual;
						samAlignments.push_back(sam_aln);
					}

					SV call(samAlignments, chr, chr, start, end, "undef",  
						totalReads, totalReadsAssembled, bamFile, localReference, 
						"small",  0, 1e-5, kmerDiversity, meanBaseQual);
					string svtype = "undetermined";

					vector<vcf_t> smallSV;
					smallSV = call.classifySmall(svtype);
					for (auto& b : smallSV) {
						vcf_v_private.push_back(b);
					}
					idxContig++;
				}		
			}
			#pragma omp critical 
			{
				++show_progress;
			}
			contigsJoint.clear();
			totalReadsAssembled = 0;
			totalReads = 0;
		}
	#pragma omp critical
		vcf_v.insert(vcf_v.end(), vcf_v_private.begin(), vcf_v_private.end());
    }
}

