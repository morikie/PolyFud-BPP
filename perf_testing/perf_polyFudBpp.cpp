#include <iostream>
#include <random>
#include <string>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
#include <boost/spirit/include/qi.hpp>
#include <seqan/seq_io.h>
#include "../src/polyFudBpp.hpp"
#include "../src/polarUtility.hpp"
#include "../src/refGeneParser.hpp"
#include "../src/polyFud.hpp"
#include "perf_polyFudBpp.hpp"


namespace fs = boost::filesystem;
namespace qi = boost::spirit::qi;

//structure needed for parsing the data sets
struct PasPosition {
	std::string seqId;
	std::string chr;
	size_t pos;
	std::string strand;

	void reset() {
		this->chr = std::string();
		this->pos = 0;
		this->strand = std::string();
		this->seqId = std::string();
	}
};


//used by boost::spirit to parse into a PasPosition object
BOOST_FUSION_ADAPT_STRUCT (
	PasPosition,
	(std::string, seqId)
	(std::string, chr)
	(size_t, pos)
	(std::string, strand)
)


PolyFudBpp::bppVectorPerTranscriptMap setNewBppDictionary(fs::path bppFile) {
	std::cout << "Initializing BPP map...";
	std::ifstream inFile(bppFile.string());
	std::string line;
	std::string transcriptId;
	PolyFudBpp::utrBppVector tempVec;
	PolyFudBpp::bppVectorPerTranscriptMap tempMap;
	size_t lineCount = 0;
	while (std::getline(inFile, line)) {
		switch (lineCount) {
		case 0:
		{	

			qi::parse(line.begin(), line.end(), 
				">" >> +qi::char_, 
				transcriptId);
			//std::cerr << transcriptId << std::endl;
			lineCount++;
			break;
		}
		case 1: 
		{	
			qi::parse(line.begin(), line.end(),
				 +(qi::double_ >> " "),
				 tempVec);
			if (tempMap.find(transcriptId) == tempMap.end()) {
				tempMap[transcriptId] = tempVec;
			}
			transcriptId.clear();
			tempVec.clear();
			lineCount = 0;
			break;
		}
		}
	}
	inFile.close();
	std::cout << "DONE" << std::endl;
	/*
	for (auto & pair : tempMap) {
		std::cerr << pair.first << ": ";
		for (auto & item : pair.second) {
			std::cerr << item << ", ";
		}
		std::cerr << std::endl;
	}
	*/
	return tempMap;	
}


/**
 * Parsing FASTA-like formats
 */
void getPositiveDataSet(std::vector<PasPosition> & returnVec, const fs::path & positiveFasta) {
	std::ifstream in(positiveFasta.string());
	std::string line;
	size_t lineCount = 0;
	PasPosition pasPos;

	while (std::getline(in, line)) {
		switch (lineCount) {
		case 0:
		{	
			qi::parse(line.begin(), line.end(), 
				">" >> +~qi::char_('.') >> qi::omit[*~qi::char_("|")] >> "|" 
				>> +~qi::char_('|') >> "|" 
				>> qi::uint_ >> "|" 
				>> +qi::char_, 
				pasPos);
			returnVec.push_back(pasPos);
			pasPos.reset();
			lineCount++;
			break;
		}
		case 1: 
		{
			lineCount = 0;
			break;
		}
		}
	}
	in.close();
}


bool testCase(RefGeneParser & parser) {
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		if (polar::utility::buildIndexFile(referenceGenome)) {
			seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str());
		} else {
			std::cerr << "could not open index file for " << referenceGenome << std::endl;
			return EXIT_FAILURE;
		}
	}
	std::string testId = "NM_004081";
	std::string testUtr = "taaattccgttgttactcaagatgactgcttcaagggtaaaagagtgcatcgctttagaagaagtttggcagtatttaaatctgttggatcctctcagctatctagtttcatgggaagttgctggttttgaatattaagctaaaagttttccactattacagaaattctgaattttggtaaatcacactgaaactttctgtataacttgtattattagactctctagttttatcttaacactgaaactgttcttcattagatgtttatttagaacctggttctgtgtttaatatatagtttaaagtaacaaataatcgagactgaaagaatgttaagatttatctgcaaggatttttaaaaaattgaaacttgcattttaagtgtttaaaagcaaatactgactttcaaaaaagtttttaaaacctgatttgaaagctaacaattttgatagtctgaacacaagcatttcacttctccaagaagtacctgtgaacagtacaatatttcagtattgagctttgcatttatgatttatctagaaatttacctcaaaagcagaatttttaaaactgcatttttaatcagtggaactcaatgtatagttagctttattgaagtcttatccaaacccagtaaaacagattctaagcaaacagtccaatcagtgagtcataatgtttattcaaagtattttatcttttatctagaatccacatatgtatgtccaatttgattgggatagtagttaggataactaaaattctgggcctaattttttaaagaatccaagacaaactaaactttactgggtatataaccttctcaatgagttaccattcttttttataaaaaaaattgttccttgaaatgctaaacttaatggctgtatgtgaaatttgcaaaatactggtattaaagaacgctgcagcttttttatgtcactcaaaggttaatcggagtatctgaaaggaattgtttttataaaaacattgaagtattagttacttgctataaatagatttttatttttgttttttagcctgttatatttccttctgtaaaataaaatatgtccagaagaggcatgttgtttctagattaggtagtgtcctcattttatattgtgaccacacagctagagcaccagagcccttttgctatactcacagtcttgttttcccagcctcttttactagtctttcaggaggtttgctcttagaactggtgatgtaaagaatggaagtagctgtatgagcagttcaaaggccaagccgtggaatggtagcaatgggatataatacctttctaagggaaacatttgtatcagtatcatttgatctgccatggacatgtgtttaaagtggctttctggcccttctttcaatggcttcttccctaaaacgtggagactctaagttaatgtcgttactatgggccatattactaatgcccactggggtctatgatttctcaaaattttcattcggaatccgaaggatacagtctttaaactttagaattcccaagaaggctttattacacctcagaaattgaaagcaccatgactttgtccattaaaaaattatccatagtttttttagtgcttttaacattccgacatacatcattctgtgattaaatctccagatttctgtaaatgatacctacattctaaagagttaattctaattattccgatatgaccttaaggaaaagtaaaggaataaatttttgtctttgttgaagtatttaatagagtaaggtaaagaagatattaagtccctttcaaaatggaaaattaattctaaactgagaaaaatgttcctactacctattgctgatactgtctttgcataaatgaataaaaataaactttttttcttcaaatgtg";
	std::string testOffset = "tttttggctttccgatgtaataatgtaaaatggtggggagttgcgtgggaactgtgtaacaaggtttaaattcgtataacaagctttagattcttaaaatgcagaagtataaagttcagtatactaatctgtctgagttagcccataaaagcaaatgtaggtacaaagataagtttaagaggtgcatcaacagcagtgcag";
	auto txPro = parser.getValueByKey(testId);
	std::string extractedUtr = polar::utility::getUtrSequence(txPro, faiIndex);
	
	size_t testLengthFromGenome = polar::utility::getUtrLength(txPro);
	//std::cerr << "utrLenFromGenome: " << testLengthFromGenome << " | utrLen NM_004081: " << testUtr.size() << std::endl;
	if (extractedUtr != testUtr) {
		std::cerr << "extractedUtr:" << std::endl << extractedUtr << std::endl << testId << std::endl << testUtr << std::endl;
		return false;
	}
	size_t zero = 25345239;
	size_t two = 25345237;
	size_t shouldBeZero = polar::utility::mapGenomePosToTxPos(txPro, zero);
	size_t shouldBeTwo = polar::utility::mapGenomePosToTxPos(txPro, two);
	if (shouldBeZero != 0 || shouldBeTwo != 2) {
		std::cerr << "shouldBeZero: " << shouldBeZero << std::endl;
		std::cerr << "shouldBeTwo: " << shouldBeTwo << std::endl;
		//return false;
	}
	
	PolyFudBpp tempObj = PolyFudBpp(testId);
	std::string tempObjUtr = tempObj.getUtrSeq();
	std::string tempObjOffsetSeq = tempObj.getOffsetSeq();
	if (testUtr != tempObjUtr) return false;
	if (tempObjOffsetSeq != testOffset) {
		std::cerr << testOffset << std::endl 
			<< tempObjOffsetSeq << std::endl;
		return false;
	}
	std::string testId2 = "NM_018836";
	std::string testUtr2 = "ctggccgaagtcttttttacctcctgggggcagggcagacgccgtgtgtctgtttcacggattccgttggtgaacctgtaaaaacaaaacaaacaaaacaaaacaaaaaagacaaaacctaaaactgagctatctaagggggagggtccccgcacctaccacttctgtttgccggtgggaaactcacagagcaggacgctctaggccaaatctatttttgtaaaaatgctcatgcctatgggtgactgccttctcccagagttttctttggagaacagaaagaagaaaggaaagaaaggaaccagaggcagagagacgaggatacccagcgaaagggacgggaggaagcatccgaaacctaggattcgtcctacgattctgaacctgtgccaataataccattatgtgccatgtactgacccgaaaggctcggccgcagagccggggcccagcgaatcacgcagagaaatcttacagaaaacaggggtgggaatctcttccgatagagtcgctatttctggttaatatacatatataaatatataaatacaaacacacacacacactttttttgtactgtagcaatttttgaagatcttaaatgttcctttttaaaaaaaagaattgtgttataggttacaaaatctgatttatttaacatgcttagtatgagcagaataaaccagtgttttctactttggcaactcacgtcacacacatattacacacatgtgcgcattacacacacacaatacacatacatgcatatagacgcatctattggaaatgcagttccacaggtgagcatgttctttctggtgacctggtattccatcaccattcaccccaggggacagcctcgaccgagacaaggaggcccttaaatgacagcctgcatttgctagacggttggtgagtggcatcaaatgtgtgacttactatcttgggccagaactaagaatgccaaggttttatatatgtgtgtgtatatatatatatatatatatatatatatatatatatatgtttgtgtgtgtatatatatatatatatatatatgtttgtgtgtgtatatatatgtttgtgtatatatatacacatatgcatacatatgatttttttttttcatttaagtgttggaagatgctacctaacagccacgttcacatttacgtagctggttgcttacaaacgggcctgagcccctggttgggtgggtggtggattcttggacgtgtgtgtcatacaagcatagactggattaaagaagttttccagttccaaaaattaaaggaatatatcctta";
	std::string testOffset2 = "tgatgtgtgtgtgtaatatcagggcagaacttagacatacgtgaagggccccggttggtttgaaaacgaaaaatagtcattctgtgtgcaaaccacaaggctgccccagtcaggcagcgccctgacctggcctgtgctgcattgccttcccttgcgcaggtgggcaggtgtggcccgctttttctagggcccaagggtgac";
	txPro = parser.getValueByKey(testId2);
	if (! (txPro == parser.emptyRefGeneProperties)) {
		extractedUtr = polar::utility::getUtrSequence(txPro, faiIndex);
	} else {
		return false;
	}
	
	zero = 4715104;
	two = 4715106;
	shouldBeZero = polar::utility::mapGenomePosToTxPos(txPro, zero);
	shouldBeTwo = polar::utility::mapGenomePosToTxPos(txPro, two);

	if (shouldBeZero != 0 || shouldBeTwo != 2) {
		std::cerr << "shouldBeZero: " << shouldBeZero << std::endl;
		std::cerr << "shouldBeTwo: " << shouldBeTwo << std::endl;
		return false;
	}
	PolyFudBpp tempObj2 = PolyFudBpp(testId2);
	std::string tempObj2Utr = tempObj2.getUtrSeq();
	std::string tempObj2OffsetSeq = tempObj2.getOffsetSeq();

	if (testUtr2 != tempObj2Utr) {
		std::cerr << testUtr2 << std::endl 
			<< tempObj2Utr << std::endl;
		return false;
	}
	if (tempObj2OffsetSeq != testOffset2) {
		std::cerr << testOffset2 << std::endl 
			<< tempObj2OffsetSeq << std::endl;
		return false;
	}
	
	testLengthFromGenome = polar::utility::getUtrLength(txPro);
	//std::cerr << "utrLenFromGenome: " << testLengthFromGenome << " | utrLen " << testId2 << ": " << testUtr2.size() << std::endl;
	if (extractedUtr != testUtr2) {
		std::cerr << "extractedUtr:" << std::endl << extractedUtr << std::endl << testId2 << std::endl << testUtr2 << std::endl;
		return false;
	}
	return true;
}


int main (int argc, char * argv[]) {
	fs::path bppOut = "../perf_testing/bppOutput.txt";
	fs::path bppOutTn = "../perf_testing/bppOutputTn.txt";
	fs::path positiveDataset = "../perf_testing/positiveSet.fa";
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	fs::path refGenPath = "ucsc_data/refGene.txt";
	fs::path bppDictTn = "../perf_testing/utrBppPerTranscriptTn.txt";	
	std::vector<PasPosition> pasVector;
	getPositiveDataSet(pasVector, positiveDataset);
	typedef size_t truePos;
	std::vector<std::pair<PolyFudBpp, truePos> > resVector;
	RefGeneParser refGen = RefGeneParser(refGenPath);		

	if (! testCase(refGen)) {
		std::cerr << "test failed!" << std::endl;
		return 1;		
	}
	std::unordered_map<PolyFud::motifSequence, unsigned int> countDiscardedPas = {
		{std::string("aataaa"), 0},
		{std::string("attaaa"), 0},
		{std::string("tataaa"), 0},
		{std::string("agtaaa"), 0},
		{std::string("aagaaa"), 0},
		{std::string("aatata"), 0},
		{std::string("aataca"), 0},
		{std::string("cataaa"), 0},
		{std::string("gataaa"), 0},
		{std::string("aatgaa"), 0},
		{std::string("actaaa"), 0},
		{std::string("aataga"), 0}
	};
	
	std::unordered_map<PolyFud::motifSequence, double> thresholdMap = {
		{std::string("aataaa"), 0.0},
		{std::string("attaaa"), 0.0},
		{std::string("tataaa"), 0.0},
		{std::string("agtaaa"), 0.0},
		{std::string("aagaaa"), 0.0},
		{std::string("aatata"), 0.0},
		{std::string("aataca"), 0.0},
		{std::string("cataaa"), 0.0},
		{std::string("gataaa"), 0.0},
		{std::string("aatgaa"), 0.0},
		{std::string("actaaa"), 0.0},
		{std::string("aataga"), 0.0}
	};
	PolyFud::setThresholdMap(thresholdMap);
	
	size_t ignoredPas = 0;	
	size_t loopCounter = 0;
	size_t evaluatedPas = 0;
	size_t maxUtrLength = 3000;
	size_t numTooLongSequences = 0;
	std::set<std::string> uniques;
	std::ofstream bppOutStream(bppOut.string());
	//Prepare the positiveSet.fa sequences and send them to PolyFudBpp for prediction
	for (auto & pas : pasVector) {
		auto txProp = refGen.getValueByKey(pas.seqId);
		size_t utrLength = polar::utility::getUtrLength(txProp);
		size_t txPos = polar::utility::mapGenomePosToTxPos(txProp, pas.pos);
		size_t txLen = polar::utility::getTxLength(txProp);
		size_t pasUtrPos = txPos - (txLen - utrLength);
		if (txPos == UINT_MAX) {
			std::cerr << "warning: couldn't map positions for " << pas.seqId << std::endl;
		}
		if (static_cast<int>(txPos) - (static_cast<int>(txLen) - static_cast<int>(utrLength)) < 0 ) {
			std::cerr << "uint overflow: " << pas.seqId << " | txPos: " << txPos
				<< " | txLen: " << txLen 
				<< " | utrLength: " << utrLength << std::endl;
		}
		if (argc > 1 &&std::string(argv[1]) == "verbose") {
			std::cerr << pas.seqId << " | " << "utrLength: " << utrLength  << " | ";
		}
		
		if (utrLength < maxUtrLength && txPos != UINT_MAX) {
			PolyFudBpp tempObj = PolyFudBpp(pas.seqId);
			resVector.push_back(std::make_pair(tempObj, pasUtrPos));
			evaluatedPas++;
			if (uniques.insert(tempObj.getTxId()).second) {
				std::vector<double> bppVector = tempObj.getMaxBppVector();
				bppOutStream << ">" << tempObj.getTxId() << "\n";
				for(double & d : bppVector) {
					bppOutStream << d << " ";
				}
				bppOutStream << "\n";
			}
		} else {
			if (utrLength >= maxUtrLength) {
				
				numTooLongSequences++;
			}
			++ignoredPas;
		}
		if (argc > 1 && std::string(argv[1]) == "verbose") {
			double process = static_cast<double>(loopCounter) / pasVector.size() * 100;
			std::cout << "Process: " <<  std::setprecision(2) << std::setw(6) << std::fixed << process << "%" << std::endl;
				//<< "\r" << std::flush;
		}
		++loopCounter;
	}
	bppOutStream.close();

	std::cout << std::endl << "ignoredPas: " << ignoredPas << std::endl;
	std::cout << "number of utr longer than " << maxUtrLength << ": " << numTooLongSequences << std::endl;
	std::cout << "evaluatedPas: " << evaluatedPas << std::endl;
	
	size_t total = resVector.size();
	double threshold = 0.0;
	double stepSize = 0.05;
	//sensitivity evaluation
	for (;threshold <= 2.0; threshold += stepSize) {
		size_t notFound = 0;
		std::ofstream sensitivityOut("sensitivityOut.csv");
		for (auto & resObject : resVector) {
			auto posObjVector = resObject.first.getResults(threshold);
			bool found = false;
			for (auto & posObj : posObjVector) {
				if (posObj.pos == resObject.second) {
					sensitivityOut << resObject.first.getTxId() << "|" 
						<< posObj.pos << ", "
						<< posObj.truthValue << std::endl;
					found = true;
					break;
				}
			}
			if (! found) {
				notFound++;
			}
		}	
		sensitivityOut.close();
		double sensitivity = (static_cast<double>(total) - notFound) / total;
		std::cerr << "Sensitivity at threshold: " << std::fixed << std::setprecision(2) << threshold << " | " << sensitivity << ", TP+FN=" << total << ", TP=" << total - notFound << std::endl;
	}
	
	//specificity calculations
	resVector.clear();	
	PolyFudBpp::utrBppMap.clear();
	PolyFudBpp::utrBppMap = setNewBppDictionary(bppDictTn);
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		std::cerr << "could not open index file for " << referenceGenome << std::endl;
	}
	std::ofstream bppOutputTn(bppOutTn.string());
	loopCounter = 0;
	for (auto & pas : pasVector) {
		auto txProp = refGen.getValueByKey(pas.seqId);
		size_t utrLength = polar::utility::getUtrLength(txProp);
		size_t txPos = polar::utility::mapGenomePosToTxPos(txProp, pas.pos);
		size_t txLen = polar::utility::getTxLength(txProp);
		size_t pasUtrPos = txPos - (txLen - utrLength);
		if (txPos == UINT_MAX) {
			std::cerr << "warning: couldn't map positions for " << pas.seqId << std::endl;
		}
		if (static_cast<int>(txPos) - (static_cast<int>(txLen) - static_cast<int>(utrLength)) < 0 ) {
			std::cerr << "uint overflow: " << pas.seqId << " | txPos: " << txPos
				<< " | txLen: " << txLen 
				<< " | utrLength: " << utrLength << std::endl;
		}	
		if (argc > 1 &&std::string(argv[1]) == "verbose") {
			std::cerr << pas.seqId + "_" + std::to_string(pasUtrPos) << " | " << "utrLength: " << utrLength << " | ";  
		}
		
		if (utrLength < maxUtrLength && txPos != UINT_MAX) {
			//obtain utr sequence
			std::string utr = polar::utility::getUtrSequence(txProp, faiIndex);
			//obtain motif
			std::string motif = std::string(utr.begin() + pasUtrPos, utr.begin() + pasUtrPos + 6);
			//obtain sequence after utr
			std::string utrOffset = polar::utility::getSeqAfterUtr(txProp, faiIndex, 200);
			//shuffle the sequence after the utr
			std::random_shuffle(utrOffset.begin(), utrOffset.end());
			//delete PAS from utr sequence
			utr.erase(utr.begin() + pasUtrPos, utr.begin() + pasUtrPos + 6);
			//shuffle the utr
			std::random_shuffle(utr.begin(), utr.end());
			//reinsert the PAS into the utr sequence
			utr.insert(pasUtrPos, motif);
			
			//std::cerr << pas.seqId + "_" + std::to_string(txLen-utrLength) << std::endl; 
			PolyFudBpp tempObj = PolyFudBpp(pas.seqId + "_" + std::to_string(txLen-utrLength), utr, utrOffset); //+ "_" + std::to_string(pasUtrPos)
			resVector.push_back(std::make_pair(tempObj, pasUtrPos));
			evaluatedPas++;
			std::vector<double> bppVector = tempObj.getMaxBppVector();
			bppOutputTn << ">" << tempObj.getTxId() << "\n";
			for(double & d : bppVector) {
				bppOutputTn << d << " ";
			}
			bppOutputTn << "\n";
		} else {
			if (utrLength >= maxUtrLength) {
				std::string utr = polar::utility::getUtrSequence(txProp, faiIndex);
				std::string motif = std::string(utr.begin() + pasUtrPos, utr.begin() + pasUtrPos + 6);
				countDiscardedPas[motif]++;
				
			}
			++ignoredPas;
		}
		if (argc > 1 && std::string(argv[1]) == "verbose") {
			double process = static_cast<double>(loopCounter) / pasVector.size() * 100;
			std::cout << "Process: " <<  std::setprecision(2) << std::setw(6) << std::fixed << process << "%" << std::endl;
				//<< "\r" << std::flush;
		}
		++loopCounter;
	}
	total = resVector.size();
	threshold = 0.00;
	//specificity evaluation
	for (;threshold <= 2.0; threshold += stepSize) {
		size_t notFound = 0;
		std::ofstream specificityOut("specificityOut.csv");
		for (auto & resObject : resVector) {
			auto posObjVector = resObject.first.getResults(threshold);
			bool found = false;
			for (auto & posObj : posObjVector) {
				if (posObj.pos == resObject.second) {
					specificityOut << resObject.first.getTxId() << "|" 
						<< posObj.pos << ", "
						<< posObj.truthValue << std::endl;
					found = true;
					break;
				}
			}
			if (! found) {
				notFound++;
			}
		}	
		specificityOut.close();
		double specificity = static_cast<double>(notFound) / total;
		std::cerr << "Specificity at threshold: " << threshold << " | " << std::fixed << std::setprecision(2) << specificity << ", TP+FP=" << total << ", TP=" << notFound << std::endl;
	}
	std::ofstream posSetBppFuzzy("positiveSetBppFuzzy.fa");
	for(auto & pas : pasVector) {
		posSetBppFuzzy << pas.seqId << ", " << pas.pos << std::endl;
	}
	std::cerr << std::endl;
	for (auto & discPas : countDiscardedPas) {
		std::cerr << discPas.first << ": " << discPas.second << std::endl;
	}
	return EXIT_SUCCESS;
}
