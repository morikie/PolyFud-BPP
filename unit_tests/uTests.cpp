#define BOOST_TEST_MODULE polar_test
#define BOOST_TEST_MAIN polar 

#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/test/included/unit_test.hpp>
#include "../src/polarUtility.hpp"
#include "../src/polyFud.hpp"
#include "../src/polyFudBpp.hpp"
#include "../src/refGeneParser.hpp"
#include "../src/seqStruct.hpp"
#include "../src/utr3Finder.hpp"

extern "C" {
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/utils.h>
}


namespace fs = boost::filesystem;
namespace spirit = boost::spirit;

/**
 * Floating point comparison. Returns true if two floats are equal up to the desired ULP (units in the last place). 
 * Copied from http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon;
 */
template<typename T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
	almost_equal(T x, T y, int ulp)
{
	return std::abs(x-y) < std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
		|| std::abs(x-y) < std::numeric_limits<T>::min();
}

/* PolyFud tests */

BOOST_AUTO_TEST_CASE ( polyFud ) {	
	std::string seq = "agggagagattaatttttcccaccccaataaaccccccccacagtataaccaagagattttatttgaaaaacc";
	SeqStruct  txTest1 {
		seq,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none,
		boost::none
	};
	size_t utr3MotifPos = 26;
	PolyFud utr3FinderFuz_test1 (txTest1, true);
	auto pairLeft1 = PolyFud::dseLocMap.find(std::string("aataaa"))->second.getLeftStraight();
	auto pairRight1 = PolyFud::dseLocMap.find(std::string("aataaa"))->second.getRightStraight();
	double slopeLeft1 = 1.0 / 15;
	double interceptLeft1 = - 2.0 / 3;
	double slopeRight1 = - 1.0 / 20;
	double interceptRight1 = 11.0 / 4;
	BOOST_CHECK(almost_equal(pairLeft1.first, slopeLeft1, 2));
	BOOST_CHECK(almost_equal(pairLeft1.second, interceptLeft1, 2));
	BOOST_CHECK(almost_equal(pairRight1.first, slopeRight1, 2));
	BOOST_CHECK(almost_equal(pairRight1.second, interceptRight1, 2));

	auto pairLeft2 = PolyFud::dseLocMap.find(std::string("attaaa"))->second.getLeftStraight();
	auto pairRight2 = PolyFud::dseLocMap.find(std::string("attaaa"))->second.getRightStraight();
	double slopeLeft2 = 1.0 / 15;
	double interceptLeft2 = - 2.0 / 3;
	double slopeRight2 = - 1.0 / 27;
	double interceptRight2 = 20.0 / 9;
	BOOST_CHECK(almost_equal(pairLeft2.first, slopeLeft2, 2));
	BOOST_CHECK(almost_equal(pairLeft2.second, interceptLeft2, 2));
	BOOST_CHECK(almost_equal(pairRight2.first, slopeRight2, 2));
	BOOST_CHECK(almost_equal(pairRight2.second, interceptRight2, 2));
	BOOST_CHECK_EQUAL(utr3MotifPos, utr3FinderFuz_test1.getPolyaMotifPos()[0].pos);
}

/* polar utility tests */

BOOST_AUTO_TEST_CASE ( polarUtility ) {
	fs::path referenceGenome = "reference_genome/hg19/reference_genome.fa";
	fs::path refGenomeIndex = "reference_genome/hg19/reference_genome.fa.fai";
	fs::path refGene = "ucsc_data/refGene.txt";
	
	RefGeneParser parser(refGene);
	
	seqan::FaiIndex faiIndex;
	if (! seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str())) {
		if (polar::utility::buildIndexFile(referenceGenome)) {
			seqan::open(faiIndex, referenceGenome.c_str(), refGenomeIndex.c_str());
		} else {
			std::cerr << "could not open index file for " << referenceGenome << std::endl;
		}
	}
	
	std::string testId = "NM_004081";
	std::string testUtr = "taaattccgttgttactcaagatgactgcttcaagggtaaaagagtgcatcgctttagaagaagtttggcagtatttaaatctgttggatcctctcagctatctagtttcatgggaagttgctggttttgaatattaagctaaaagttttccactattacagaaattctgaattttggtaaatcacactgaaactttctgtataacttgtattattagactctctagttttatcttaacactgaaactgttcttcattagatgtttatttagaacctggttctgtgtttaatatatagtttaaagtaacaaataatcgagactgaaagaatgttaagatttatctgcaaggatttttaaaaaattgaaacttgcattttaagtgtttaaaagcaaatactgactttcaaaaaagtttttaaaacctgatttgaaagctaacaattttgatagtctgaacacaagcatttcacttctccaagaagtacctgtgaacagtacaatatttcagtattgagctttgcatttatgatttatctagaaatttacctcaaaagcagaatttttaaaactgcatttttaatcagtggaactcaatgtatagttagctttattgaagtcttatccaaacccagtaaaacagattctaagcaaacagtccaatcagtgagtcataatgtttattcaaagtattttatcttttatctagaatccacatatgtatgtccaatttgattgggatagtagttaggataactaaaattctgggcctaattttttaaagaatccaagacaaactaaactttactgggtatataaccttctcaatgagttaccattcttttttataaaaaaaattgttccttgaaatgctaaacttaatggctgtatgtgaaatttgcaaaatactggtattaaagaacgctgcagcttttttatgtcactcaaaggttaatcggagtatctgaaaggaattgtttttataaaaacattgaagtattagttacttgctataaatagatttttatttttgttttttagcctgttatatttccttctgtaaaataaaatatgtccagaagaggcatgttgtttctagattaggtagtgtcctcattttatattgtgaccacacagctagagcaccagagcccttttgctatactcacagtcttgttttcccagcctcttttactagtctttcaggaggtttgctcttagaactggtgatgtaaagaatggaagtagctgtatgagcagttcaaaggccaagccgtggaatggtagcaatgggatataatacctttctaagggaaacatttgtatcagtatcatttgatctgccatggacatgtgtttaaagtggctttctggcccttctttcaatggcttcttccctaaaacgtggagactctaagttaatgtcgttactatgggccatattactaatgcccactggggtctatgatttctcaaaattttcattcggaatccgaaggatacagtctttaaactttagaattcccaagaaggctttattacacctcagaaattgaaagcaccatgactttgtccattaaaaaattatccatagtttttttagtgcttttaacattccgacatacatcattctgtgattaaatctccagatttctgtaaatgatacctacattctaaagagttaattctaattattccgatatgaccttaaggaaaagtaaaggaataaatttttgtctttgttgaagtatttaatagagtaaggtaaagaagatattaagtccctttcaaaatggaaaattaattctaaactgagaaaaatgttcctactacctattgctgatactgtctttgcataaatgaataaaaataaactttttttcttcaaatgtg";
	std::string testOffset = "tttttggctttccgatgtaataatgtaaaatggtggggagttgcgtgggaactgtgtaacaaggtttaaattcgtataacaagctttagattcttaaaatgcagaagtataaagttcagtatactaatctgtctgagttagcccataaaagcaaatgtaggtacaaagataagtttaagaggtgcatcaacagcagtgcag";
	auto txPro = parser.getValueByKey(testId);
	std::string extractedUtr = polar::utility::getUtrSequence(txPro, faiIndex);
	size_t testLengthFromGenome = polar::utility::getUtrLength(txPro);
	BOOST_CHECK_EQUAL(extractedUtr, testUtr);
	BOOST_CHECK_EQUAL(testLengthFromGenome, testUtr.size());

	size_t zero = 25345239;
	size_t two = 25345237;
	size_t shouldBeZero = polar::utility::mapGenomePosToTxPos(txPro, zero);
	size_t shouldBeTwo = polar::utility::mapGenomePosToTxPos(txPro, two);
	BOOST_CHECK_EQUAL(shouldBeZero, 0u); 
	BOOST_CHECK_EQUAL(shouldBeTwo, 2u);
	
	
	std::string testId2 = "NM_018836";
	std::string testUtr2 = "ctggccgaagtcttttttacctcctgggggcagggcagacgccgtgtgtctgtttcacggattccgttggtgaacctgtaaaaacaaaacaaacaaaacaaaacaaaaaagacaaaacctaaaactgagctatctaagggggagggtccccgcacctaccacttctgtttgccggtgggaaactcacagagcaggacgctctaggccaaatctatttttgtaaaaatgctcatgcctatgggtgactgccttctcccagagttttctttggagaacagaaagaagaaaggaaagaaaggaaccagaggcagagagacgaggatacccagcgaaagggacgggaggaagcatccgaaacctaggattcgtcctacgattctgaacctgtgccaataataccattatgtgccatgtactgacccgaaaggctcggccgcagagccggggcccagcgaatcacgcagagaaatcttacagaaaacaggggtgggaatctcttccgatagagtcgctatttctggttaatatacatatataaatatataaatacaaacacacacacacactttttttgtactgtagcaatttttgaagatcttaaatgttcctttttaaaaaaaagaattgtgttataggttacaaaatctgatttatttaacatgcttagtatgagcagaataaaccagtgttttctactttggcaactcacgtcacacacatattacacacatgtgcgcattacacacacacaatacacatacatgcatatagacgcatctattggaaatgcagttccacaggtgagcatgttctttctggtgacctggtattccatcaccattcaccccaggggacagcctcgaccgagacaaggaggcccttaaatgacagcctgcatttgctagacggttggtgagtggcatcaaatgtgtgacttactatcttgggccagaactaagaatgccaaggttttatatatgtgtgtgtatatatatatatatatatatatatatatatatatatatgtttgtgtgtgtatatatatatatatatatatatgtttgtgtgtgtatatatatgtttgtgtatatatatacacatatgcatacatatgatttttttttttcatttaagtgttggaagatgctacctaacagccacgttcacatttacgtagctggttgcttacaaacgggcctgagcccctggttgggtgggtggtggattcttggacgtgtgtgtcatacaagcatagactggattaaagaagttttccagttccaaaaattaaaggaatatatcctta";
	std::string testOffset2 = "tgatgtgtgtgtgtaatatcagggcagaacttagacatacgtgaagggccccggttggtttgaaaacgaaaaatagtcattctgtgtgcaaaccacaaggctgccccagtcaggcagcgccctgacctggcctgtgctgcattgccttcccttgcgcaggtgggcaggtgtggcccgctttttctagggcccaagggtgac";
	txPro = parser.getValueByKey(testId2);
	extractedUtr = polar::utility::getUtrSequence(txPro, faiIndex);
	testLengthFromGenome = polar::utility::getUtrLength(txPro);
	BOOST_CHECK_EQUAL(testLengthFromGenome, testUtr2.size());
	BOOST_CHECK_EQUAL(extractedUtr, testUtr2);

	zero = 4715104;
	two = 4715106;
	shouldBeZero = polar::utility::mapGenomePosToTxPos(txPro, zero);
	shouldBeTwo = polar::utility::mapGenomePosToTxPos(txPro, two);
	BOOST_CHECK_EQUAL(shouldBeZero, 0u);
	BOOST_CHECK_EQUAL(shouldBeTwo, 2u);
		
}

BOOST_AUTO_TEST_CASE( polyFudBpp ) {

	//Test if data retrieved from reference genome is the same. This transcript is classified as "-"-strand sequence. Returned sequence should be reverse complemented.
	std::string testId = "NM_004081";
	std::string testUtr = "taaattccgttgttactcaagatgactgcttcaagggtaaaagagtgcatcgctttagaagaagtttggcagtatttaaatctgttggatcctctcagctatctagtttcatgggaagttgctggttttgaatattaagctaaaagttttccactattacagaaattctgaattttggtaaatcacactgaaactttctgtataacttgtattattagactctctagttttatcttaacactgaaactgttcttcattagatgtttatttagaacctggttctgtgtttaatatatagtttaaagtaacaaataatcgagactgaaagaatgttaagatttatctgcaaggatttttaaaaaattgaaacttgcattttaagtgtttaaaagcaaatactgactttcaaaaaagtttttaaaacctgatttgaaagctaacaattttgatagtctgaacacaagcatttcacttctccaagaagtacctgtgaacagtacaatatttcagtattgagctttgcatttatgatttatctagaaatttacctcaaaagcagaatttttaaaactgcatttttaatcagtggaactcaatgtatagttagctttattgaagtcttatccaaacccagtaaaacagattctaagcaaacagtccaatcagtgagtcataatgtttattcaaagtattttatcttttatctagaatccacatatgtatgtccaatttgattgggatagtagttaggataactaaaattctgggcctaattttttaaagaatccaagacaaactaaactttactgggtatataaccttctcaatgagttaccattcttttttataaaaaaaattgttccttgaaatgctaaacttaatggctgtatgtgaaatttgcaaaatactggtattaaagaacgctgcagcttttttatgtcactcaaaggttaatcggagtatctgaaaggaattgtttttataaaaacattgaagtattagttacttgctataaatagatttttatttttgttttttagcctgttatatttccttctgtaaaataaaatatgtccagaagaggcatgttgtttctagattaggtagtgtcctcattttatattgtgaccacacagctagagcaccagagcccttttgctatactcacagtcttgttttcccagcctcttttactagtctttcaggaggtttgctcttagaactggtgatgtaaagaatggaagtagctgtatgagcagttcaaaggccaagccgtggaatggtagcaatgggatataatacctttctaagggaaacatttgtatcagtatcatttgatctgccatggacatgtgtttaaagtggctttctggcccttctttcaatggcttcttccctaaaacgtggagactctaagttaatgtcgttactatgggccatattactaatgcccactggggtctatgatttctcaaaattttcattcggaatccgaaggatacagtctttaaactttagaattcccaagaaggctttattacacctcagaaattgaaagcaccatgactttgtccattaaaaaattatccatagtttttttagtgcttttaacattccgacatacatcattctgtgattaaatctccagatttctgtaaatgatacctacattctaaagagttaattctaattattccgatatgaccttaaggaaaagtaaaggaataaatttttgtctttgttgaagtatttaatagagtaaggtaaagaagatattaagtccctttcaaaatggaaaattaattctaaactgagaaaaatgttcctactacctattgctgatactgtctttgcataaatgaataaaaataaactttttttcttcaaatgtg";
	std::string testOffset = "tttttggctttccgatgtaataatgtaaaatggtggggagttgcgtgggaactgtgtaacaaggtttaaattcgtataacaagctttagattcttaaaatgcagaagtataaagttcagtatactaatctgtctgagttagcccataaaagcaaatgtaggtacaaagataagtttaagaggtgcatcaacagcagtgcag";
	
	PolyFudBpp tempObj = PolyFudBpp(testId);
	std::string tempObjUtr = tempObj.getUtrSeq();
	std::string tempObjOffsetSeq = tempObj.getOffsetSeq();
	BOOST_CHECK_EQUAL(testUtr, tempObjUtr);
	BOOST_CHECK_EQUAL(tempObjOffsetSeq, testOffset);

	//Test if positions in the result object are correct. The positions are the real PAS positions, they must be found by the naive search.
	//Naturally, the truth value might classify them as non-functional but that is not for the test to check.
	//Sequence contains two known PASes.
	size_t pasPos1 = 1861;
	size_t pasPos2 = 1867;
	auto res1 = tempObj.getResults(0.0);
	bool found1=false, found2=false;
	for (auto & obj : res1) {
		if (obj.pos == pasPos1) found1 = true;
		if (obj.pos == pasPos2) found2 = true;
	}
	BOOST_CHECK(found1 && found2);

	//Same tests as above. Different sequence, this time from "+"-strand.
	std::string testId2 = "NM_018836";
	std::string testUtr2 = "ctggccgaagtcttttttacctcctgggggcagggcagacgccgtgtgtctgtttcacggattccgttggtgaacctgtaaaaacaaaacaaacaaaacaaaacaaaaaagacaaaacctaaaactgagctatctaagggggagggtccccgcacctaccacttctgtttgccggtgggaaactcacagagcaggacgctctaggccaaatctatttttgtaaaaatgctcatgcctatgggtgactgccttctcccagagttttctttggagaacagaaagaagaaaggaaagaaaggaaccagaggcagagagacgaggatacccagcgaaagggacgggaggaagcatccgaaacctaggattcgtcctacgattctgaacctgtgccaataataccattatgtgccatgtactgacccgaaaggctcggccgcagagccggggcccagcgaatcacgcagagaaatcttacagaaaacaggggtgggaatctcttccgatagagtcgctatttctggttaatatacatatataaatatataaatacaaacacacacacacactttttttgtactgtagcaatttttgaagatcttaaatgttcctttttaaaaaaaagaattgtgttataggttacaaaatctgatttatttaacatgcttagtatgagcagaataaaccagtgttttctactttggcaactcacgtcacacacatattacacacatgtgcgcattacacacacacaatacacatacatgcatatagacgcatctattggaaatgcagttccacaggtgagcatgttctttctggtgacctggtattccatcaccattcaccccaggggacagcctcgaccgagacaaggaggcccttaaatgacagcctgcatttgctagacggttggtgagtggcatcaaatgtgtgacttactatcttgggccagaactaagaatgccaaggttttatatatgtgtgtgtatatatatatatatatatatatatatatatatatatatgtttgtgtgtgtatatatatatatatatatatatgtttgtgtgtgtatatatatgtttgtgtatatatatacacatatgcatacatatgatttttttttttcatttaagtgttggaagatgctacctaacagccacgttcacatttacgtagctggttgcttacaaacgggcctgagcccctggttgggtgggtggtggattcttggacgtgtgtgtcatacaagcatagactggattaaagaagttttccagttccaaaaattaaaggaatatatcctta";
	std::string testOffset2 = "tgatgtgtgtgtgtaatatcagggcagaacttagacatacgtgaagggccccggttggtttgaaaacgaaaaatagtcattctgtgtgcaaaccacaaggctgccccagtcaggcagcgccctgacctggcctgtgctgcattgccttcccttgcgcaggtgggcaggtgtggcccgctttttctagggcccaagggtgac";

	PolyFudBpp tempObj2 = PolyFudBpp(testId2);
	std::string tempObj2Utr = tempObj2.getUtrSeq();
	std::string tempObj2OffsetSeq = tempObj2.getOffsetSeq();
	BOOST_CHECK_EQUAL(testUtr2, tempObj2Utr);
	BOOST_CHECK_EQUAL(tempObj2OffsetSeq, testOffset2);

	//Positional test as above. Sequence contains three known PASes.
	pasPos1 = 676;
	pasPos2 = 1260;
	size_t pasPos3 = 1286;
	auto res2 = tempObj2.getResults(0.0);
	found1=false, found2=false;
	bool found3=false;
	for (auto & obj : res2) {
		if (obj.pos == pasPos1) found1 = true;
		if (obj.pos == pasPos2) found2 = true;
		if (obj.pos == pasPos3) found3 = true;
	}
	BOOST_CHECK(found1 && found2 && found3);
}

