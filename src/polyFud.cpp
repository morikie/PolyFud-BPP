#include <algorithm>
#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <boost/foreach.hpp>
#include "polarUtility.hpp"
#include "seqStruct.hpp"
#include "polyFud.hpp"


/**
 * Different TV for the dse location map.
 */
PolyFud::pasToDseLocMap PolyFud::dseLocMap = { 
	{std::string("aataaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 55))},
	{std::string("attaaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(33, 60))},
	{std::string("tataaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("agtaaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("aagaaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("aatata"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("aataca"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("cataaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("gataaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("aatgaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("actaaa"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))},
	{std::string("aataga"), PolyFud::DseLocation(std::make_pair(10, 25), std::make_pair(35, 60))}
};


/**
 * Initialization of static member that holds value for the lower/upper bound of the DSE uracil content.
 * Taken from the original paper ("Prediction of non-canonical polyadenylation signals..."; doi: 10.1016/j.jbiosc.2009.01.001)
 */
PolyFud::pasToUcontentMap PolyFud::dseUracilMap = {
	{std::string("aataaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("attaaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("tataaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("agtaaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aagaaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aatata"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aataca"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("cataaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("gataaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aatgaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("actaaa"), PolyFud::UracilContent(0.33, 0.78, 1.0)},
	{std::string("aataga"), PolyFud::UracilContent(0.33, 0.78, 1.0)}
};


/**
 * Initialization of USE uracil parameters.
 */
PolyFud::pasToUcontentMap PolyFud::useUracilMap = {
	{std::string("aataaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("attaaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("tataaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("agtaaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("aagaaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("aatata"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("aataca"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("cataaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("gataaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("aatgaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("actaaa"), PolyFud::UracilContent(0.26, 0.78, 0.4)},
	{std::string("aataga"), PolyFud::UracilContent(0.26, 0.78, 0.4)}
};


/**
 * Thresholds for the listed PAS hexamers.
 */
std::unordered_map<PolyFud::motifSequence, double> PolyFud::thresholdMap = {
	{std::string("aataaa"), 0.90},
	{std::string("attaaa"), 0.90},
	{std::string("tataaa"), 0.70},
	{std::string("agtaaa"), 0.70},
	{std::string("aagaaa"), 0.70},
	{std::string("aatata"), 0.70},
	{std::string("aataca"), 0.70},
	{std::string("cataaa"), 0.70},
	{std::string("gataaa"), 0.70},
	{std::string("aatgaa"), 0.70},
	{std::string("actaaa"), 0.70},
	{std::string("aataga"), 0.70}
};


/**
 * Constructor.
 */
PolyFud::PolyFud(const SeqStruct & sSt, const bool & seaBack):
	Utr3Finder(sSt),
	searchBackward(seaBack)	
{
	if (this->searchBackward) {
		const std::string & seq = this->seqStruct.seq;
		std::transform(seq.rbegin(), seq.rend(), std::back_inserter(this->reversedSeq), polar::utility::complement);
	}
	this->findPolyaMotif();

}


/**
 * Destructor.
 */
PolyFud::~PolyFud() {}


/**
 * Predicts the Poly(A) motives in a sequence. Entry point after constructor call.
 */
void PolyFud::findPolyaMotif() {
	const std::string & seq = this->seqStruct.seq;
	const std::string & revSeq = this->reversedSeq;
	std::vector<size_t> candidatePositions;
	std::vector<size_t> revCandidatePositions;
	//searching forward strand
	BOOST_FOREACH (const PolyFud::pasToDseLocMap::value_type & v, PolyFud::dseLocMap) {
		auto posIt = seq.begin();
		while ((posIt = std::search(posIt, seq.end(), v.first.begin(), v.first.end())) != seq.end()) {
			candidatePositions.push_back(std::distance(seq.begin(), posIt));
			posIt++;
		}
	}
	
	//searching backward strand
	if (this->searchBackward) {
		BOOST_FOREACH (const PolyFud::pasToDseLocMap::value_type & v, PolyFud::dseLocMap) {
			auto posIt = revSeq.begin();
			while ((posIt = std::search(posIt, revSeq.end(), v.first.begin(), v.first.end())) != revSeq.end()) {
				revCandidatePositions.push_back(std::distance(revSeq.begin(), posIt));
				posIt++;
			}
		}
	}
	//verifying forward candidates
	BOOST_FOREACH(const size_t & candPos, candidatePositions) {
		std::string motif(seq.begin() + candPos, seq.begin() + candPos + 6);
		double combDseTValue;
		double useTvalue;
		double finalTvalue;
		
		combDseTValue = calcCombinedDseTvalue(candPos, seq);
		useTvalue = calcUseTvalue(candPos, seq);
		finalTvalue = combDseTValue + useTvalue;
		if (finalTvalue >= thresholdMap.find(motif)->second) {
			polyaPosVector.push_back(Utr3FinderResult{
				candPos, 
				finalTvalue,
				(this->searchBackward) ? "+" : "0"
				});
		}
	}
	//verifying backward candidates
	if (this->searchBackward) {
		BOOST_FOREACH(const size_t & candPos, revCandidatePositions) {
			std::string motif(revSeq.begin() + candPos, revSeq.begin() + candPos + 6);
			double combDseTValue;
			double useTvalue;
			double finalTvalue;
			
			combDseTValue = calcCombinedDseTvalue(candPos, revSeq);
			useTvalue = calcUseTvalue(candPos, revSeq);
			finalTvalue = combDseTValue + useTvalue;
			if (finalTvalue >= thresholdMap.find(motif)->second) {
				polyaPosVector.push_back(Utr3FinderResult{
					seq.size() - 1 - candPos,
					finalTvalue,
					"-"
					});
			}
		}
	}
}


/*
 * Scans downstream of a PAS candidate(i.e. its position)  for an uracil-rich region (the DSE).
 * Returns the combined truth value of uracil content and location for a certain window.
 * 
 */
double PolyFud::calcCombinedDseTvalue(const size_t & pos, const std::string & seq) {		
	std::string motif(seq.begin() + pos, seq.begin() + pos + 6);
	DseLocation & dseLoc = PolyFud::dseLocMap.find(motif)->second;
		
	auto start = seq.begin() + pos + dseLoc.getLeftRange().first + 6;
	std::list<std::string::value_type> slidingWindow(start, start + this->windowSize);
	
	double dseLocTvalue = 0.0;
	double dseUcontentTvalue = 0.0;
	double truthValue = 0.0;
	double maxTruthValue = 0.0;
	size_t uracilCounter = std::count(slidingWindow.begin(),slidingWindow.end(), 't');
	
	std::string::const_iterator end;
	
	if (std::distance(start, seq.end()) <static_cast<int>(this->windowSize + 1)) return 0.0;
	
	if (dseLoc.getRightRange().second - dseLoc.getLeftRange().first < static_cast<size_t>(std::distance(start, seq.end()))) {
		end = start + dseLoc.getRightRange().second - dseLoc.getLeftRange().first;
	} else {
		end = seq.end() - (this->windowSize);
	}

	for (auto posIt = start; posIt != end; posIt++) {
		double uContent = static_cast<double>(uracilCounter) / this->windowSize;
		size_t distanceToPas = std::distance(start - dseLoc.getLeftRange().first, posIt);
		dseLocTvalue = locationTvalue(motif, distanceToPas, this->dseLocMap);
		dseUcontentTvalue = nucleotideContentTvalue(motif, uContent, this->dseUracilMap);
		truthValue = std::min(dseLocTvalue, dseUcontentTvalue);
		if (truthValue > maxTruthValue) maxTruthValue = truthValue;
		
		if (slidingWindow.front() == 't' && uracilCounter > 0) uracilCounter--;
		slidingWindow.pop_front();
		slidingWindow.push_back(*(posIt + this->windowSize));
		if (slidingWindow.back() == 't') uracilCounter++;
		
	}
	return maxTruthValue;
}


/**
 * Scans upstream of a potential PAS for a uracil-rich region.
 * Returns the max truth value.
 */
double PolyFud::calcUseTvalue(const size_t & pos, const std::string & seq) {
	size_t searchRange = 20;
	std::string motif(seq.begin() + pos, seq.begin() + pos + 6);
	std::string::const_reverse_iterator start = seq.rend() - pos;
	std::list<std::string::value_type> slidingWindow(start, start + this->windowSize);
	
	double useUcontentTvalue = 0.0;
	double maxTvalue = 0.0;
	size_t uracilCounter = std::count(slidingWindow.begin(), slidingWindow.end(), 't');

	if (static_cast<size_t>(std::distance(start, seq.rend())) < this->windowSize + 1) return 0.0;
	
	std::string::const_reverse_iterator end;
	if (static_cast<size_t>(std::distance(start, seq.rend())) < searchRange) {
		end = seq.rend() - (this->windowSize);
	} else {
		end = start + searchRange - (this->windowSize);
	}

	for (auto posIt = start; posIt != end; posIt++) {
		double uContent = static_cast<double>(uracilCounter) / this->windowSize;
		useUcontentTvalue = nucleotideContentTvalue(motif, uContent, this->useUracilMap);
		if (useUcontentTvalue > maxTvalue) maxTvalue = useUcontentTvalue;
		
		if (slidingWindow.front() == 't' && uracilCounter > 0) uracilCounter--;
		slidingWindow.pop_front();
		slidingWindow.push_back(*(posIt + this->windowSize));
		if (slidingWindow.back() == 't') uracilCounter++;
	}
	return maxTvalue;
}


/**
 * Set the threshold map to determine a correct prediction.
 */
void PolyFud::setThresholdMap(std::unordered_map<std::string, double> & map) {
	PolyFud::thresholdMap = map;
}


/**
 * Checks if a variant hits the Poly(A) motif.
 */
bool PolyFud::isMutationInMotif() const {
	return false;
}


/**
 * Returns the sequence that was searched on.
 */
std::string PolyFud::getSequence() const {
	return this->seqStruct.seq;
}


/**
 * Returns the reverse complemented sequence that was searched on.
 */
std::string PolyFud::getRevComplementSeq() const {
	return this->reversedSeq;
}


/**
 * Returns all authentic PAS found in the sequence.
 * [Probably needs to be deleted. Redundant/not needed?]
 */
std::string PolyFud::getMotifSequence(const Utr3FinderResult & result) const {
	if (result.pos == Utr3Finder::noHitPos || result.pos > this->seqStruct.seq.size() - 6) return std::string();
	if (result.strand =="+") {
		auto motifStart = this->seqStruct.seq.begin() + result.pos;
		auto motifEnd = this->seqStruct.seq.begin() + result.pos + 6;
		return std::string(motifStart, motifEnd);	
	} else if (result.strand == "-") {
		auto motifStart = this->reversedSeq.begin() + result.pos;
		auto motifEnd = this->reversedSeq.begin() + result.pos + 6;
		return std::string(motifStart, motifEnd);	
	} else {
		return std::string();
	}
}


/**
 * Returns the positions of the found PAS.
 */
std::vector<Utr3Finder::Utr3FinderResult> PolyFud::getPolyaMotifPos() const {
	return this->polyaPosVector;
}


/**
 * Writes information about the variant (from VCF e.g.) and the predicted PASs.
 */
void PolyFud::writeInfo() const {
	std::stringstream ss;
	ss << "Poly(A) pos: ";
	if (this->polyaPosVector.empty()) {
		ss << "none found";
		return;
	}
	for (auto it = this->polyaPosVector.begin(); it != polyaPosVector.end(); it++) {
		ss << it->pos << "(" << it->strand << ") ";
	}
	if (this->seqStruct.genomicPos) ss << ", gPos: " << *(this->seqStruct.genomicPos);
	if (this->seqStruct.chrom) ss << ", chrom: " << *(this->seqStruct.chrom);
	
	std::cerr << ss.str() << std::endl;
}


/**
 * Returns the truth value for a given position of a potential DSE.
 */
double PolyFud::locationTvalue(const std::string & pas, const size_t & pos, const PolyFud::pasToDseLocMap & map) const {
	const DseLocation & dseLoc = map.find(pas)->second;
	PolyFud::DseLocation::range leftRange = dseLoc.getLeftRange();
	PolyFud::DseLocation::range rightRange = dseLoc.getRightRange();
	PolyFud::DseLocation::straight leftStraight = dseLoc.getLeftStraight();
	PolyFud::DseLocation::straight rightStraight = dseLoc.getRightStraight();

	if (pos >= leftRange.second && pos <= rightRange.first) {
		return 1.0;
	} else if (pos < leftRange.first || pos > rightRange.second) {
		return 0.0;
	} else if (pos > leftRange.first && pos < leftRange.second) {
		return leftStraight.first * static_cast<double>(pos) + leftStraight.second;
	} else {
		return rightStraight.first * static_cast<double>(pos) + rightStraight.second;
	}
}


/**
 * Returns the truth value for a given uracil content for a DSE.
 */
double PolyFud::nucleotideContentTvalue(const std::string & pas, const double & uContent, const PolyFud::pasToUcontentMap & map) const {	
	const UracilContent & nucleotideContent = map.find(pas)->second;
	PolyFud::UracilContent::straight interStraight = nucleotideContent.getStraight();
	double uB = nucleotideContent.getUpperBound();
	double lB = nucleotideContent.getLowerBound();
	double maxTv = nucleotideContent.getMaxTruthValue();
	
	if (uContent >= uB) {
		return maxTv;
	} else if (uContent <= lB) {
		return 0.0;
	} else {
		return interStraight.first * uContent + interStraight.second;
	}
}


/**
 * Constructor DseLocation.
 */
PolyFud::DseLocation::DseLocation(PolyFud::DseLocation::range p, PolyFud::DseLocation::range n):
	positiveIntermediate(p),
	negativeIntermediate(n)
{
	this->calcStraights();	
}


/**
 * Destructor
 */
PolyFud::DseLocation::~DseLocation() {}


/**
 * Calculates the slopes and intercepts of the two straights (the two sides of the trapezoid) for the values between zero and one.
 */
void PolyFud::DseLocation::calcStraights() {
	if (positiveIntermediate.second < positiveIntermediate.first) 
		throw std::invalid_argument("DseLocation: invalid range arguments");
	if (negativeIntermediate.second < negativeIntermediate.first) 
		throw std::invalid_argument("DseLocation: invalid range arguments");
	if (positiveIntermediate.first > negativeIntermediate.first || positiveIntermediate.first > negativeIntermediate.second) 
		throw std::invalid_argument("DseLocation: first range has larger values than second range");

	double slopePositiveIntermediate = 1.0 / 
		(static_cast<double>(positiveIntermediate.second) - static_cast<double>(positiveIntermediate.first));
	double slopeNegativeIntermediate = - 1.0 / 
		(static_cast<double>(negativeIntermediate.second) - static_cast<double>(negativeIntermediate.first));
	
	double interceptPositiveIntermediate = 1.0 - slopePositiveIntermediate * static_cast<double>(positiveIntermediate.second);
	double interceptNegativeIntermediate = 1.0 - slopeNegativeIntermediate * static_cast<double>(negativeIntermediate.first);

	this->positiveStraight = std::make_pair(slopePositiveIntermediate, interceptPositiveIntermediate);
	this->negativeStraight = std::make_pair(slopeNegativeIntermediate, interceptNegativeIntermediate);
}


/**
 * Getter.
 */
PolyFud::DseLocation::range PolyFud::DseLocation::getLeftRange() const {
	return this->positiveIntermediate;
}


/**
 * Getter.
 */
PolyFud::DseLocation::range PolyFud::DseLocation::getRightRange() const {
	return this->negativeIntermediate;
}


/**
 * Getter.
 */
PolyFud::DseLocation::straight PolyFud::DseLocation::getLeftStraight() const {
	return this->positiveStraight;
}


/**
 * Getter.
 */
PolyFud::DseLocation::straight PolyFud::DseLocation::getRightStraight() const  {
	return this->negativeStraight;
}


/**
 * Constructor UracilContent.
 */
PolyFud::UracilContent::UracilContent(double lB, double uB, double mTv):
	lowerBound(lB),
	upperBound(uB),
	maxTruthValue(mTv) 
{
	this->calcStraight();	
}


/**
 * Destructor
 */
PolyFud::UracilContent::~UracilContent() {}


/**
 * Calculates the slope and interecept of the straight for the values between zero and one.
 */
void PolyFud::UracilContent::calcStraight() {
	if (upperBound < lowerBound) 
		throw std::invalid_argument("UracilContent: lower bound bigger than upper bound");
	if (this->maxTruthValue < 0.0) 
		throw std::invalid_argument("UracilContent: maximum truth value is lower than zero");
	
	double slopeIntermediate = this->maxTruthValue / (static_cast<double>(upperBound) - static_cast<double>(lowerBound));
	double interceptIntermediate = this->maxTruthValue - slopeIntermediate * static_cast<double>(upperBound);

	this->intermediate = std::make_pair(slopeIntermediate, interceptIntermediate);
}	


/**
 * Getter.
 */
double PolyFud::UracilContent::getLowerBound() const {
	return this->lowerBound;
}


/**
 * Getter.
 */
double PolyFud::UracilContent::getUpperBound() const {
	return this->upperBound;
}


/**
 * Getter.
 */
double PolyFud::UracilContent::getMaxTruthValue() const {
	return this->maxTruthValue;
}


/**
 * Getter.
 */
PolyFud::UracilContent::straight PolyFud::UracilContent::getStraight() const {
	return this->intermediate;
}

