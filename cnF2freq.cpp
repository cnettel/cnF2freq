// cnF2freq, (c) Carl Nettelblad, Department of Information Technology, Uppsala University
// 2008-2016
//
// PlantImpute 1.5, with support for forward-backward and old "true" tree-style cnF2freq
// algorithm. This reduces mamoery requirements and speeds up imputation. A lot. Still
// som e kinks to work out. Please contact the author regarding proper references.
//
// carl.nettelblad@it.uu.se
//
// This code is allowed to be freely used for any commercial or research purpose. If the code is
// used or integrated into another project largely unchanged, attribution to the original author
// should be included. No warranties are given.


#define NDEBUG
// These defines fixed an error in one particular site installation of the Portland compiler.
#define _STLP_EXPOSE_GLOBALS_IMPLEMENTATION 1
#define _REENTRANT 1
#define _SECURE_SCL 0
// For MSCVC
// Note: MSVC OpenMP support is insufficient in current release. Disable OpenMP for compilation
// in MSVC.
//
// Recent releases of g++ on MacOS X and Linux, as well as the Intel C++ compiler on
// Linux (x86 and x64) and Windows have been tested. The Portland compiler collection works only
// with some optimization settings, some race conditions in STL, although some
// workarounds are used.
#define _CRT_SECURE_NO_WARNINGS

#include <vector>
#include <string.h>
#include <stdio.h>
#include <omp.h>


const int ANALYZE_FLAG_FORWARD = 0;
const int ANALYZE_FLAG_BACKWARD = 16;
const int ANALYZE_FLAG_STORE = 32;

float templgeno[8] = { -1, -0.5,
	0,  0.5,
	0,  0.5,
	1, -0.5 };



// _MSC_VER is here to be interpreted as any compiler providing TR1 C++ headers
//#ifdef _MSC_VER
//#include <array>
//#else
// Boost also provides an array implementation, that is largely compatible
//#include <boost/array.hpp>
//#endif
#include <array>

#include <boost/math/distributions/binomial.hpp>
#include <boost/static_assert.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
/*#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>*/
#include <string>

#include <boost/container/flat_map.hpp>
#include <boost/program_options.hpp>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/ast/variant.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/include/std_tuple.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <iostream>
#include <fstream>

#include <vector>

using namespace boost::spirit;
namespace po = boost::program_options;

#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <set>
#include <algorithm>
#include <math.h>
#include <bitset>
#include <map>
#include <float.h>
#include <numeric> // use these libraries
#include <type_traits>

// Will be in C++17...
template <class T>
constexpr std::add_const_t<T>& as_const(T& t) noexcept
{
  return t;
}

using namespace std; // use functions that are part of the standard library
#ifdef _MSC_VER
//using namespace tr1;
#else
using namespace boost;
#endif

using namespace boost;
//using namespace boost::mpi;
using namespace boost::random;
#define flat_map boost::container::flat_map


#define none cnF2freqNONE


#ifndef _MSC_VER
#define _isnan isnan
#define _finite isfinite
#endif


boost::random::mt19937 rng;
int myrand(int max)
{
	uniform_int<> orig(0, max - 1);

	return orig(rng);
}



struct negtrue
{
	const int value;
	negtrue(bool neg, int value) : value(neg ? -1 : (value >= 0 ? value : 0))
	{
	}

	operator bool() const
	{
		return value < 0;
	}
};

struct MarkerValueType {};

constexpr MarkerValueType MarkerValue;

struct MarkerVal
{
private:
	int val;
	
public:
    constexpr explicit MarkerVal(int val) : val(val) {}
	constexpr MarkerVal() : val(0) {}
	constexpr int value() const
	{
		return val;
	}

        constexpr bool operator == (const MarkerVal& rhs) const
	{
		return val == rhs.val;
	}

        constexpr bool operator != (const MarkerVal& rhs) const
	{
		return val != rhs.val;
	}

	constexpr bool operator < (const MarkerVal& rhs) const
	{
		return val < rhs.val;
	}
};

constexpr const MarkerVal operator* (const int& val, const MarkerValueType& rhs)
{
	return MarkerVal(val);
}

typedef pair<MarkerVal, MarkerVal> MarkerValPair;


const MarkerVal UnknownMarkerVal = (MarkerVal)0;
const MarkerVal sexmarkerval = 9 * MarkerValue;

const float maxdiff = 0.000005;

#include "settings.h"

bool early = false;

vector<double> markerposes;
vector<double> actrec[2];
//int selfgen = 0;
double genrec[3];
vector<unsigned int> chromstarts;

vector<int> markertranslation;
typedef vector<std::array<double, 5 > > MWTYPE;

template<class T> class vectorplus;

template<> class vectorplus<MWTYPE >
{
public:
	MWTYPE operator () (const MWTYPE& l, const MWTYPE& r)
	{
		MWTYPE result;
		result.resize(l.size());

		for (int i = 0; i < l.size(); i++)
		{
			for (int j = 0; j < 5; j++)
			{
				result[i][j] = l[i][j] + r[i][j];
			}
		}

		return result;
	}
};

template<class T> class vectorplus<vector<T> >
{
public:
	vector<T> operator () (const vector<T>& l, const vector<T>& r)
	{
		vector<T> result;
		result.resize(l.size());

		for (int i = 0; i < l.size(); i++)
		{
			result[i] = l[i] + r[i];
		}

		return result;
	}
};

/*namespace boost { namespace mpi {

template<class T>
struct is_commutative<vectorplus<vector<T> >, vector<T> >
: mpl::true_ { };

} }*/

MWTYPE markerweight;
MWTYPE markerweight2;

double discstep = 0.1;
double baserec[2];
int sexc = 1;


// Is a specific marker value a admissible as a match to marker b
// The logic we use currently represents "0" as unknown, anything else as a known
// marker value.
// NOTE, VALUE CHANGED FOR ZEROS!
template<bool zeropropagate> bool markermiss(MarkerVal& a, const MarkerVal b)
{
	// This is the real logic; we do not bind anything at all when zeropropagate is true
	if (zeropropagate) return false;

	if (a == UnknownMarkerVal)
	{
		if (!zeropropagate) a = b;
		return false;
	}
	if (b == UnknownMarkerVal && a != sexmarkerval) return false;

	return a != b;
}


// Move a tree-based binary flag up a generation. The structure of bit flags might look like
// (1)(3)(3), where each group of (3) is in itself (1)(2), where (2) is of course (1)(1).
int upflagit(int flag, int parnum, int genwidth)
{

	if (flag < 0) return flag;
	flag >>= parnum * (genwidth - 1);
	flag &= ((1 << (genwidth - 1)) - 1);

	return flag;
}

struct individ;

int generation = 1;
int shiftflagmode;
int lockpos[NUMSHIFTS];
int quickmark[NUMSHIFTS];
int quickgen[NUMSHIFTS];

template<class T2> class PerStateArray
{
public:
	typedef std::array<T2, NUMTYPES> T;
};

template<class T2> class StateToStateMatrix
{
public:
	typedef typename PerStateArray<
		typename PerStateArray<T2>::T >
		::T T;
};

#if !DOFB
// The quick prefixes are caches that retain the last invocation.
// Only used fully when we stop at marker positions exactly, i.e. not a general grid search.
double quickfactor[NUMSHIFTS];
PerStateArray<int>::T quickendmarker[NUMSHIFTS];
PerStateArray<double>::T quickendfactor[NUMSHIFTS];
StateToStateMatrix<double>::T quickendprobs[NUMSHIFTS];
PerStateArray<double>::T quickmem[NUMSHIFTS];
#endif

// A hashed store of inheritance pathway branches that are known to be impossible.
// Since we can track the two branches that make up the state in the F_2 individual independently,
// this optimization can reduce part of the cost by sqrt(number of states).
typedef std::array<std::array<std::array<std::array<std::array<std::array<std::array<int, 4>, HALFNUMSHIFTS>, HALFNUMPATHS + 1>, HALFNUMTYPES>, 2>, 2>, 2> IAT;

EXTERNFORGCC IAT impossible;

// A memory structure storing haplo information for later update.
// By keeping essentially thread-independent copies, no critical sections have to
// be acquired during the updates.
EXTERNFORGCC std::array<std::array<float, 2>, INDCOUNT> haplos;
EXTERNFORGCC std::array<std::array<flat_map<MarkerVal, float>, 2>, INDCOUNT> infprobs;

#if !DOFB
// done, factors and cacheprobs all keep track of the same data
// done indicates that a specific index (in the binary tree of blocks of multi-step transitions) is done
// with a "generation id" that's semi-unique, meaning no active clearing of the data structure is performed
EXTERNFORGCC vector<int> done[NUMSHIFTS];
// factors contain the mantissas of the extended floating-point representation
EXTERNFORGCC vector<PerStateArray<double>::T > factors[NUMSHIFTS];
// cacheprobs contain actual transitions from every possible state to every possible other state
EXTERNFORGCC vector<StateToStateMatrix<double>::T > cacheprobs[NUMSHIFTS];
#else
EXTERNFORGCC vector<std::array<PerStateArray<double>::T, 2> > fwbw[NUMSHIFTS];
EXTERNFORGCC vector<std::array<double, 2> > fwbwfactors[NUMSHIFTS];
int fwbwdone[NUMSHIFTS];
#endif

EXTERNFORGCC vector<individ*> reltree;
EXTERNFORGCC flat_map<individ*, int> relmap;
EXTERNFORGCC flat_map<individ*, int> relmapshift;

//#pragma omp threadprivate(realdone, realfactors, realcacheprobs)
#pragma omp threadprivate(generation, shiftflagmode, impossible, haplos, lockpos, reltree, relmap, relmapshift, infprobs)
#if !DOFB
#pragma omp threadprivate(quickmark, quickgen, quickmem, quickfactor, quickendfactor, quickendprobs, done, factors, cacheprobs)
#else
#pragma omp threadprivate(fwbw,fwbwfactors,fwbwdone)
#endif

#ifdef DOEXTERNFORGCC
IAT impossible;
std::array<std::array<float, 2>, INDCOUNT> haplos;
vector<PerStateArray<double>::T > factors[NUMSHIFTS];
vector<individ*> reltree;
flat_map<individ*, int> relmap; //containing flag2 indices
flat_map<individ*, int> relmapshift; //containing flag2 indices
std::array<std::array<flat_map<MarkerVal, float>, 2>, INDCOUNT> infprobs;
#if !DOFB
vector<int> done[NUMSHIFTS];
vector<StateToStateMatrix<double>::T > cacheprobs[NUMSHIFTS];
#else
vector<std::array<PerStateArray<double>::T, 2> > fwbw[NUMSHIFTS];
vector<std::array<double, 2> > fwbwfactors[NUMSHIFTS];
#endif
#endif




// We put all thread-local structures in our own separate struct. This is because many compilers implementing OpenMP
// use a relatively expensive call for determining thread-local global data, which is not cached over function calls.
// The simple pointer-arithmetic for a lookup within a struct is cheaper.
struct threadblock
{
	int* const generation;
	int* const shiftflagmode;
#if !DOFB
	int* const quickmark;
	int* const quickgen;
#endif
	int* const lockpos;
#if !DOFB
	double* const quickfactor;
	PerStateArray<double>::T* const quickmem;
	PerStateArray<int>::T* const quickendmarker;
#endif
	IAT* const impossible;
	std::array<std::array<float, 2>, INDCOUNT>* const haplos;
	std::array<std::array<flat_map<MarkerVal, float>, 2>, INDCOUNT>* infprobs;
#if !DOFB
	vector<int>* const done;
	vector<PerStateArray<double>::T >* const factors;
	vector<StateToStateMatrix<double>::T >* const cacheprobs;
	PerStateArray<double>::T* const quickendfactor;
	StateToStateMatrix<double>::T* const quickendprobs;
#else
	vector<std::array<PerStateArray<double>::T, 2> >* fwbw;
	vector<std::array<double, 2> >* fwbwfactors;
	int* fwbwdone;
#endif	

	threadblock() : generation(&::generation), shiftflagmode(&::shiftflagmode), impossible(&::impossible),
		haplos(&::haplos), lockpos(::lockpos), infprobs(&::infprobs),
#if !DOFB
		done(::done), factors(::factors), cacheprobs(::cacheprobs),
		quickmark(::quickmark), quickgen(::quickgen), quickmem(::quickmem),
		quickfactor(::quickfactor), quickendfactor(::quickendfactor), quickendprobs(::quickendprobs),
			quickendmarker(::quickendmarker),
#else
		fwbw(::fwbw), fwbwfactors(::fwbwfactors), fwbwdone(::fwbwdone)
#endif
	{
	};
};


// Turners are used as mix-ins when probabilities are filtered at the "fixed" marker
// (the one we are actually probing).
// This functionality is used to invert the state in different manner, to test for the
// possibility that the complete haplotype assignment from an arbitrary point until the end
// has been mixed up. This can easily happen as the assignment of haplotype numbers is
// basically arbitrary, so the only thing keeping it in place is the linkage to the original
// defining locus.
class noneturner
{
public:
	void operator() (const PerStateArray<double>::T& probs) const
	{
	};

	bool canquickend() const
	{
		return true;
	}
} none;

class aroundturner
{
	int turn;
	int flagmodeshift;

public:
	aroundturner(int turn)
	{
		if (NUMGEN == 3)
		{
			this->turn = turn & 54;
			flagmodeshift = (turn >> TYPEBITS) | ((turn & 1) ? 2 : 0) |
				((turn & 8) ? 4 : 0);
			// Mis-guided (?) attempt to correct for relskews
			if (RELSKEWSTATES && false)
			  {
			    this->turn |= (flagmodeshift & 1) << BITS_W_SELF;
			  }
		}
		else
		{
			this->turn = turn & 3;
			flagmodeshift = (turn >> TYPEBITS);
		}
	}

	bool canquickend() const
	{
		return false;
	}


	void operator() (PerStateArray<double>::T& probs) const
	{
		PerStateArray<double>::T probs2;

		for (int i = 0; i < NUMTYPES; i++)
		{
			int newval = i ^ turn;
			if (flagmodeshift)
			{
				// Emulating actually switching parental genotypes around
				//		      newval = ((newval & 1) << 1) | ((newval & 2) >> 1);
			}
			probs2[newval] = probs[i];
		}

		for (int i = 0; i < NUMTYPES; i++)
		{
			probs[i] = probs2[i];
		}

		// Not *tb out of laziness, the number of calls is limited
		shiftflagmode ^= flagmodeshift;
		// Hacky, tacky
	};
};

// A struct containing some auxiliary arguments to the trackpossible family of functions.
// These parameters are not known at compile-time, hence not template arguments, but they do not
// change over the series of recursive calls.
const struct trackpossibleparams
{
	const float updateval;
	int* const gstr;

	trackpossibleparams() : updateval(0.0f), gstr(0)
	{
	}

	trackpossibleparams(float updateval, int* gstr, MarkerVal markerval = UnknownMarkerVal) : updateval(updateval), gstr(gstr)
	{
	}
} tpdefault;



struct classicstop
{
	int lockpos;
	int genotype;

	classicstop(int lockpos, int genotype) : lockpos(lockpos), genotype(genotype)
	{}

	explicit operator int() const
	{
		return lockpos;
	}


	const int getgenotype(int startmark) const
	{
		return genotype;
	}



	negtrue okstep(const int startmark, const int endmark) const
	{
		bool found = (lockpos <= -1000 - startmark && lockpos > -1000 - endmark) ||
			(lockpos >= markerposes[startmark] && lockpos <= markerposes[endmark]);

		return negtrue(!found, genotype);
	}

	bool fixtofind(int& genotype, double& startpos, double& endpos, int j) const
	{
		bool tofind = (lockpos >= startpos && lockpos <= endpos);
		if (tofind)
		{
			endpos = lockpos;
		}


		if (lockpos == -1000 - j + 1)
		{
			tofind = true;
			endpos = startpos;
		}

		if (tofind)
		{
			genotype = this->genotype;
		}

		return tofind;
	}

};


struct smnonecommon
{
	const int getgenotype(int startmark) const
	{
		return -1;
	}
};

struct nonestop : smnonecommon
{
	operator int() const
	{
		assert(false);

		return -1;
	}

	negtrue okstep(const int startmark, const int endmark) const
	{
		return negtrue(true, 0);
	}

	bool fixtofind(int& genotype, double& startpos, double& endpos, int j) const
	{
		return false;
	}
} NONESTOP;

struct tssmcommon
{
	int lockpos;
	int genotype[2];


	operator int() const
	{
		return lockpos;
	}

	negtrue okstep(const int startmark, const int endmark) const
	{
		// twicestop is only used for markers
		bool found = (lockpos <= -1000 - startmark + 1 && lockpos > -1000 - endmark);

		if (found)
		{
			int index = -lockpos - startmark - 1000;

			if (index == 0 || index == 1) return negtrue(!found, genotype[index]);
		}

		return negtrue(!found, 0);
	}

	bool fixtofind(int& genotype, double& startpos, double& endpos, int j) const
	{
		bool two = lockpos == -1000 - j + 2;
		if (lockpos == -1000 - j + 1 || two)
		{
			endpos = startpos;
			genotype = this->genotype[two];
			return true;
		}

		return false;
	}

protected: tssmcommon(int lockpos) : lockpos(lockpos)
{
}
};

struct twicestop : tssmcommon
{
	twicestop(int lockpos, int* genotype) : tssmcommon(lockpos)
	{
		this->genotype[0] = genotype[0];
		this->genotype[1] = genotype[1];
	}

	twicestop(int lockpos, int g1, int g2) : tssmcommon(lockpos)
	{
		this->genotype[0] = g1;
		this->genotype[1] = g2;
	}

	const int getgenotype(int startmark) const
	{
		return genotype[-1000 - startmark != lockpos];
	}

};

struct stopmodpair : tssmcommon, smnonecommon
{
	typedef std::array<std::array<float, 2>, 2> miniactrecT;
	miniactrecT actrec;

	stopmodpair(int lockpos, miniactrecT actrec) : tssmcommon(lockpos), actrec(actrec)
	{
		genotype[0] = -1;
		genotype[1] = -1;
	}
};

template<class G> const bool canquickend(int startmark, const G &stopdata)
{
	return false;
}

/*template<> const bool canquickend<smnonecommon>(int startmark, const smnonecommon &stopdata)
{
return false;
}*/

template<> const bool canquickend<classicstop>(int startmark, const classicstop &stopdata)
{
	return stopdata.lockpos <= -1000 && stopdata.genotype >= 0;
}

template<> const bool canquickend<twicestop>(int startmark, const twicestop &stopdata)
{
	return stopdata.lockpos == -1000 - startmark + 1;
}


// Should only be called for markers where locking is really done



template<class G> bool firstlockmatch(int lockpos, const G& stopdata)
{
	return stopdata.lockpos == lockpos;
}


template<> bool firstlockmatch<nonestop>(const int lockpos, const nonestop& stopdata)
{
	return false;
}

#ifdef PERMARKERACTREC
template<class G> double getactrec(const G& stopdata, double startpos, double endpos, int k, int j, int gen)
{
	return actrec[k][j];
}
#else
template<class G> double getactrec(const G& stopdata, double startpos, double endpos, int k, int j, int gen)
{
	return genrec[gen];
}
#endif

template<> double getactrec<stopmodpair>(const stopmodpair& stopdata, double startpos, double endpos, int k, int j, int gen)
{
	int index = j - 1 + stopdata.lockpos + 1000;

	if (index == 0 || index == 1) return stopdata.actrec[k][index];

	return actrec[k][j];
}

const int HAPLOS = 1;
const int GENOS = 2;

struct clause
{
	long long weight;
	vector<int> cinds;
	//vector<individ*> individuals;

	string toString() {
		string s = boost::lexical_cast<std::string>(weight);
		//std::stringstream ints;
		//std::copy(cInds.begin(), cInds.end(), std::ostream_iterator<int>(ints, " "));
		for (int i = 0; i < cinds.size(); i++) {
			if (cinds[i]) {
				s = s + " " + boost::lexical_cast<std::string>(cinds[i]);
			}
		}
		//s = s + " " + ints.str();
		return s;// note each line ends in a space
	}

	string clausetostring() {
		string s = "";
		for (int i = 0; i < cinds.size(); i++) {
			if (cinds[i]) {
				s = s + " " + boost::lexical_cast<std::string>(cinds[i]);
			}
		}
		return s;
	}

	string weighttostring() {
		return boost::lexical_cast<std::string>(weight);
	}

};

// A structure containing most of the information on an individual
struct individ
{
	// The individual #.
	int n;
	string name;
	// Generation number. Convenient, while not strictly needed.
	int gen;
	// Number of children.
	int children;
	// Am I explicitly created to not contain genotypes?
	bool empty;
	// Do I lack (relevant) parents?
	bool founder;
	// Parents.
	individ* pars[2];
	//
	vector<individ*> kids;
	// Sex.
	bool sex;
	// Line or strain of origin, should only exist in founders.
	int strain;
	// Marker data as a list of pairs. No specific ordering assumed.
	vector<MarkerValPair > markerdata;
	vector<pair<double, double> > markersure;

	vector<MarkerValPair > priormarkerdata;
	vector<pair<double, double> > priormarkersure;

	// Temporary storage of all possible marker values, used in fixparents.
	vector<flat_map<MarkerVal, pair<int, double> > > markervals;
	// The haplotype weight, or skewness. Introducing an actual ordering of the value in markerdata.
	vector<float> haploweight;
	// Relative skewness, i.e. shifts between adjacent markers.
	vector<float> relhaplo;
	// The cost-benefit value of inverting the haplotype assignment from an arbitrary marker point on.
	vector<double> negshift;
	vector<int> lastinved;
	vector<unsigned int> lockstart;

	vector<std::array<flat_map<MarkerVal, double>, 2> > infprobs;

	vector<int> genotypegrid;

	// Accumulators for haplotype skewness updates.
	vector<double> haplobase;
	vector<double> haplocount;

	individ()
	{
		pars[0] = 0;
		pars[1] = 0;
		strain = 0;
		sex = false;
		empty = false;
		// False is the safe option, other code will be take shortcuts if founder is true.
		founder = false;
	}

	bool arerelated(individ* b, vector<individ*> stack = vector<individ*>(), int gens = 0)
	{
		if (gens > 2) return false;
		if (!this) abort();

		if (b == this) return true;
		if (find(stack.begin(), stack.end(), this) != stack.end())
		{
			return false;
		}

		stack.push_back(this);
		if (stack.size() == 1)
		{
			if (b && b->arerelated(this, stack, gens + 1)) return true;
		}

		if (pars[0] && pars[0]->arerelated(b, stack, gens + 1)) return true;
		if (pars[1] && pars[1]->arerelated(b, stack, gens + 1)) return true;

		for (int i = 0; i < kids.size(); i++)
		{
			for (int j = 0; j < b->kids.size(); j++)
			{
				if (b->kids[j] == kids[i]) return true;
			}
		}
		//stack.pop_back();

		return false;
	}


	// A wrapper for calling trackpossible. Additional logic introduced to do lookup in the "impossible" tables of branches
	// not allowed at the current point. This will generally mean most branches and reduces the number of states visited considerably
	// in typical cases.
	//
	// The double overload means that the class can be treated as a conventional function call, returning a double, when the pre-lookup
	// is not needed. In the "real" trackpossible method, one instance is called in that way, while the other one is pre-looked up.
	template<int update, bool zeropropagate> struct recursetrackpossible
	{
		individ* mother;
		int upflagr;
		int upflag2r;
		int upshiftr;
		double secondval;
		const trackpossibleparams& extparams;
		const int genwidth;
		const MarkerVal markerval;
		const int marker;
		const threadblock& tb;
		int firstpar;
		int* impossibleref;
		int impossibleval;
		bool prelok;

		recursetrackpossible(individ* mother, const threadblock& tb, MarkerVal markerval, double secondval, int marker, int upflag, int upflag2, int upshift, int genwidth, int f2n, int firstpar, int numrealpar, const trackpossibleparams& extparams) :
			genwidth(genwidth), extparams(extparams), markerval(markerval), marker(marker), tb(tb), firstpar(firstpar), mother(mother), secondval(secondval)

		{
			upflagr = upflagit(upflag, firstpar, genwidth);
			upflag2r = upflagit(upflag2, firstpar, genwidth >> (NUMGEN - NUMFLAG2GEN));
			upshiftr = upflagit(upshift, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));

			prelok = true;
			if (!zeropropagate && genwidth == (1 << (NUMGEN - 1)))
			{
				if (DOIMPOSSIBLE)
				{
					impossibleref = &(*tb.impossible)[*(tb.shiftflagmode) & 1][firstpar][f2n][upflagr][upflag2r + 1][upshiftr][marker & 3];
					impossibleval = (*tb.generation) * markerposes.size() + marker;
				}
				else
				{
					impossibleref = 0;
				}

				if (DOIMPOSSIBLE && impossibleref && *impossibleref == impossibleval)
				{
					prelok = false;
				}
			}
		}

		operator double()
		{
			if (!prelok)
			{
				return 0;
			}

			if (!mother || !mother->pars[firstpar])
			{
				return 1 + secondval;
			}

			double baseval =
				mother->pars[firstpar]->trackpossible<update, zeropropagate>(tb, markerval, secondval, marker,
					upflagr,
					upflag2r,
					upshiftr, extparams, genwidth >> 1);

			if (DOIMPOSSIBLE && impossibleref && !zeropropagate && !update && genwidth == (1 << (NUMGEN - 1)) && !baseval)
			{
				*impossibleref = impossibleval;
			}

			return baseval;
		}
	};

	// zeropropagate also implies that there is a gstr value, we propagate zeros to find any possible source strain
	// The main logic of tracking a specific inheritance pathway, computing haplotype weights, and overall feasibility.
	// update: Should haplotype weights be updated?
	// zeropropagate: Should zero marker values be kept, or changed to the actual values they are matched against.
	// threadblock: Reference to all thread-local data.
	// inmarkerval: The marker val we are matching against.
	// marker: The marker number.
	// flag: The actual genotype flag. Note that this is zero-extended to be one bit more than the actual enumeration of flags
	// (lowest bit always 0).
	// flag99: What's mostly called flag2. A specific shift state of marker numbers, or -1 to indicate that all shifts are allowed.
	// localshift: Based on shiftflagmode, the overall mapping *for the complete sequence* of strand numbers to actual parentage for all
	// individuals in the analysis.
	// extparams: external parameters.
	// genwidth: The width of the generation flags.
	template<int update, bool zeropropagate> double trackpossible(const threadblock& tb, MarkerVal inmarkerval, double secondval, const unsigned int marker,
		const unsigned int flag, const int flag99, int localshift = 0, const trackpossibleparams& extparams = tpdefault,
		const int genwidth = 1 << (NUMGEN - 1)) /*const*/
	{
		// This used to be a nice silent null check. Compilers don't like that, so we abort in order to try to detect those cases.
		if (this == NULL) abort();

		// TYPEBITS are the ordinary TYPEBITS. Anything set beyond those indicates selfing. g is multiplied by 2 to become flag, hence TYPEBITS + 1
		const bool rootgen = (genwidth == (1 << (NUMGEN - 1)));
		bool selfingNOW = false;
		bool relskewingNOW = false;

		const bool attopnow = (genwidth == HAPLOTYPING) || founder;

		const int selfval = (flag >> (TYPEBITS + 1)) & SELFMASK;

		if (rootgen)
		{
			selfingNOW = SELFING && selfval;
			relskewingNOW = RELSKEWSTATES;
		}

		MarkerVal selfmarker[2];
		double selfsure[2];

		int upflag2 = -1;
		const int upflag = flag >> 1;
		const int upshift = localshift >> 1;
		int f2s = 0;
		const MarkerVal* const themarker = selfingNOW ? selfmarker : &markerdata[marker].first;
		const double* const themarkersure = selfingNOW ? selfsure : &markersure[marker].first;
		int f2end = 2;

		if (flag99 != -1 && genwidth >> (NUMGEN - NUMFLAG2GEN) > 0)
		{
			upflag2 = flag99 >> 1;
			  f2s = (flag99 & 1);
			    f2end = (flag99 & 1) + 1;
		}

		if (relskewingNOW)
		{
		  const int relskewval = (flag >> (BITS_W_SELF + 1));
		    //std::cout << relskewval << std::endl;
			f2s = max(f2s, relskewval);
			f2end = min(f2end, relskewval + 1);
		}

		int firstpar = flag & 1;
		double ok = 0;

		//#pragma ivdep
		// This ivdep is quite complicated, since we actually change markerval, but not in
		// conditions where this alters the results of the other iteration.

		// flag2 determines which value in the tuple to check. flag and localshift determine the haplotype weight value
		// assigned to that test, and which parent to match against that value.
		for (int flag2 = f2s; flag2 < f2end && (HAPLOTYPING || !ok); flag2++)
		{
			int f2n = (flag2 & 1);
			if (selfingNOW)
			{
				int selfindex = (selfval >> 1) ^ f2n;
				selfmarker[0] = UnknownMarkerVal;
				selfmarker[1] = UnknownMarkerVal;
				selfsure[0] = 0;
				selfsure[1] = 0;
				MarkerVal first = markerdata[marker].first;

				if (!markermiss<false>(first, markerdata[marker].second))
				{
					selfmarker[selfindex] = first;
					selfsure[selfindex] = 1 - (1 - markersure[marker].first) * (1 - markersure[marker].second);
				}
				else
				{
					selfmarker[selfindex] = markerdata[marker].second;
					if (markersure[marker].first == 0) { return 0; }
					selfsure[selfindex] = 1 - (markersure[marker].first) * (1 - markersure[marker].second);
				}
			}

			bool allthesame = themarker[0] == themarker[1];
			MarkerVal markerval = inmarkerval;
			double baseval;

			int realf2n = f2n;
			double mainsecondval = 0;
			// If this marker value is not compatible, there is no point in trying.
			if (markermiss<zeropropagate>(markerval, themarker[f2n]))
			{
				baseval = themarkersure[f2n];
				if (themarkersure[f2n] && secondval) mainsecondval = (1.0 - themarkersure[f2n]) * secondval;
			}
			else
			{
				double effectivesecondval = (inmarkerval == UnknownMarkerVal) ? 1 : secondval;				
				baseval = 1.0 - themarkersure[f2n];
				if (themarkersure[f2n] && effectivesecondval) mainsecondval = (themarkersure[f2n] * effectivesecondval);
			}
			 
			// Include it all in one big thing, or 
			if (attopnow)
			{
				baseval += mainsecondval;
				mainsecondval = 0;
			}
			else
			{
				mainsecondval /= baseval;
			}

			if (!baseval) continue;
			bool doupdatehaplo = true;

			// Normalize, in some sense.
			f2n ^= ((firstpar ^ localshift) & 1);

			if (zeropropagate || !genwidth)
			{
				baseval *= 0.5;
				doupdatehaplo = false;
			}
			//			else if (/*!empty &&*/ (allthesame && (CORRECTIONINFERENCE) || (themarker[0] == UnknownMarkerVal && themarker[1] == UnknownMarkerVal && themarkersure[0] + themarkersure[1] == 0)))
			else if (/*!empty &&*/ ((!relskewingNOW) && allthesame && ((CORRECTIONINFERENCE) || (themarkersure[0] == themarkersure[1]))) || selfingNOW)
			{
				baseval *= ((f2n) ? 1.0 : 0.0);
				doupdatehaplo = false;
			}
			else
			{
				if (HAPLOTYPING)
				{
				  // TODO interaction relskew/selfing
					baseval *= fabs((f2n ? 1.0 : 0.0) - haploweight[marker]);
				}
				else
				{
					// No haplotype weights, all interpretations allowed.
					baseval *= 0.5;
				}
			}

			if (!baseval)
			{
				continue;
			}

			// If we are at maximum depth, by depth limit or by lack of ancestor
			if (attopnow || !pars[firstpar])
			{
				// TODO: If pars[firstpar] exists and is empty, things are messy
				// The empty one could, in turn, have parents with genotypes.
				if (zeropropagate && extparams.gstr)
				{
					*(extparams.gstr) += (themarker[realf2n] == (2 * MarkerValue));
				}
			}

			// There should be some other flag for the actual search depth
			if (attopnow)
			{
			}
			else
			{
				// Track to the parent generation, creating evaluation objects first
				// These do a lookup in a special hash for combinations known to be 0, to avoid unnecessary calls
				// Both are checked before either branch of the pedigree is actually traced.
				recursetrackpossible<update, zeropropagate> subtrack1 =
					recursetrackpossible<update, zeropropagate>(this, tb, markerval, mainsecondval, marker,
						upflag,
						upflag2,
						upshift,
						genwidth,
						f2n,
						firstpar,
						0,
						extparams);

				if (subtrack1.prelok && (!zeropropagate || rootgen) && !(update & GENOS))
				{
					double secsecondval = 0;
					if (themarkersure[!realf2n])
					{
						baseval *= (1 - themarkersure[!realf2n]);
						secsecondval = themarkersure[!realf2n] / (1 - themarkersure[!realf2n]);
					}
					baseval *= recursetrackpossible<update, zeropropagate>(this, tb, themarker[!realf2n], secsecondval,
						marker,
						upflag,
						upflag2,
						upshift,
						genwidth,
						f2n,
						!firstpar,
						1,
						extparams);

					if (selfingNOW && extparams.gstr && !(selfval & (1 << !firstpar))) *extparams.gstr = 0;
				}
				if (!baseval) continue;

				if (!(selfingNOW && extparams.gstr && !(selfval & (1 << firstpar))))
					baseval *= subtrack1;
			}

			if (!baseval) continue;

			ok += baseval;

			if ((update & HAPLOS) /*&& !allthesame*/ && doupdatehaplo)
			{
				(*tb.haplos)[n][f2n] += extparams.updateval;				
			}
			if (update & GENOS)
			{
				(*tb.infprobs)[n][realf2n][markerval] += extparams.updateval;
			}
		}
		if (selfingNOW && extparams.gstr) *extparams.gstr *= 2;
		return ok;
	}


	// calltrackpossible is a slight wrapper that hides at least some of the internal parameters needed for the recursion from outside callers
	template<int update, bool zeropropagate> double calltrackpossible(const threadblock& tb, const MarkerVal* const markervals, const unsigned int marker,
		const int genotype, const unsigned int offset, const int flag2, const double updateval = 0.0)
	{
		return trackpossible<update, zeropropagate>(tb, UnknownMarkerVal, 0,
			marker, genotype * 2, flag2, *(tb.shiftflagmode), trackpossibleparams(updateval, 0));
	}

	// "Fix" parents, i.e. infer correct marker values for any zero values existing.
	// Depending on some specifics, this can either assume that all existing marker values will be reflected in some offspring individuals,
	// or indeed only infer specifically those things that can be said from certain, i.e. if a zero value exists, and the other branch of
	// the pedigree doesn't allow a marker value found in the offspring, then it can be assumed to have been transmitted from the individual
	// with the zero value. 
	void fixparents(unsigned int marker, bool latephase)
	{
		MarkerValPair& themarker = markerdata[marker];
		individ* parp[3] = { pars[0], pars[1], this };

#pragma omp critical (parmarkerval)
		for (int i = 0; i < 3; i++)
		{
			if (parp[i])
			{
				while (parp[i]->markervals.size() <= marker)
				{
					parp[i]->markervals.resize(markerdata.size());
				}
			}
		}

		double okvals[2] = { 0 };
		// We only accept an interpretation when it is by exclusion the only possible one. As soon as one intepretation has gained acceptance,
		// no need to retest it.
		for (shiftflagmode = 0; shiftflagmode < 1; shiftflagmode += 1)
		{
			threadblock tb;
			for (int i = 0; i < NUMTYPES; i++)
			{
				for (int flag2 = 0; flag2 < NUMPATHS; flag2++)
				{
					if (okvals[flag2 & 1]) continue;

					double ok = calltrackpossible<false, false>(tb, 0, marker, i, 0, flag2);
					if (!ok) continue;
					okvals[flag2 & 1] += ok;
					if (okvals[0] && okvals[1]) break;
				}
				if (okvals[0] && okvals[1]) break;
			}
		}

		if (!okvals[0] && !okvals[1])
		{
			fprintf(stderr, "Clearing %d:%d (was %d,%d)\n", this->n, marker, markerdata[marker].first.value(), markerdata[marker].second.value());
			markerdata[marker] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
			markersure[marker] = make_pair(0.0, 0.0);
		}

		if ((((bool)okvals[0]) ^ ((bool)okvals[1])) || latephase)
		{
			for (int flag2 = 0; flag2 < 2; flag2++)
			{
				if (!okvals[flag2]) continue;

				for (int k = 0; k < 2; k++)
				{
					if (pars[k])
					{
						int u = ((k ^ flag2) & 1);
						if (latephase || (&themarker.first)[u] != UnknownMarkerVal)
#pragma omp critical (parmarkerval)
						{
							double old1 = 1;
							int old2 = 0;
							auto i = pars[k]->markervals[marker].find((&themarker.first)[u]);
							if (i != pars[k]->markervals[marker].end())
							{
								old1 = i->second.second;
								old2 = i->second.first;
							}
							double probit = (markersure[marker].first + markersure[marker].second);
							probit /= (1 - probit);
							pars[k]->markervals[marker][(&themarker.first)[u]] = make_pair(old2 + 1, old1 * probit);
						}
					}
				}
			}
		}
	}

	// Verifies that an update should take place and performs it. Just a wrapper to calltrackpossible.
	void updatehaplo(const threadblock& tb, const unsigned int marker, const unsigned int i, const int flag2, const double updateval)
	{
		MarkerValPair& themarker = markerdata[marker];

		double ok = calltrackpossible<false, false>(tb, &themarker.first, marker, i, 0, flag2, updateval);

		if (ok)
		{
 			calltrackpossible<true, false>(tb, &themarker.first, marker, i, 0, flag2, updateval);
		}
		else
		{

		}
	}

	// Adjust the probability, i.e. filter all probability values based on the haplotype weights and overall admissibility
	// for the different states.
	void adjustprobs(const threadblock& tb, PerStateArray<double>::T& probs, const unsigned int marker, double& factor,
		const bool oldruleout, int flag99)
	{
		double sum = 0;
		//PerStateArray<double>::T probs2;

		const MarkerValPair& themarker = markerdata[marker];
		const bool ruleout = true; // TODO

		for (int q = 0; q <= (int)!ruleout; q++)
		{
			int f2start = 0;
			int f2end = NUMPATHS;

			// Negative values other than -1 just cause trouble

			// Can we really trust the logic to expand correctly even for zero values?
			//if (flag99 >= 0 || !somezero[marker])
			{
				f2start = flag99;
				f2end = flag99 + 1;
			}

			for (unsigned int i = 0; i < NUMTYPES; i++) // genotype
			{
			  double val = probs[i];

				// We will multiply this already small number with an even smaller number... let's assume it's zero and be done with it.
				if (val < 1e-200)
				{
					probs[i] = 0;
					continue;
				}

				double realok = 0;

				for (int flag2 = f2start; flag2 < f2end; flag2++)
				{
					realok += calltrackpossible<false, false>(tb, &themarker.first, marker, i, 0, flag2);
				}

				
 				// TODO UGLY CAPPING
				if (HAPLOTYPING)
				{
					val *= (double)realok;
				}
				else
				{
					val *= (bool)realok;
				}
				sum += val;
				probs[i] = val;
			}

			if (!ruleout && sum == 0)
			{
			  assert(false && "Bring back probs2");
				for (int i = 0; i < NUMTYPES; i++)
				{
				  //					probs[i] = probs2[i];
				}

				// The code sees: This doesn't make sense, ignore this marker!
				// NOTE: If the parent-grandparent genotypes are completely incompatible, the problem
				// is NOT eliminated.

				// Current caching scheme does not execute full elimination. Efficiency gains possible by not allowing
				// for genotyping error in every single marker. A small epsilon should otherwise be introduced whenever
				// any probs[i] is zero.
				markerdata[marker] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
				fprintf(stderr, "Error in %x, marker %d, impossible path\n", this, marker);
				sum = 1;
			}
			else
				break;
		}

		if (sum <= 0)
		{
			factor = MINFACTOR;
		}
		else
		{
			// Normalize, update auxiliary exponent
			for (int i = 0; i < NUMTYPES; i++)
			{
				probs[i] /= sum;
			}
			factor += log(sum);
		}
	}

	// Append a "multi-step" transition. If the cached values (essentially a N * N transition matrix for the steps from startmark to
	// endmark) are missing, calculate them first.
#if !DOFB
	double fillortake(const threadblock& tb, const int index, const unsigned int startmark, const unsigned int endmark, PerStateArray<double>::T& probs)
	{
		if ((tb.done[*(tb.shiftflagmode)])[index] != (*tb.generation))
		{

			for (unsigned int i = 0; i < NUMTYPES; i++)
			{
				PerStateArray<double>::T probs2 = { {0} };
				probs2[i] = 1;

				// Note ruleout here, could be set to false if we "prime" with an evenly distributed probs first
				(tb.factors[*tb.shiftflagmode])[index][i] = quickanalyze<false, noneturner>(tb, none, startmark,
					endmark,
					NONESTOP,
					-1,
					true,
					probs2);

				double sum = 0;

#pragma ivdep:back
				for (int j = 0; j < NUMTYPES; j++)
				{
					if (!_finite(probs2[j])) probs2[j] = 0.0;
					sum += probs2[j];
				}

				sum *= NUMTYPES;
				if (sum <= 0)
				{
					(tb.factors[*tb.shiftflagmode])[index][i] = MINFACTOR;
					sum = 0;
				}
				else
				{
					(tb.factors[*tb.shiftflagmode])[index][i] += log(sum);
					sum = 1 / sum;
				}
				double* probdest = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];

#pragma ivdep
				for (int j = 0; j < NUMTYPES; j++)
				{
					probdest[j] = probs2[j] * sum;
				}
			}
			(tb.done[*tb.shiftflagmode])[index] = (*tb.generation);
		}

		double factor = MINFACTOR;
		for (int i = 0; i < NUMTYPES; i++)
		{
			if (probs[i] == 0.0) continue;

			factor = max(factor, (tb.factors[*tb.shiftflagmode])[index][i]);
		}

		PerStateArray<double>::T probs2 = { {0} };

		for (int i = 0; i < NUMTYPES; i++)
		{
			double step = (tb.factors[*tb.shiftflagmode])[index][i] - factor;
			if (probs[i] == 0.0 || step <= -100.0f) continue;
			double basef = exp((double)step) * probs[i];
			if (basef == 0.0) continue;
			const double* probsource = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];
#pragma ivdep
			for (int j = 0; j < NUMTYPES; j++)
			{
				const double term = basef * probsource[j];
				probs2[j] += term;
			}
		}

		double sum = 0;

		for (int i = 0; i < NUMTYPES; i++)
		{
			sum += probs2[i];
		}

		if (sum <= 0)
		{
			factor = MINFACTOR;
			sum = 0;
		}
		else
		{
			factor += log(sum);
			sum = 1 / sum;
		}

		for (int i = 0; i < NUMTYPES; i++)
		{
			probs[i] = probs2[i] * sum;
		}


		return factor;
	}

	// Analyze for a specific range, including a possible fixed specific state at some position (determined by stopdata)
	template<bool inclusive, class T, class G> double quickanalyze(const threadblock& tb, const T& turner, unsigned int startmark,
		const unsigned int endmark, const G &stopdata, const int flag2, bool ruleout, PerStateArray<double>::T& probs,
		float minfactor = MINFACTOR)
	{
		unsigned int stepsize;
		double factor = 0;
		bool allowfull = inclusive;
		bool frommem = false;
		if (inclusive && tb.quickgen[*tb.shiftflagmode] == *tb.generation && firstlockmatch(tb.lockpos[*tb.shiftflagmode], stopdata))
		{
			allowfull = true;
			frommem = true;
			factor = tb.quickfactor[*tb.shiftflagmode];
			if (factor <= minfactor)
			{
				return MINFACTOR;
			}

			startmark = tb.quickmark[*tb.shiftflagmode];
			//			printf("STARTING FROM MARKER: %d\n", startmark);
			for (int i = 0; i < NUMTYPES; i++)
			{
				probs[i] = tb.quickmem[*tb.shiftflagmode][i];
			}
		}

		// Loop until we have reached the end, with blocks of varying sizes.
		while (startmark < endmark)
		{
			for (stepsize = 1; stepsize < (endmark - startmark + allowfull) &&
				stopdata.okstep(startmark, startmark + stepsize) &&
				!(startmark & (stepsize - 1)); stepsize *= 2);

				// A single step, either due to the fixated genotypee being within this range, or simply because we've gone all the
				// way down the tree.
				if (stepsize <= 2)
				{
					stepsize = 1;

					if (!frommem && !stopdata.okstep(startmark, startmark + 1))
					{
						// If we have a fixated genotype at marker x in one call, it is very likely that the next call will also
						// be a fixated genotype at marker x. Possibly another one, but still fixated. The fixated genotype does not
						// change the probability values leading up to this position, so they are cached.
						tb.quickgen[*tb.shiftflagmode] = *tb.generation;
						tb.quickmark[*tb.shiftflagmode] = startmark;
						tb.lockpos[*tb.shiftflagmode] = (int)stopdata;
						tb.quickfactor[*tb.shiftflagmode] = factor;

						for (int i = 0; i < NUMTYPES; i++)
						{
							tb.quickmem[*tb.shiftflagmode][i] = probs[i];
							tb.quickendmarker[*tb.shiftflagmode][i] = -1;
							tb.quickendfactor[*tb.shiftflagmode][i] = 1.0;
						}
						frommem = true;
					}
					bool willquickend = (frommem && turner.canquickend() && canquickend(startmark, stopdata));
					int genotype = stopdata.getgenotype(startmark);

					if (genotype < 0) willquickend = false;

					if (willquickend)
					{
						if (factor + tb.quickendfactor[*tb.shiftflagmode][genotype] <= minfactor)
						{
							return MINFACTOR;
						}
						// If we are doing a quick end
						factor += realanalyze<4, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
					}
					else
					{
						factor += realanalyze<0, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
					}

					if (!_finite(factor) || factor <= minfactor)
					{						
						if (!_finite(factor)) fprintf(stderr, "Non-finiteA %d %d\n", n, startmark);
						return MINFACTOR;
					}
					if (willquickend)
					{
						double wasval = probs[genotype];
						if (tb.quickendfactor[*tb.shiftflagmode][genotype] > 0.0 || tb.quickendmarker[*tb.shiftflagmode][genotype] != startmark)
						{
							// Precalc the full set of transitions from this very marker, for all states, all the way to the end
							// That way, any new analysis starting from this marker can be solved with no recomputation at all.
							// Note that only the latest marker used is stored here, but as we generally loop over all states and possibly
							// multiple flag values, it still helps.
							for (int i = 0; i < NUMTYPES; i++)
							{
								tb.quickendprobs[*tb.shiftflagmode][genotype][i] = 0;
							}
							tb.quickendprobs[*tb.shiftflagmode][genotype][genotype] = 1.0;
							double superfactor = realanalyze<0 | 2, T>(tb, turner, startmark, startmark + stepsize, NONESTOP, -1, ruleout, &tb.quickendprobs[*tb.shiftflagmode][genotype]);

							superfactor += quickanalyze<true, T>(tb, turner, startmark + stepsize, endmark, NONESTOP, flag2, ruleout, tb.quickendprobs[*tb.shiftflagmode][genotype],
								// all uses of this precalced data will have the
								// quickfactor component in common, so taking that
								// into account for the limit is not problematic
								// any strong filtering at this specific locus can be
								// handled by the return MINFACTOR line above
								minfactor - tb.quickfactor[*tb.shiftflagmode]);

							tb.quickendmarker[*tb.shiftflagmode][genotype] = startmark;
							tb.quickendfactor[*tb.shiftflagmode][genotype] = superfactor;
						}

						factor += tb.quickendfactor[*tb.shiftflagmode][genotype];
						factor += log(wasval);
						if (factor <= minfactor) return MINFACTOR;

						for (int i = 0; i < NUMTYPES; i++)
						{
							probs[i] = tb.quickendprobs[*tb.shiftflagmode][genotype][i];
						}

						return factor;
					}
				}
				else
				{
					// Use a stored multi-step transition.
					stepsize /= 2;

					int base = 0;
					for (int q = stepsize / 2; q > 1; q /= 2)
					{
						base += markerposes.size() / q;
					}

					base += startmark / stepsize;

					factor += fillortake(tb, base, startmark, startmark + stepsize, probs);
				}
			startmark += stepsize;
			allowfull |= true;

			if (!_finite(factor) || factor <= minfactor)
			{
				if (!_finite(factor)) fprintf(stderr, "Non-finiteB %d %d\n", n, startmark);
				if (!frommem && !stopdata.okstep(startmark, endmark))
				{
					tb.quickgen[*tb.shiftflagmode] = *tb.generation;
					tb.lockpos[*tb.shiftflagmode] = (int)stopdata;
					tb.quickfactor[*tb.shiftflagmode] = MINFACTOR;
				}

				return MINFACTOR;
			}
		}

		return factor;
	}
#else
// Analyze for a specific range, including a possible fixed specific state at some position (determined by stopdata)
	template<bool inclusive, class T, class G> double quickanalyze(const threadblock& tb, const T& turner, unsigned int startmark,
		const unsigned int endmark, const G &stopdata, const int flag2, bool ruleout, PerStateArray<double>::T& probs,
		float minfactor = MINFACTOR)
	{
		// Probe the distance to test
		int newstart = startmark;
		int origstart = startmark;
		bool allowfull = inclusive;

		while (stopdata.okstep(startmark, newstart + 1))
		{
			int stepsize;
			for (stepsize = 1; stepsize < (endmark - newstart + allowfull) &&
				stopdata.okstep(newstart, newstart + stepsize); stepsize *= 2)
			{
			}

			if (stepsize == 1) break;

			stepsize /= 2;
			newstart += stepsize;
		}

		double factor = tb.fwbwfactors[*tb.shiftflagmode][newstart][0];
		probs = tb.fwbw[*tb.shiftflagmode][newstart][0];
		startmark = newstart;

		if (factor < minfactor) return factor;

		while (startmark < endmark)
		{
			int stepsize = 1;

			bool willquickend = /*(turner.canquickend() && canquickend(startmark, stopdata))*/ stopdata.okstep(startmark + 1, endmark);
			int genotype = stopdata.getgenotype(startmark);

			if (willquickend)
			{
				// If we are doing a quick end
				factor += realanalyze<4, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
				
				// We might have turned, this value might not exist
				initfwbw(tb, origstart, endmark);

				double sum = 0;

				if (factor < minfactor) return factor;

				for (int k = 0; k < NUMTYPES; k++)
				{
					probs[k] *= tb.fwbw[*tb.shiftflagmode][startmark][1][k];
					sum += probs[k];
				}

				factor += tb.fwbwfactors[*tb.shiftflagmode][startmark][1];

				if (sum <= 0)
				{
					factor = MINFACTOR;
				}
				else
				{
					// Normalize, update auxiliary exponent
					sum = 1 / sum;
					for (int i = 0; i < NUMTYPES; i++)
					{
						probs[i] *= sum;
					}
					factor -= log(sum);
				}

				return factor;
			}
			else
			{
				std::cerr << "Incomplete at " << startmark << std::endl;
				factor += realanalyze<0, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
			}

			startmark += stepsize;
			allowfull |= true;
		}

		return factor;
	}
#endif

#if !DOFB
	// A wrapper to quickanalyze, preparing the start probability vector.
	template<class T, class G> double doanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const G& stopdata,
		const int flag2, bool ruleout = false, PerStateArray<double>::T* realprobs = 0, float minfactor = MINFACTOR)
	{
		PerStateArray<double>::T fakeprobs;
		PerStateArray<double>::T& probs = realprobs ? *realprobs : fakeprobs;

		if (realprobs)
		{
			//probs = realprobs;
		}
		else
		{
			//probs = fakeprobs;
			double selfingfactors[4];
			if (SELFING)
			{
				int selfgen = gen - 2;
				selfingfactors[0] = 1.0 / (1 << selfgen);
				selfingfactors[1] = (1 - selfingfactors[0]) * 0.5;
				selfingfactors[2] = (1 - selfingfactors[0]) * 0.5;
				selfingfactors[3] = 0;
			}
			for (int i = 0; i < NUMTYPES; i++)
			{
				fakeprobs[i] = EVENGEN * (SELFING ?
					selfingfactors[i >> TYPEBITS] : 1.0);
			}
		}

		double factor = quickanalyze<true, T>(tb, turner, startmark, endmark, stopdata, flag2, ruleout, probs, minfactor);
		bool small = !_finite(factor) || minfactor >= factor;

		if (!small) adjustprobs(tb, probs, endmark, factor, ruleout, -1); // TODO, the very last marker can be somewhat distorted

		return factor;
	}
#else
	void initfwbw(const threadblock& tb, const int startmark, const int endmark)
	{
		if (tb.fwbwdone[*(tb.shiftflagmode)] != *(tb.generation))
		{
			PerStateArray<double>::T probs;

			// Initialize forward-backward matrices in one big go.
			//probs = fakeprobs;
			double selfingfactors[4];
			if (SELFING)
			{
				int selfgen = gen - 2;
				selfingfactors[0] = 1.0 / (1 << selfgen);
				selfingfactors[1] = (1 - selfingfactors[0]) * 0.5;
				selfingfactors[2] = (1 - selfingfactors[0]) * 0.5;
				selfingfactors[3] = 0;
			}
			for (int i = 0; i < NUMTYPES; i++)
			{
				probs[i] = EVENGEN * (SELFING ?
					selfingfactors[(i >> TYPEBITS) & SELFMASK] : 1.0);
			}

			realanalyze<ANALYZE_FLAG_STORE | ANALYZE_FLAG_FORWARD | 1, noneturner>(tb, noneturner(), startmark, endmark, NONESTOP, -1, false, &probs);

			for (int i = 0; i < NUMTYPES; i++)
			{
				probs[i] = 1.0;
			}
			realanalyze<ANALYZE_FLAG_STORE | ANALYZE_FLAG_BACKWARD | 1, noneturner>(tb, noneturner(), startmark, endmark, NONESTOP, -1, false, &probs);

			tb.fwbwdone[*(tb.shiftflagmode)] = *(tb.generation);
		}
	}

	template<class T, class G> double doanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const G& stopdata,
		const int flag2, bool ruleout = false, PerStateArray<double>::T* realprobs = 0, float minfactor = MINFACTOR)
	{
		if (realprobs != 0) fprintf(stderr, "THIS WAS NOT EXPECTED. WHEN IS THIS USED?\n");
		PerStateArray<double>::T probs;

		initfwbw(tb, startmark, endmark);

		return quickanalyze<true, T>(tb, turner, startmark, endmark, stopdata, flag2, ruleout, probs, minfactor);
	}
#endif

	// This is the actual analyzing code. It works with no caches, and can function independently, but is generally only used to patch in those
	// elements not found in caches by quickanalyze and fillortake.
	//
	// Both transition and emission (through adjustprobs) probabilities are handled here.
	//
	// first bit in updateend signals whether the interval is end-inclusive at endmark
	// the second bit in updateend will be 0 if the interval is end-inclusive at startmark, and 1 IF NOT
	// the third bit will cause the code to quit early, after processing the genotype and turner condition
	// the fourth bit will cause the code to go backwards, if set
	// the fifth bit will cause whe code to store probs in fwbw structures
	// Be careful how things switch meaning when going backwards! Could be refactored!
	template<int updateend, class T, class G> double realanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const G& stopdata,
		const int flag2, const bool ruleout = false, PerStateArray<double>::T* realprobs = 0)
	{
		PerStateArray<double>::T fakeprobs;
		PerStateArray<double>::T& probs = realprobs ? *realprobs : fakeprobs;

		if (realprobs)
		{
			//probs = realprobs;
		}
		else
		{
			//probs = fakeprobs;
			for (int i = 0; i < NUMTYPES; i++)
			{
				fakeprobs[i] = EVENGEN;
			}
		}

		// The *logged* normalization factor
		double factor = 0;

		// Walk over all markers.
		int d = 1;
		int firstmark = startmark;
		int lastmark = endmark;
		if (updateend & ANALYZE_FLAG_BACKWARD)
		{
			d = -1;
			swap(firstmark, lastmark);
#if DOFB
			if ((updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				copy(probs.begin(), probs.end(),
					fwbw[*tb.shiftflagmode][firstmark][(bool)(updateend & ANALYZE_FLAG_BACKWARD)].begin());

				fwbwfactors[*tb.shiftflagmode][firstmark][(bool)(updateend & ANALYZE_FLAG_BACKWARD)] = factor;
			}
#endif
		}

		for (int j = firstmark + d; j != lastmark + d; j += d)
		{
			double startpos = markerposes[j - d];
			double endpos = markerposes[j];
			int genotype = -1;

			if (updateend & ANALYZE_FLAG_BACKWARD) swap(startpos, endpos);

			bool tofind = stopdata.fixtofind(genotype, startpos, endpos, j);

			int f2use = -1;

			if (tofind)
			{
				f2use = flag2;
			}

#if DOFB
			if (!(updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				copy(probs.begin(), probs.end(),
					fwbw[*tb.shiftflagmode][j - d][(bool)(updateend & ANALYZE_FLAG_BACKWARD)].begin());

				fwbwfactors[*tb.shiftflagmode][j - d][(bool)(updateend & ANALYZE_FLAG_BACKWARD)] = factor;
			}
#endif

			if (genotype != -2)
			{
				// If we are at the very first position, and the specific flag was set, include the emission probabilities for
				// the previous marker. Used to maximize the caching.
				if (!((updateend & 2) && (j == firstmark + d))) adjustprobs(tb, probs, j - d, factor, ruleout, f2use);
			}
			else
			{
				// We could do some stuff here to avoid excessive adjustprobs calls
				// a -2 genotype does not only mean that all genotypes are allowed, but indeed that the marker data at this marker
				// is ignored!
			}

			// For a specific intra-marker region, we have two cases: the case of a fixated position between the two markers,
			// and the simple case of no fixated position.
			for (int iter = 0; iter <= (int)tofind; iter++)
			{
				// If iter is 1, we have currently handled the transition all the way to the fixated position. Now filter to keep
				// only a single state value positive.
				if (iter)
				{
					turner(probs);
					if (genotype >= 0)
					{
						for (int i = 0; i < NUMTYPES; i++)
						{
							probs[i] = probs[i] * (i == genotype);
						}
					}

					// Were we asked to stop at this very state, in the middle of things?
					if (updateend & 4) return factor;
				}


				double dist = (endpos - startpos);

				// Compute transitions, assuming there is any distance to transfer over.
				if (dist > 0)
				{
					PerStateArray<double>::T probs2 = { {0} };
					double recprob[2 + SELFING][2];
					// Compute recombination probabilities for this specific distance, for the two sexes.
					// (as the sex-dependent marker distance might not be a simple transformation, the actrec
					// data comes into play).
					const int selfgen = gen - 2;
#pragma ivdep
					for (int gen = 0; gen < 2 + SELFING; gen++)
					{
						for (int k = 0; k < 2; k++)
						{
							recprob[gen][k] = 0.5 * (1.0 - exp((gen == 2 ? selfgen : 1) * getactrec(stopdata, startpos, endpos, k, j - d + 1, gen) * (dist)));
							//					if (iter == tofind) recprob[k] = max(recprob[k], 1e-5);
						}
					}

					// Precompute values for presence (/lack of) a crossover event for either sex.
					double other[2][2][2];
#pragma ivdep
					for (int gen = 0; gen < 2; gen++)
					{
						for (int m = 0; m < 2; m++)
						{
							for (int k = 0; k < 2; k++)
							{
								double prob = recprob[gen][k];
								if (m) prob = 1.0 - prob;

								other[gen][m][k] = prob;
							}
						}
					}

					std::array<double, NONSELFNUMTYPES> recombprec;

#pragma ivdep
					for (int index = 0; index < NONSELFNUMTYPES; index++)
					{
						recombprec[index] = 1;
					}

					double selfprec[4 * SELFING + 1][4 * SELFING + 1] = { 0 };
					if (SELFING)
					{
						selfprec[0][1] = selfprec[0][2] = recprob[2][0];
						selfprec[0][0] = 1 - 2 * recprob[2][0];
						selfprec[1][0] = selfgen ? selfprec[0][1] * 2.0 / ((1 << selfgen) - 1) : 1;
						selfprec[1][2] = selfprec[1][0] * selfprec[0][1]; // TODO, UNDERESTIMATED!
						selfprec[2][0] = selfprec[1][0];
						selfprec[2][1] = selfprec[1][2];
						selfprec[2][2] = selfprec[1][1] = 1 - selfprec[1][0] - selfprec[1][2];
					}

					// Compute probabilities for arbitrary xor values of current and future state
					for (int t = 0; t < TYPEBITS; t++)
					{
						int sex = TYPESEXES[t];
						int gen = TYPEGENS[t];
#pragma ivdep
						for (int index = 0; index < NONSELFNUMTYPES; index++)
						{
							int val = !((index >> t) & 1);
							recombprec[index] *= other[gen][val][sex];
						}
					}

					float relscore[2] = { 1, 1 };
					if (RELSKEWSTATES && iter == tofind)
					{
					  relscore[0] = relhaplo[min(j, j - d)];
					  relscore[1] = 1 - relhaplo[min(j, j - d)];
					}

					// Use those xor values
					// For the 4-state model, this is an inefficient way to go about it, but it is quite a bit more efficient for
					// the 64-state model (or beyond).
					for (int from = 0; from < VALIDSELFNUMTYPES; from++) // SELFING condition removes the double-bit set case, which is not allowed
					{
						if (probs[from] < MINFACTOR || !probs[from]) continue;
#pragma ivdep
						for (int to = 0; to < VALIDSELFNUMTYPES; to++)
						{
							int xored = from ^ to;
							probs2[to] += probs[from] * recombprec[xored & (NONSELFNUMTYPES - 1)] * (SELFING ? selfprec[(from >> TYPEBITS) & SELFMASK][(to >> TYPEBITS) & SELFMASK] : 1) *
								(RELSKEWSTATES ? relscore[(xored >> BITS_W_SELF) & 1] : 1);
						}
					}

					for (int c = 0; c < NUMTYPES; c++)
					{
						probs[c] = probs2[c];
					}
				}
				//				else
				{
					//				if (iter == tofind)
					{
						/*for (int c = 0; c < VALIDSELFNUMTYPES; c++)
						{
							if (probs[c] < 1e-200) probs[c] = 1e-200;
						}*/
					}
				}

				// startpos and endpos are always defined in the forward sense, i.e. startpos being upstream of endpos
				if (updateend & ANALYZE_FLAG_BACKWARD)
				{
					endpos = startpos;
					startpos = markerposes[j];
				}
				else
				{
					startpos = endpos;
					endpos = markerposes[j];
				}
			}

#if DOFB
			if ((updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				copy(probs.begin(), probs.end(),
					fwbw[*tb.shiftflagmode][j][(bool)(updateend & ANALYZE_FLAG_BACKWARD)].begin());

				fwbwfactors[*tb.shiftflagmode][j][(bool)(updateend & ANALYZE_FLAG_BACKWARD)] = factor;
			}
#endif
		}

		if (updateend & 1)
		{
			adjustprobs(tb, probs, lastmark, factor, ruleout, -1); // TODO
#if DOFB
			if (!(updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				copy(probs.begin(), probs.end(),
					fwbw[*tb.shiftflagmode][lastmark][(bool)(updateend & ANALYZE_FLAG_BACKWARD)].begin());

				fwbwfactors[*tb.shiftflagmode][lastmark][(bool)(updateend & ANALYZE_FLAG_BACKWARD)] = factor;
			}
#endif
		}

		return factor;
	}
};

// Oh, how we waste memory, in a pseudo-O(1) manner
individ* individer[INDCOUNT];
// dous contains those individuals that really should be analyzed
vector<individ*> dous;

// retrieve the individual with a specific number
#ifdef _MSC_VER
individ* const getind(int n, bool real = false)
#else
individ* const __attribute__((const)) getind(int n, bool real = false)
#endif
{
	if (!real && !individer[n]) n = 0;
	if (n <= 0) return 0;

	if (!individer[n])
	{
		printf("Creating %d\n", n);
		individer[n] = new individ();
		individer[n]->n = n;
		individ* ind = individer[n];
		ind->empty = false;
		ind->markerdata.resize(markerposes.size());
		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());

		ind->infprobs.resize(markerposes.size());
		ind->markersure.resize(markerposes.size());

		if (RELSKEWS)
		{
			ind->relhaplo.resize(markerposes.size());
		}

		for (int x = 0; x < markerposes.size(); x++)
		{
			ind->markerdata[x] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
			ind->haplobase[x] = 0;
			ind->haplocount[x] = 0;
			ind->haploweight[x] = 0.5;
			ind->negshift[x] = 0;
			ind->markersure[x] = make_pair(0, 0);
			if (RELSKEWS)
			{
				ind->relhaplo[x] = 0.5;
			}
		}

		//		ind->semishift.resize(5000);
		ind->lastinved.resize(chromstarts.size());
		ind->lockstart.resize(chromstarts.size());

		for (int i = 0; i < chromstarts.size(); i++)
		{
			ind->lastinved[i] = -1;
			ind->lockstart[i] = 0;
		}


	}

	return individer[n];
}

// read qtlmas sample data, code used for testing haplotypes, quite a bit of hardcoding present in this function
void readqtlmas()
{
	FILE* inpheno = fopen("phenotype.txt", "r");
	FILE* ingeno = fopen("genotype_cor.txt", "r");
	markertranslation.resize(6000);
	for (int i = 0; i < 6; i++)
	{
		chromstarts.push_back(i * 1000);
		for (int j = 0; j < 1000; j++)
		{
			markertranslation[i * 1000 + j] = i * 1000 + j;
			markerposes.push_back(j * 0.01 * 2);
			for (int t = 0; t < 2; t++)
			{
				actrec[t].push_back(baserec[t]);
			}
		}
	}
	chromstarts.push_back(6000);

	char tlf[16384];
	fgets(tlf, 16384, inpheno);

	for (int indn = 1; indn < 7000; indn++)
	{
		int fa = 0;
		int mo = 0;
		individ* ind = getind(indn);
		ind->pars[0] = getind(fa);
		ind->pars[1] = getind(mo);

		ind->gen = 5;

		ind->sex = 0;
		ind->strain = 1;

		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());
		ind->markerdata.resize(markerposes.size());

		for (int i = 0; i < 6000; i++)
		{
			ind->haploweight[i] = 0.0;
		}
	}

	while (fgets(tlf, 16384, inpheno))
	{
		int indn, fa, mo, sex, gen;
		if (sscanf(tlf, "%d %d %d %d %d", &indn, &fa, &mo, &sex, &gen) < 4) break;


		individ* ind = getind(indn);

		ind->pars[0] = getind(fa);
		ind->pars[1] = getind(mo);


		if (fa && mo && (!ind->pars[0]->haplobase.size() ||
			!ind->pars[1]->haplobase.size()))
		{
			printf("PROBLEM %d %d %d\n", indn, fa, mo);
		}
		ind->gen = gen;

		ind->sex = sex - 1;
		ind->strain = (indn % 2) + 1;

		if (gen >= 2) dous.push_back(ind);

		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());
		ind->markerdata.resize(markerposes.size());
	}

	int indn;
	while (fscanf(ingeno, "%d", &indn) == 1 && indn)
	{
		individ* ind = getind(indn);
		if (!ind->markerdata.size())
		{
			printf("Stopping at %d\n", indn);
			break;
		}

		for (int i = 0; i < 6000; i++)
		{
			int a, b;
			fscanf(ingeno, "%d %d", &a, &b);
			if (i % 10)
			{
				a = 0;
				b = 0;
			}

			ind->markerdata[i] = make_pair(a * MarkerValue, b * MarkerValue);

		}
	}
}

void readqtlmas14()
{
	FILE* ingeno = fopen("genotypes.txt", "r");
	FILE* inminfo = fopen("marker-info.txt", "r");
	FILE* inped = fopen("pedigree.txt", "r");
	markertranslation.resize(10033);
	int n, c, bppos;
	int n2 = 0;
	int oldc = -1;

	fprintf(stderr, "Input: %d %d %d\n", ingeno, inminfo, inped);

	while (fscanf(inminfo, "%d %d %d", &n, &c, &bppos) == 3)
	{
		markertranslation[n - 1] = n - 1;
		if (c != oldc)
		{
			chromstarts.push_back(n - 1);			
			n2 = 0;
			oldc = c;
		}
		markerposes.push_back(bppos / 1000000.0);
		for (int t = 0; t < 2; t++)
		{
			actrec[t].push_back(baserec[t]);
		}
		n2++;
	}
	chromstarts.push_back(n);
	markertranslation[n] = n;
	markertranslation[n + 1] = n + 1;
	markerposes.push_back(0);
	markerposes.push_back(1);
	for (int t = 0; t < 2; t++)
	{
		actrec[t].push_back(baserec[t]);
	}
	for (int t = 0; t < 2; t++)
	{
		actrec[t].push_back(baserec[t]);
	}



	int fa, mo;
	int indn;
	char sex[3];
	while (fscanf(inped, "%d %d %d %s", &indn, &fa, &mo, sex) == 4)
	{
		individ* ind = getind(indn, true);
		ind->pars[0] = getind(fa, true);
		ind->pars[1] = getind(mo, true);

		ind->gen = 0;
		ind->gen += indn > 20;
		ind->gen += indn > 438;
		ind->gen += indn > 1349;
		ind->gen += indn > 2326;

		ind->n = indn;

		if (ind->gen > 0) dous.push_back(ind);

		ind->sex = (sex[0] == 'F');
		ind->strain = 1;
		ind->markerdata.resize(markerposes.size());

		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());

		ind->infprobs.resize(markerposes.size());
		ind->markersure.resize(markerposes.size());
		//		ind->semishift.resize(5000);
		ind->lastinved.resize(chromstarts.size());
		ind->lockstart.resize(chromstarts.size());

		for (int i = 0; i < chromstarts.size(); i++)
		{
			ind->lastinved[i] = -1;
			ind->lockstart[i] = 0;
		}

		for (int i = 0; i < markerposes.size(); i++)
		{
			ind->haploweight[i] = 0.5;
		}
	}

	int indread = 0;
	while (fscanf(ingeno, "%d", &indn) == 1 && indn)
	{
		individ* ind = getind(indn, true);
		indread++;
		if (!ind->markerdata.size())
		{
			printf("Stopping at %d\n", indn);
			break;
		}

		for (int i = 0; i < /*markerposes.size() - 2*/ 10031; i++)
		{
			int a, b;
			fscanf(ingeno, "%d %d", &a, &b);
			if (a > 5 || b > 5) fprintf(stderr, "Marker offset mismatch: %d %d %d %d\n", indn, i, a, b);

			if (i < markerposes.size())
			{
				ind->markerdata[i] = make_pair(a * MarkerValue, b * MarkerValue);

				if (ind->gen < 4)
				{
					ind->markerdata[i] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
				}
				else
					if (i == 1181) fprintf(stderr, "Input check: %d %d %d\n", indn, a, b);
			}

		}

		//		if (indn > 2360) break;
	}
	fprintf(stderr, "Individuals read: %d\n", indread);
	getind(1633)->markerdata[4] = make_pair(2 * MarkerValue, 2 * MarkerValue);
	getind(1351)->markerdata[4] = make_pair(2 * MarkerValue, 2 * MarkerValue);
}


// read marker info in a format similar to that accepted by ccoeff
void readmarkerinfo(FILE* in)
{
	int n, m;
	fscanf(in, "%d %d", &n, &m);
	vector<int> count;
	markertranslation.resize(m);

	for (int i = 0, j = 0; i < n; i++)
	{
		fscanf(in, "%d", &m);
		count.push_back(m);

		printf("Number on ch %d: %d \n", i + 1, m);

		for (int k = 0; k < count.end()[-1]; k++)
		{
			fscanf(in, "%d", &m);
			markertranslation[m - 1] = ++j;
		}
	}

	int pos = 0;
	for (int i = 0; i < n; i++)
	{
		chromstarts.push_back(pos);

		vector<double> part[2];
		for (int t = 0; t < sexc; t++)
		{
			fscanf(in, "%d", &m);

			double j = 0;
			for (int k = 0; k < count[i]; k++)
			{
				double v;
				fscanf(in, "%lf", &v);
				j += v;

				for (int p = 0; p < 2 / sexc; p++)
				{
					part[t + p].push_back(j / discstep);
				}
			}
		}

		for (int k = 0; k < count[i]; k++)
		{
			double sum = 0;
			for (int t = 0; t < 2; t++)
			{
				sum += part[t][k];
			}
			sum /= 2.0;

			markerposes.push_back(sum);
			printf("%lf\n", sum);

			for (int t = 0; t < 2; t++)
			{
				double avgdist = 0;
				if (k)
				{
					avgdist = markerposes[pos] - markerposes[pos - 1];
				}

				if (avgdist)
				{
					actrec[t].push_back(baserec[t] * (part[t][k] - part[t][k - 1]) / avgdist);
				}
				else
				{
					actrec[t].push_back(-1.0);
				}
			}
			pos++;
		}
	}
	chromstarts.push_back(pos);

	printf("%d chromosomes, %d markers, really %d/%d\n", n, m, markerposes.size(), chromstarts.size() - 1);
}

// read a pedigree in a format similar to the one accepted by ccoeff
void readped(FILE* in)
{
	int famsize;

	// The format is a set of full sibships, each started by four founders (possibly some identical), two parents
	while (fscanf(in, "%d", &famsize) == 1)
	{
		for (int i = 0; i < famsize + 6; i++)
		{
			int indn;
			int fa, mo;
			int sex = 0;

			// strain here could just as well have been called line
			// OTOH, in a line-based file format, animal lines can be confusing
			int strain = -1;
			char tlf[255];
			while (fgets(tlf, 255, in))
			{
				if (sscanf(tlf, "%d %d %d %d %d", &indn, &fa, &mo, &sex, &strain) >= 3) break;
			}

			individ* ind = getind(indn, true);
			ind->pars[1] = getind(fa, true);
			ind->pars[0] = getind(mo, true);

			if (ind->pars[0] && ind->pars[1] && ind->pars[0]->sex == 1 && ind->pars[1]->sex == 0) {
				swap(ind->pars[0], ind->pars[1]);
				printf("%d\n", indn);
			}

			ind->sex = sex - 1;
			ind->strain = strain;
			// The generation is defined by the index, the first few are always the first two generations
			ind->gen = 0;
			ind->gen += i >= 4;
			ind->gen += i >= 6;

			// Only analyze the F2s
			if (i >= 6) dous.push_back(ind);
		}
	}
	printf("Pedigree containing %d F2 individuals\n", dous.size());
}

// Read marker data and dimension several data fields dependent on marker data
// These are expected to be in a format similar to the one accepted by ccoeff
void readmarkerdata(FILE* in)
{
	int indn;
	int n = 0;
	while (fscanf(in, "%d", &indn) == 1)
	{
		individ* ind = getind(indn, true);
		n++;
		ind->markerdata.resize(markerposes.size());
		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());

		for (unsigned int i = 0; i < markerposes.size(); i++)
		{
			ind->haploweight[i] = 0.5;
		}

		for (unsigned int i = 0; i < markertranslation.size(); i++)
		{
			int a, b;
			if (fscanf(in, "%d %d", &a, &b) != 2)
			{
				fprintf(stderr, "Marker mismatch: %d\n", indn);
			}
			if (markertranslation[i])
			{
				ind->markerdata[markertranslation[i] - 1] =
					make_pair(a * MarkerValue, b * MarkerValue);
			}
		}
	}
	printf("Marker data parsed for %d individuals\n", n);
}

void readhaploweights(FILE* in)
{
	int indn;
	int end = markerposes[chromstarts[1] - 1];

	while (fscanf(in, "%d", &indn) == 1)
	{
		// TODO chromcount?
		individ* ind = getind(indn, true);
		double num1, num2;
		int a, b;
		for (int i = chromstarts[0]; i < chromstarts[1]; i++)
		{
			fscanf(in, "%lf %d %d %lf", &num1, &a, &b, &num2);
			ind->haploweight[i] = num1;

			if (ind->haploweight[i] <= 0.5)
			{
				swap(ind->markerdata[i].first, ind->markerdata[i].second);
				ind->haploweight[i] = 1 - ind->haploweight[i];
			}
			if (ind->haploweight[i]) ind->haploweight[i] = 0.5;
		}
		ind->genotypegrid.resize(end + 1);

		if (!ind->pars[0] && !ind->pars[1])
		{
			for (int i = 0; i <= end; i++)
			{
				ind->genotypegrid[i] = (ind->strain - 1) * 3;
			}
		}
	}

	FILE* indgrid = fopen("indgrid.out", "w");
	FILE* indmarkers = fopen("indmarkers.out", "w");
	for (int gen = 0; gen <= 2; gen++)
	{
		for (int i = 0; i < INDCOUNT; i++)
		{
			individ* ind = getind(i, true);

			if ((ind && ind->gen == gen && ind->genotypegrid.size()) /*&& ind->pars[0] && ind->pars[1]*/)
			{
				fprintf(indmarkers, "%d ", i);
				fprintf(indgrid, "%d ", i);

				int phase = myrand(4);
				int marker = 0;

				if (gen)
				{
					for (int marker = 0; marker < chromstarts[1]; marker++)
					{
						if (myrand(2))
						{
							swap(ind->markerdata[marker].first, ind->markerdata[marker].second);
						}
					}
				}


				for (int i = 0; i <= end; i++)
				{
					int gval = 0;
					int startmarker = marker;
					if (gen)
					{
						for (int p = 0; p < 2 && gen; p++)
						{
							marker = startmarker;

							int mask = 1 << p;
							bool parhalf = phase & mask;

							gval |= mask * ((bool)(ind->pars[p]->genotypegrid[i] & (1 << parhalf)));

							for (; markerposes[marker] < i + 1 && marker < chromstarts[1]; marker++)
							{
								if ((&ind->markerdata[marker].first)[p] != UnknownMarkerVal)
								{
									(&ind->markerdata[marker].first)[p] = (&ind->pars[p]->markerdata[marker].first)[parhalf];
								}
							}

							phase &= 65535 - mask;

							int high = 100;
							int lo = 0;

							if (marker < chromstarts[1])
							{
								high = 10000;
								if (actrec[p][marker])
									lo = actrec[p][marker] * -5000;
							}

							if (myrand(high) <= lo) parhalf = !parhalf;

							phase |= parhalf * mask;
						}
						ind->genotypegrid[i] = gval;
					}

					fprintf(indgrid, " %d", ind->genotypegrid[i]);
				}
				for (int marker = 0; marker < chromstarts[1]; marker++)
				{
					fprintf(indmarkers, " %d %d", ind->markerdata[marker].first.value(), ind->markerdata[marker].second.value());
				}

				fprintf(indmarkers, "\n");
				fprintf(indgrid, "\n");
			}
			/*else
			{
			for (int i = 0; i <= end; i++)
			{
			fprintf(indgrid, " %d", ind->genotypegrid[i]);
			}
			}*/


		}
	}
}

void lockhaplos(individ* ind, int i)
{
	unsigned int j;

	if (ind->lockstart[i] >= chromstarts[i + 1])
	{
		ind->lockstart[i] = 0;
	}
	int bestpos = -1;
	double bestcount = -1;

	for (j = max(chromstarts[i], ind->lockstart[i]); j != chromstarts[i + 1]; j++)
	{
		if (ind->markerdata[j].first == ind->markerdata[j].second) continue;
		if (fabs(ind->haploweight[j] - 0.5) > 0.49999) continue;

		individ* pars[3] = { ind, ind->pars[0], ind->pars[1] };
		double count = 0;

		for (int p = 0; p < 3; p++)
		{
			if (!pars[p] || pars[p]->markerdata.size() <= j) continue;

			for (int i = 0; i < 2; i++)
			{
				if ((&(pars[p]->markerdata[j].first))[i] != UnknownMarkerVal)
				{
					count += 1.0 - (&(pars[p]->markersure[j].first))[i];
				}
			}
			if (p == 0) count *= 2;
		}

		if (count > bestcount)
		{
			bestcount = count;
			bestpos = j;
		}
	}

	if (bestpos == -1)
	{
		return;
	}
	j = bestpos;

	printf("Fixing point: %d %d %d\n", ind->n, i + 1, j);
	fflush(stdout);

	if (j != chromstarts[i + 1])
	{
		ind->haploweight[j] = 0;
		ind->lockstart[i] = j + 1;
	}
}


double dosureval(int what, pair<int, double> val)
{
	if (val.second == 0) return 0;

	/*  double toret = log(val.second);
	toret += log(0.5) * (what * 0.5 - val.first);
	toret *= 2;
	toret /= what;*/
	double toret = log(val.second);
	toret /= what;
	toret *= 4;
	toret = exp(toret);
	toret /= (1 + toret);

	return toret;
}

// Some operations performed when marker data has been read, independent of format.
void postmarkerdata()
{
	int any, anyrem;
	bool latephase = false;


	markerweight.resize(markerposes.size());
	// If inference is active, add "new" marker data until all data has been found.
	if (CORRECTIONINFERENCE) do
	{
#pragma omp parallel for schedule(dynamic,32)
		// all haploweights MUST be non-zero at this point, as we do not explore all shiftflagmode values
		for (int i = 1; i < INDCOUNT; i++)
		{
			individ* ind = getind(i);
			if (ind) ind->children = 0;
		}
#pragma omp parallel for schedule(dynamic,32)
		for (int i = 1; i < INDCOUNT; i++)
		{
			individ* ind = getind(i);
			if (!ind) continue;

			for (int j = 0; j < 2; j++)
			{
				if (ind->pars[j])
				{
					ind->pars[j]->children++;
				}
			}

			if (ind->markerdata.size())
			{
				generation++;

				for (int g = 0; g < ind->markerdata.size(); g++)
				{
					ind->fixparents(g, latephase);
				}
			}
		}

		any = 0;
		anyrem = 0;
		for (int i = 1; i < INDCOUNT; i++)
		{
			individ* ind = getind(i);
			if (!ind) continue;
			// Only run for sex 2, tailored to half sibships
			if (!ind->sex) continue;

			for (int g = 0; g < (int)ind->markervals.size(); g++)
			{
				//				if (g == 1) continue;

				int startsize = ind->markervals[g].size();
				const int known =
					(ind->markerdata[g].first != UnknownMarkerVal) +
					(ind->markerdata[g].second != UnknownMarkerVal);

				const int oldany = any;

				if (latephase && known == 2 && !startsize && ind->gen < 1)
				{
					ind->markerdata[g].second = UnknownMarkerVal;
					ind->markerdata[g].first = UnknownMarkerVal;
					ind->markersure[g] = make_pair(0.0, 0.0);
					any++;
					anyrem++;
				}
				else
					if (known == 2) continue;

				ind->markervals[g].erase(UnknownMarkerVal);
				startsize = ind->markervals[g].size();

				// map insert preserves existing value if present
				if (ind->markerdata[g].first != UnknownMarkerVal) ind->markervals[g].insert(make_pair(ind->markerdata[g].first, make_pair(ind->children, ind->markersure[g].first)));
				if (ind->markerdata[g].second != UnknownMarkerVal) ind->markervals[g].insert(make_pair(ind->markerdata[g].second, make_pair(ind->children, ind->markersure[g].second)));

				if (ind->markervals[g].size() >= 3)
				{
					fprintf(stderr, "Error, too many matches: %d\t%d\n", i, g);
				}
				if (!latephase && ind->markervals[g].size() == 2)
				{
					int knowncount = ind->markervals[g].begin()->second.first + (++ind->markervals[g].begin())->second.first;
					ind->markerdata[g] = make_pair(ind->markervals[g].begin()->first, (++ind->markervals[g].begin())->first);
					ind->markersure[g] = make_pair(dosureval(knowncount, ind->markervals[g].begin()->second), dosureval(knowncount, (++ind->markervals[g].begin())->second));
					any++;
				}
				if (latephase && ind->markervals[g].size() == 1 && startsize == 1 && known == 1 && !ind->pars[0] && !ind->pars[1])
				{
					if (ind->markerdata[g].first == UnknownMarkerVal || ind->markerdata[g].second == UnknownMarkerVal) any++;
					ind->markerdata[g] = make_pair(ind->markervals[g].begin()->first, ind->markervals[g].begin()->first);
					ind->markersure[g] = make_pair(dosureval(ind->children, ind->markervals[g].begin()->second),
						dosureval(ind->children, ind->markervals[g].begin()->second));
				} // DANGEROUS ASSUMPTIONS
				else if (!latephase && ind->markervals[g].size() == 1 && known == 0)
				{
					any++;
					ind->markerdata[g] = make_pair(ind->markervals[g].begin()->first, UnknownMarkerVal);
					ind->markersure[g] = make_pair(dosureval(ind->children, ind->markervals[g].begin()->second), 0.0);
				}


				if (any != oldany) printf("Correction at %d, marker %d (%d;%d) (%lf;%lf)\n", i, g,
					ind->markerdata[g].first.value(), ind->markerdata[g].second.value(),
					ind->markersure[g].first, ind->markersure[g].second);

			}

			for (int g = 0; g < (int)ind->markervals.size(); g++)
			{
				if (ind->markerdata[g].first == sexmarkerval) {
					ind->markerdata[g] = make_pair(ind->markerdata[g].second, ind->markerdata[g].first);
				}
				ind->markervals[g].clear();
			}
		}
		fprintf(stderr, "Number of corrected genotypes: %d\n", any);
		if (latephase)
		{
			latephase = false;
		}
		else
		{
			/*			if (!any && !latephase)
			  {
				  any++;
				  latephase = true;
				  }*/
		}
	} while (any > anyrem);

	for (int i = 1; i < INDCOUNT; i++)
	{
		individ* ind = getind(i);

		// Lock the first position
		if (HAPLOTYPING && ind && ind->haploweight.size())
			// These days we do the locking in all generations
		{
			ind->markervals.clear();
			// Lock the haplotype (in an arbitrary manner) for the first marker in each linkage groupfl
			// Propagation would be more efficient if we locked the mid-position (which should really be determined in cM)
			// Doing so would be much more opaque, though...
			for (unsigned int i = 0; i < chromstarts.size() - 1; i++)
			{
				//if (!ind->pars[0] && !ind->pars[1])
					lockhaplos(ind, i);
			}
		}
	}
}

typedef boost::tuple<individ*, double, int> negshiftcand;
struct inferiorrelated
{
	negshiftcand ourtuple;
	bool anymatch;

	inferiorrelated(negshiftcand init) : ourtuple(init), anymatch(false)
	{
	}

	bool operator() (negshiftcand b)
	{
		if (b.get<0>()->arerelated(ourtuple.get<0>()))
		{
			if (b.get<1>() > ourtuple.get<1>()) return true;
			else
				anymatch = true;
		}

		return false;
	}
};

struct negshifter
{
	int c;
	negshifter(int c) : c(c)
	{
	}

	void operator() (negshiftcand b)
	{
		int minstart = b.get<2>();
		individ* ind = b.get<0>();
		ind->lastinved[c] = minstart;
		for (int p = minstart + 1; p < (int)chromstarts[c + 1]; p++)
		{
			if (p == minstart + 1) fprintf(stdout, "Inv: %d %d\n", ind->n, p);
			ind->haploweight[p] = 1.0f - ind->haploweight[p];
		}
	}
};

bool ignoreflag2(int flag2, int g, int shiftflagmode, int q, int flag2ignore, const flat_map<individ*, int>& relmap, const flat_map<individ*, int>& relmapshift)
{
	int flag2filter = (1 << 30) - 1;
	// Below lines relied on incorrect assumption of remapping of inheritance
	/*  const int selfval = g >> TYPEBITS;
	  if (SELFING && selfval)
		{
	  int basefilter = (HALFNUMPATHS - 1) << 1;
		  if ((((flag2 ^ (g * 2)) / HALFNUMPATHS) ^ (flag2 ^ (g * 2))) & basefilter)
		{
		  return true;
		}
	  flag2filter = basefilter * (selfval == 1 ? HALFNUMPATHS : 1);
	  flag2filter |= 1;
	  //	printf("Selfval is %d, flag2filter is %d, flag2ignore is %d\n", selfval, flag2filter, flag2ignore);
	  }*/
	if (flag2 & (flag2ignore & flag2filter)) return true;

	int marker = -q - 1000;
	for (auto i = as_const(relmap).begin(), j = as_const(relmapshift).begin(); i != relmap.end(); i++, j++)
	{
		int currfilter = (i->second & flag2filter);
		int filtered = ((flag2 ^ (g * 2)) & currfilter);
		// Require ALL bits in the flag to be set, if at least one is set
		if (filtered && filtered != currfilter) return true;
		//if (marker >= 0 && i->first->markerdata[marker].first == UnknownMarkerVal && i->first->markerdata[marker].second == UnknownMarkerVal && (!(flag2 & i->second)))
		if (marker >= 0 && i->first->markerdata[marker].first == i->first->markerdata[marker].second && i->first->markersure[marker].first == i->first->markersure[marker].second && !(((bool)filtered) ^ ((bool)(shiftflagmode & j->second))) &&
			((!RELSKEWSTATES || currfilter != 1 ) && (!SELFING/* || selfgen == 0*/)))
		{
		  //			return false;
		  return true;
		}
	}
	return false;
}


template<int N> struct valuereporter
{
	std::array<double, N> probs;

  valuereporter()
  {
  		for (int k = 0; k < N; k++)
		{
			probs[k] = 0;
		}
	}

	void report(vector<char>& outqueue)
	{
		double probsum = 0;
		for (int i = 0; i < N; i++)
		{
			probsum += probs[i];
		}
		probsum = 1 / probsum;
		for (int i = 0; i < N; i++)
		{
			char string[255];
			int val;
			sprintf(string, "%.5lf%c%n", probs[i] /** probsum*/, i == N - 1 ? '\n' : '\t', &val);
			for (int k = 0; k < val; k++)
			{
				outqueue.push_back(string[k]);
			}
		}
	}
};

struct genotypereporter : valuereporter<3>
{
	void addval(int q, int mapval, int g, int flag2, double val)
	{
		probs[mapval] += val;
	}
};

struct statereporter : valuereporter<NUMTYPES>
{
	void addval(int q, int mapval, int g, int flag2, double val)
	{
		probs[g] += val;
	}
};

template<typename T> T square(T v)
{
	return v * v;
}

void resizecaches()
{
	// Some heaps are not properly synchronized. Putting a critical section here makes the operations not safe,
	// but *safer*.
#pragma omp critical(uglynewhack)
	for (int t = 0; t < NUMSHIFTS; t++)
	{
#if !DOFB 
		factors[t].resize(markerposes.size());
		cacheprobs[t].resize(markerposes.size());
		done[t].resize(markerposes.size());
#else
		fwbwfactors[t].resize(markerposes.size());
		fwbw[t].resize(markerposes.size());
		fwbwdone[t] = 0;
#endif
	}
}

// Global scale factor, 1.0 meaning "use EM estimate".
double scalefactor = 0.002;

pair<int, int> fixtrees(int j)
{
	int flag2ignore = 0;
	int shiftignore = 0;

	reltree.clear();
	relmap.clear();
	relmapshift.clear();
	reltree.push_back(dous[j]);
	relmap[dous[j]] = 1;
	relmapshift[dous[j]] = 1;

	if (HAPLOTYPING)
	{
		flag2ignore = 1;
		shiftignore = 0;
		bool anylev1 = false;
		for (int lev1 = 0; lev1 < 2; lev1++)
		{
			individ* lev1i = dous[j]->pars[lev1];
			if (!lev1i) continue;
			int flag2base = 1 << (1 + lev1 * ((1 << (NUMFLAG2GEN - 1)) - 1));
			int shiftval = (NUMGEN == 3) ? (2 << lev1) : 0;

			if (!lev1i->empty)
			{
				flag2ignore |= flag2base;
				relmap[lev1i] |= flag2base;
				relmapshift[lev1i] |= shiftval;
			}

			bool anypars = false;
			reltree.push_back(lev1i);
			if (NUMGEN > 2)
			{
				for (int lev2 = 0; lev2 < 2; lev2++)
				{
					individ* lev2i = lev1i->pars[lev2];
					if (!lev2i) continue;

					if (!lev2i->empty)
					{
						flag2ignore |= (flag2base << (lev2 + 1));
						relmap[lev2i] |= (flag2base << (lev2 + 1));
						relmapshift[lev2i] |= 0;
						anypars = true;
					}
					reltree.push_back(lev2i);
				}
			}

			if (anypars)
			{
				shiftignore |= shiftval;
			}
			else
			{
				//lev1i->founder = true;
			}
			// Any information of relevance in parents
			if (anypars || !lev1i->empty)
			{
				anylev1 = true;
			}
		}
		if (anylev1)
		{
			shiftignore |= 1;
		}
		else
		{
			dous[j]->founder = true;
		}
		flag2ignore ^= (NUMPATHS - 1);
		shiftignore ^= (NUMSHIFTS - 1);
	}

	sort(reltree.begin(), reltree.end());
	reltree.resize(unique(reltree.begin(), reltree.end()) - reltree.begin());

	return make_pair(shiftignore, flag2ignore);
}

void moveinfprobs(int i, int k, int marker)
{
	// TODO: Sanitize for zero, inf, nan.
	for (int side = 0; side < 2; side++)
	{
		MarkerVal priorval = UnknownMarkerVal;
		std::array<double, 2> compfactors = { 1, 1 };
		//TODO REMOVE PRIOR HANDLING FROM HERE
		/*
		if (reltree[k]->priormarkerdata.size() > marker)
		{
			priorval = (&reltree[k]->priormarkerdata[marker].first)[side];
		}

		if (priorval != UnknownMarkerVal)
		{
			double priorprob = 1.0 - (&reltree[k]->priormarkersure[marker].first)[side];

			MarkerVal nowval = (&reltree[k]->markerdata[marker].first)[side];
			double nowprob = 1.0 - (&reltree[k]->priormarkersure[marker].first)[side];
			if (nowval != priorval)
			{
				nowprob = 1.0 - nowprob;
			}
			//compfactors = { (1.0 - priorprob) * (1.0 - priorprob) / (1.0 - nowprob), priorprob * priorprob / nowprob};
			compfactors = { (1.0 - priorprob), priorprob };
		}*/

		double sum = 0;
		for (auto infval : infprobs[i][side])
		{
			sum += infval.second * compfactors[infval.first == priorval];
		}

		sum = 1 / sum;

		for (auto infval : infprobs[i][side])
		{
		  reltree[k]->infprobs[marker][side][infval.first] += infval.second * compfactors[(int) (infval.first == priorval)] /* * sum*/;
		  if ((reltree[k]->n == 433 && marker >= 4086 && marker <= 4087)) fprintf(stdout, "INFPROBS: %d %d %d %d %lf %lf (%d)\n", reltree[k]->n, marker, side, infval.first, infval.second, sum, shiftflagmode);
		}
		infprobs[i][side].clear();
	}
}

void movehaplos(int i, int k, int marker)
{
	if (haplos[i][0] || haplos[i][1])
	{
		if (fabs(reltree[k]->haploweight[marker] - 0.5) < 0.4999999)
		{
			double b1 = (haplos[i][0] + exp(-400) * maxdiff * maxdiff * 0.5) /*/ reltree[k]->haploweight[marker] /** (1 - reltree[k]->haploweight[marker])*/;// * (1 + 1e-10 - rhfactor);
			double b2 = (haplos[i][1] + exp(-400) * maxdiff * maxdiff * 0.5) /*/ (1 - reltree[k]->haploweight[marker]) /** reltree[k]->haploweight[marker]*/;// * (rhfactor + 1e-10);

			double intended = (b1 - b2) / min(reltree[k]->haploweight[marker], 1 - reltree[k]->haploweight[marker]);
			//intended -= reltree[k]->haploweight[marker];

			bool neg = intended < 0;

			//intended /= sqrt(fabs(intended) + 1.0);
			// if (neg) intended = -intended;

			{
			  reltree[k]->haplobase[marker] += /*log(b1 / b2)*/b1 / (b1 + b2);
				reltree[k]->haplocount[marker] += 1;
			}
		}
		haplos[i][0] = 0;
		haplos[i][1] = 0;
	}
}

void calcdistancecolrowsums(double mwvals[1][1], double rowsums[NUMTYPES], double colsums[NUMTYPES], double &acc3, double &acc4, double mwfvals[1])
{
	for (int g = 0; g < NUMTYPES; g++)
	{
		double acc1 = 0;
		double acc2 = 0;
		for (int g2 = 0; g2 < NUMTYPES; g2++)
		{
			acc1 += mwvals[g][g2];
			rowsums[g] += mwvals[g][g2] * mwvals[g][g2];

			acc2 += mwvals[g2][g];
			colsums[g] += mwvals[g2][g] * mwvals[g2][g];
		}
		acc3 += acc1;
		acc4 += mwfvals[g];

		colsums[g] -= (acc2 * acc2) / NUMTYPES;
		rowsums[g] -= (acc1 * acc1) / NUMTYPES;
	}
}

void updatenegshifts(const bool validg[NUMTURNS], int shifts, int shiftend, int shiftignore, const double rawvals[NUMTURNS][NUMSHIFTS], int j, int marker)
{
	for (int g = 0; g < NUMTURNS; g++)
	{
		if (!validg[g]) continue;

		double val = 0;

		for (int s = shifts; s < shiftend; s++)
		{
			if (s & shiftignore) continue;
			if (rawvals[g][s] < 0) continue;
			val += rawvals[g][s];
		}

		// TODO: Capping val here, but not in sumnegval, leads to some
		// strange results, indeed. Maybe?
		if (!_finite(val) || val < 1e-174) val = 1e-174;
		{
			int g2 = g;
			if (!g) g2 = (1 << 15) - 1;

			// This is hardcoded for the generation count of 3.
			if (NUMGEN == 3)
			{
				dous[j]->negshift[marker] += log(val) * (1.0 - ((g >> 6) & 1) * 2) * ((g2 >> 6) & 1) /*/ sumnegval[6]*/;

				if (dous[j]->pars[0])
					dous[j]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 0) & 1) * 2) * ((g2 >> 0) & 1) /* / sumnegval[0] */;

				if (dous[j]->pars[1])
					dous[j]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 3) & 1) * 2) * ((g2 >> 3) & 1) /* / sumnegval[3] */;

				if (dous[j]->gen >= 2)
				{
					if (dous[j]->pars[0] && dous[j]->pars[0]->pars[0])
						dous[j]->pars[0]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 1) & 1) * 2) * ((g2 >> 1) & 1) /* / sumnegval[1]*/ / dous[j]->pars[0]->children;

					if (dous[j]->pars[0] && dous[j]->pars[0]->pars[1])
						dous[j]->pars[0]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 2) & 1) * 2) * ((g2 >> 2) & 1) /*/ sumnegval[2]*/ / dous[j]->pars[0]->children;

					if (dous[j]->pars[1] && dous[j]->pars[1]->pars[0])
						dous[j]->pars[1]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 4) & 1) * 2) * ((g2 >> 4) & 1) /*/ sumnegval[4]*/ / dous[j]->pars[1]->children;

					if (dous[j]->pars[1] && dous[j]->pars[1]->pars[1])
						dous[j]->pars[1]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 5) & 1) * 2) * ((g2 >> 5) & 1) /*/ sumnegval[5]*/ / dous[j]->pars[1]->children;
				}
			}
#if false
			else
			{
				dous[j]->negshift[marker] += val * (1.0 - ((g >> 2) & 1) * 2) /* * ((g2 >> 2) & 1) */ / (sumnegval[2]);

				if (dous[j]->pars[0])
					dous[j]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 0) & 1) * 2) * ((g2 >> 0) & 1) /*/ (sumnegval[0])*/;

				if (dous[j]->pars[1])
					dous[j]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 1) & 1) * 2) * ((g2 >> 1) & 1) /*/ (sumnegval[1])*/;
				if (dous[j]->pars[0] && dous[j]->pars[1])
				{
					nsm[make_pair(dous[j]->pars[0], dous[j]->pars[1])][marker][g] +=
						log(val) - log(rawvals[0][s]);
				}

				if (g == 1)
					nsm[make_pair(dous[j]->pars[0], dous[j]->pars[0])][marker][2] +=
					log(val) - log(rawvals[0][s]);

				if (g == 4)
					nsm[make_pair(dous[j]->pars[0], dous[j])][marker][2] +=
					log(val) - log(rawvals[0][s]);
			}
#endif
		}
	}
}

#if 0
void oldinfprobslogic(individ * ind, unsigned int j, int iter, int cno, FILE * out)
{
	MarkerVal bestvals[2] = { UnknownMarkerVal, UnknownMarkerVal };
	double bestsure[2] = { ind->markersure[j].first, ind->markersure[j].second };
	bool foundbest = true;
	bool surefound = false;

	flat_map<MarkerVal, double> surenesses;

	if (DOINFPROBS)
	{
		flat_map<MarkerVal, double> sums[2];
		for (int a = 0; a < 2; a++)
		{
			for (auto i = ind->infprobs[j][a].begin(); i != ind->infprobs[j][a].end(); i++)
			{
				sums[a][i->first.second] += i->second;
			}
		}

		sums[0][UnknownMarkerVal] = 1;
		sums[1][UnknownMarkerVal] = 1;

		for (int a = 0; a < 2; a++)
		{
			// Clear old infprobs
			for (int b = 0; b < 2; b++)
			{
				for (int c = 0; c < 2; c++)
				{
					for (int d = 0; d < 2; d++)
					{
						ind->parinfprobs[j][a][b][c][d] = 0;
					}
				}
			}

			double sum = 0;
			MarkerVal bestmarker;
			double bestval = 0;
			double bestval2 = 0;
			flat_map<MarkerVal, double> infprobs;


			for (flat_map<pair<MarkerVal, MarkerVal>, double>::iterator i = ind->infprobs[j][a].begin(); i != ind->infprobs[j][a].end(); i++)
			{
				sum += i->second;
				infprobs[i->first.first] += i->second;
			}
			if (sum <= 1e-12)
			{
				sum = 500;
				if ((&(ind->markerdata[j].first))[a] != UnknownMarkerVal)
				{
					bestval = 1 - (&(ind->markersure[j].first))[a];
					bestval2 = 1 - (&(ind->markersure[j].first))[a];
					//									if (bestval2 < 0.99) bestval2 += 0.01;
					bestmarker = (&(ind->markerdata[j].first))[a];
				}
			}

			for (flat_map<MarkerVal, double>::iterator i = infprobs.begin(); i != infprobs.end(); i++)
			{
				double sureness = i->second / sum;
				double origsureness = sureness;

				if (!((&(ind->markerdata[j].first))[a] == UnknownMarkerVal))
				{
					if (i->first != (&(ind->markerdata[j].first))[a])
					{
						if (sureness > 0.9999)
						{
							sureness = 0.9999;
						}
					}
				}

				if (origsureness > 0.9999) origsureness = 0.9999;
				surenesses[i->first] += 1 - origsureness;

				if (sureness > bestval)
				{
					bestval = sureness;
					bestval2 = origsureness;
					bestmarker = i->first;
				}
			}

			double origsum = sum;
			double sureness = bestval;
			if (sureness > ((iter % 30 == 19 && false) ? 0.99 : 0.49))
			{
				bestvals[a] = bestmarker;
				bestsure[a] = 1 - bestval2;
			}
			else
			{
				printf("Foundbest now false, with sureness %lf, marker %d, pair-half %d for ind %d\n", sureness, j, a, ind->n);
				foundbest = false;
			}
		}
	}

	if (DOINFPROBS)
	{
		if (ind->sex) foundbest = false;

		if (!foundbest && ind->markerdata[j].first != UnknownMarkerVal && ind->markerdata[j].second != UnknownMarkerVal && (bestvals[0] != UnknownMarkerVal || bestvals[1] != UnknownMarkerVal))
		{
			foundbest = true;
			bestvals[0] = UnknownMarkerVal;
			bestvals[1] = UnknownMarkerVal;

			bestsure[0] = 0;
			bestsure[1] = 0;
		}

		if (ind->lastinved[cno] == -1 || true)
		{
			if (foundbest)
			{
				/*bestsure[0] *= 0.99;
				bestsure[1] *= 0.99;*/
				if (bestvals[0] == bestvals[1] && bestvals[0] != UnknownMarkerVal)
				{
					if (bestsure[0] + bestsure[1] > 0.5 && false)
					{
						double bsum = bestsure[0] + bestsure[1];

						int index = bestsure[0] > bestsure[1];
						bestsure[index] /= 2;

						//if (bestsure[index] < 0.01) bestsure[index] = 0.01;

						MarkerVal bestwhat = UnknownMarkerVal;
						double bestsury = 4;

						for (auto i = surenesses.begin(); i != surenesses.end(); i++)
						{
							if (i->first == bestvals[0]) continue;

							if (i->second < bestsury)
							{
								bestsury = i->second;
								bestwhat = i->first;
							}
						}

						if (bestwhat == UnknownMarkerVal) bestsury = 0;
						/*else
						{
						bestsury /= 2;
						bestsury = 1 / (1 + bestsury);
						bestsury *= 2;
						bestsury = (1 - bestsury) / bestsury;
						}*/

						if (bestsury < 0 || bestsury > 1) fprintf(out, "Bestsury problem %d %d %lf\n", ind->n, j, bestsury);

						bestsure[!index] = (2 - bestsury) / 2;
						bestvals[!index] = bestwhat;
					}

					double sum = (bestsure[0] + bestsure[1]) / 2;
					double diff = fabs(bestsure[0] / sum - 1) + fabs(bestsure[1] / sum - 1);
				}
				if (fabs(0.5 - ind->haploweight[j]) == 0.5)
				{
					int min = 0;
					if (bestsure[1] < bestsure[0]) min = 1;
					if (bestvals[0] == bestvals[1] && bestvals[0] != UnknownMarkerVal)
					{
						bestsure[min] = 0;
					}
					else
					{
						bestsure[!min] /*-= bestsure[min]*/ = 0;
						bestsure[min] = 0;
					}
				}
				/*				  if (bestvals[0] == bestvals[1] && bestvals[0] != UnknownMarkerVal)
				{
				double delta = -0.1;

				if (ind->haploweight[j] > 0.5) delta = -delta;

				ind->haploweight[j] += delta;

				if (ind->haploweight[j] < 0.001) ind->haploweight[j] = 0.001;
				if (ind->haploweight[j] > 0.999) ind->haploweight[j] = 0.999;
				}*/

				if ((bestvals[0] == ind->markerdata[j].second || bestvals[1] == ind->markerdata[j].first) && bestvals[0] != bestvals[1])
				{
					/*swap(bestvals[0], bestvals[1]);
					swap(bestsure[0], bestsure[1]);*/
					//ind->haploweight[j] = 1 - ind->haploweight[j];
				}
				bool nochange = true;
				for (int i = 0; i < 2; i++)
				{
					if ((&(ind->markerdata[j].first))[i] != bestvals[i]) nochange = false;
				}

				for (int i = 0; i < 2; i++)
				{
					double old = (&(ind->markersure[j].first))[i];

					if (bestvals[i] == UnknownMarkerVal)
					{
						bestsure[i] = 0;
						continue;
					}



					if ((&(ind->markerdata[j].first))[i] != bestvals[i])
					{
						if (bestvals[i] == bestvals[!i] && (&(ind->markerdata[j].first))[i] == UnknownMarkerVal) old = 0.5;
						else
							old = 1 - old;
					}

					if (old == 0)
					{
						bestsure[i] = 0;
					}

					/*old *= 2;
					bestsure[i] *= old;

					if (bestsure[i] < old * 0.35) bestsure[i] = old * 0.35;*/

					//if (bestsure[i] <= old * 1.021) bestsure[i] *= 0.98;
					//								double factor = fabs(ind->haploweight[j] - 0.5) + 0.01;
					//								bestsure[i] *= factor;
					bestsure[i] += old/* * (0.51 - factor)*/;
					bestsure[i] /= 2 /* 0.51 */;


					/*bestsure[i] *= bestsure[i];
					bestsure[i] *= 2;*/
				}

				/*							if (nochange && bestvals[0] != bestvals[1])
				{
				bestsure[0] *= 0.99;
				bestsure[1] *= 0.99;
				}*/

				if (nochange && bestvals[0] == bestvals[1] && bestvals[0] != UnknownMarkerVal)
				{
					double ratio = bestsure[0] / bestsure[1];
					if (ratio < 1) ratio = 1 / ratio;

					if (ratio < 1.0001)
					{
						double val = bestsure[0] * 0.01;
						bestsure[0] -= val;
						bestsure[1] += val;
						if (bestsure[0] > 0.01)
						{
							fprintf(out, "Bringing aside %d %d %lf\n", ind->n, j, bestsure[0]);
						}
					}
				}

				/*							if (bestvals[0] != ind->markerdata[j].first || bestvals[1] != ind->markerdata[j].second || bestsure[0] + bestsure[1] + ind->markersure[j].first + ind->markersure[j].second > 0.01 || bestvals[0] == UnknownMarkerVal || ind->n == 1633)
				{
				fprintf(out, "Individual %d stochastic fix at marker %d, %d:%d, was %d:%d %c %lf %lf\n", ind->n, j, bestvals[0], bestvals[1], ind->markerdata[j].first, ind->markerdata[j].second, surefound ? 'S' : 'I', bestsure[0], bestsure[1]);
				}*/

				ind->markerdata[j] = make_pair(bestvals[0], bestvals[1]);
				ind->markersure[j] = make_pair(bestsure[0], bestsure[1]);
			}
		}
		for (int a = 0; a < 2; a++)
		{
			ind->sureinfprobs[j][a].clear();
			ind->infprobs[j][a].clear();
			ind->unknowninfprobs[j][a] = 0;
		}
	}
}
#endif

double caplogitchange(double intended, double orig, double epsilon, bool& hitnnn)
{
	double nnn = 100;
	if (nnn < 1.0) nnn = 1.0;

	double limn = (nnn - 1.0) * orig * (-1 + orig);

	double limd1 = -1 - (nnn - 1.0) * orig;
	double limd2 = (nnn - 1.0) * orig - nnn;

	double lim = min(limn / limd1, limn / limd2);

	intended = min((double)intended, 1.0 - epsilon);
	intended = max((double)intended, epsilon);
	double diff = intended - orig;


	if (diff > limn / limd1)
	{
	  if (!hitnnn)  fprintf(stderr, "CAP: Exceeding limit %lf > %lf, intended %lf, orig %lf\n", diff, limn/limd1, intended, orig);
		intended = orig + limn / limd1;
		hitnnn = true;
	}

	if (diff < -limn / limd2)
	{
	  if (!hitnnn) fprintf(stderr, "CAP: Underflowing limit %lf < %lf, intended %lf, orig %lf\n", diff, -limn/limd2, intended, orig);
		intended = orig - limn / limd2;
		hitnnn = true;
	}

	return intended;
}

// The actual walking over all chromosomes for all individuals in "dous"
// If "full" is set to false, we assume that haplotype inference should be done, over marker positions.
// A full scan is thus not the iteration that takes the most time, but the scan that goes over the full genome grid, not only
// marker positions.
template<bool full, typename reporterclass> void doit(FILE* out, bool printalot
#ifdef F2MPI
	, mpi::communicator& world
#endif
	)
{
	const bool doprint = full;
#ifdef F2MPI
	broadcast(world, actrec, 2, 0);
#endif

	int count = 0;
	vector<vector<std::array<float, 2> > > realgeno;

	realgeno.resize(dous.size());

	for (int a = 1; a <= 20; a++)
	{
		individ* ind2 = getind(a);
		/*	    ind2->semishift.resize(5000);
		for (int q = 0; q < 5000; q++)
		{
		for (int z = 0; z < 40; z++)
		{
		for (int p = 0; p < 2; p++)
		{
		for (int r = 0; r < 2; r++)
		{
		ind2->semishift[q][p * 2 + r][z] = (z == (a - 1) * 2 + p);
		}
		}
		}
		*/
	}

	for (unsigned int j = 0; j < dous.size(); j++)
	{
		if (dous[j]->markerdata.size())
		{
			count++;
		}
		else
		{
			dous.erase(dous.begin() + j);
		}
	}

	static int iter = 0;
	iter++;


	map<pair<individ*, individ*>, map<int, std::array<double, 8> > > nsm;

	for (int i = 0; i < INDCOUNT; i++)
	{
		individ* ind = getind(i);
		if (!ind) continue;
		ind->children = 0;
		ind->kids.clear();
		if (!ind || !ind->haplocount.size()) continue;

		for (unsigned int j = 0; j < ind->haplocount.size(); j++)
		{
			ind->haplocount[j] = 0.0f;
			ind->haplobase[j] = 0.0f;
		}

#ifdef F2MPI
		broadcast(world, ind->haploweight, 0);

		//		fprintf(out, "C%d:%d\n", world.rank(), i);

		world.barrier();
#endif

	}

	for (int j = 0; j < (int)dous.size(); j++)
	{

		for (int i = 0; i < 2; i++)
		{
			individ* pnow = dous[j]->pars[i];

			if (pnow)
			{
				pnow->children++;
				pnow->kids.push_back(dous[j]);
			}
		}
	}

	vector<set<negshiftcand> > negshiftcands;
	negshiftcands.resize(chromstarts.size());
	// Create a vector where each element corresponds to a marker and
	//contains a referense to a vector containing all the clauses for said marker
	//EBBA also: Here starts the parallels, investigate names
	vector<vector<clause>> toulInput;
	toulInput.resize(markerposes.size()); //EBBA

	for (unsigned int i = 0; i < chromstarts.size() - 1; i++)
	{
		//printf("Chromosome %d\n", i + 1);

		// The output for all individuals in a specific iteration is stored, as we have parallelized the logic and 
		// want the output to be done in order.
		vector<vector<char> > outqueue;
		outqueue.resize(dous.size());

		std::set<int> indnumbers;// to count individuals
		long long maxweight = 0;

#pragma omp parallel for schedule(dynamic,1)
		for (int j = 0; j < (int)dous.size(); j++)
		{
#ifdef F2MPI
			if (j % world.size() != world.rank()) continue;
#endif
			//			realgeno[j].resize(markerposes[chromstarts[1] - 1] + 1);

			generation++;
			threadblock tborig;
			threadblock tb = tborig;

			resizecaches();

			if (dous[j]->markerdata.size())
			{
				int qstart = -1000 - chromstarts[i];
				int qend = -1000 - chromstarts[i + 1];
												 //TODO: WRONG RANGE
				int qd = -1;
				int f2s = 0;
				int f2end = NUMPATHS;

				if (!HAPLOTYPING)
				{
					f2s = -1;
					f2end = 0;
				}

				int shifts = 0;
				int shiftend = NUMSHIFTS;

				// Special optimization hardcoded for this population structure, eagerly skipping flags that do not
				// correspond to any inheritance, i.e. if not the full pedigree of 6 individuals back is present.
				int shiftignore;
				int flag2ignore;

				std::tie(shiftignore, flag2ignore) = fixtrees(j);

				bool skipsome = false;
				for (int u = 0; u < reltree.size() && !skipsome; u++)
				{
					for (int p = 0; p < 2; p++)
					{
						if ((&(reltree[u]->markerdata[chromstarts[i]].first))[p] == sexmarkerval)
							skipsome = true;
					}
				}

				if (!HAPLOTYPING)
				{
					reltree.resize(0);
					flag2ignore = 0;
				}

				if (full)
				{
					qstart = (int)markerposes[chromstarts[i]];
					qend = (int)markerposes[chromstarts[i + 1] - 1] + 1;
					qd = 1;
					f2s = -1;
					f2end = 0;
					/*shifts =s 0;
					shiftend = 1;*/
				}

				if (dous[j]->gen < 2) shiftend = min(2, shiftend);

				// BIG SHIFT HACK, removed again as it makes resolving a shared 12 12 ancestry very hard
				bool anygood = false;
				anygood |= dous[j]->pars[0] && !dous[j]->pars[0]->empty;
				anygood |= dous[j]->pars[1] && !dous[j]->pars[1]->empty;

				//if (!anygood) shiftend = 2;
				/*if (!anygood)
				{
					shiftignore = 6;
					flag2ignore = 0;
					}*/

				double factor = -1e15;
				double factors[NUMSHIFTS];
				for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
				{
					if (shiftflagmode & shiftignore) //continue;
						factors[shiftflagmode] = -1e30;
					else
						factors[shiftflagmode] = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, NONESTOP, -1, false, 0, -10000 + factor);
					factor = max(factor, factors[shiftflagmode]);
				}

				// Normalize!
				double realfactor = 0;
				for (int s = shifts; s < shiftend; s++)
				{
					if (s & shiftignore) continue;
					realfactor += exp(factors[s] - factor);
				}
				int countf2i = 0;
				for (int k = 0; k < TYPEBITS + 1; k++)
				{
					if (flag2ignore & (1 << k)) countf2i++;
				}
				//factor -= log(1 << countf2i);

				// This output can get ugly due to race conditions. One shouldn't rely on it.
				printf("%d,%03d,%03d: %lf\t", dous[j]->n, flag2ignore, shiftignore, factor);
				factor += log(realfactor);
				printf("%lf %d\n", factor, shiftend);
				fflush(stdout);
				if (_isnan(factor) || factor < MINFACTOR) continue;

				// Walk over all chromosome positions, whether it be markers (negative q values <= -1000) or grid positions
				for (int q = qstart; q != qend; q += qd)
				{
					reporterclass reporter;
					//double mwvals[NUMTYPES][NUMTYPES] = {0};
					//double mwfvals[NUMTYPES] = {0};
					double mwvals[1][1];
					double mwfvals[1];
					double mwval[4] = { 0 };

					for (int g = 0; g < NUMTYPES; g++)
					{
						for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
						{
							if (shiftflagmode & shiftignore) continue;
							if (factor - factors[shiftflagmode] > 40) continue;
							if (q <= -1000 && DOREMAPDISTANCES)
							{
								double val;
								for (int g2 = 0; g2 < NUMTYPES; g2++)
								{
									val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, twicestop(q, g, g2),
										-1, true, 0, -5000.0 + factor) - factor;

									val = exp(val);
									mwvals[g][g2] += val;
								}

								val = 0;
								/*								if (g == 0)
								{
								stopmodpair::miniactrecT newactrec;
								for (int gs = 0; gs < 2; gs++)
								{
								for (int k = 0; k < 2; k++)
								{
								newactrec[k][gs] = actrec[k][-q - 1000 + 1 + gs];
								}
								}

								double sum = newactrec[0][0] + newactrec[0][1];
								sum /= 2;
								newactrec[0][0] = 0;
								newactrec[0][1] = sum * 2;

								val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, stopmodpair(q, newactrec),
								-1, true, 0, -20.0 + factor) - factor;
								}*/


								mwfvals[g] += val;
							}

							for (int flag2 = f2s; flag2 < f2end; flag2++)
							{
								if (ignoreflag2(flag2, g, shiftflagmode, q, flag2ignore, relmap, relmapshift)) continue;
								//if (flag2 & (flag2ignore)) continue;

								int firstpar = 0;
								double val;

								if (DOIMPOSSIBLE && q <= -1000 && flag2 != -1)
								{
									// Do a lookup in the impossible structure. If we are at a marker position, those branches that are impossible will have an exact
									// score of 0.
									// If we are NOT at a marker, this information can still be used in trackpossible, but not in here, as there is some recombination
									// added between the marker and fixated position we're checking.
									for (int a = 0; a < 2; a++)
									{
										if (a) firstpar = !firstpar;
										const int genwidth = (1 << (NUMGEN - 1));

										int f2n = ((flag2 ^ shiftflagmode) & 1);

										int upflagr = upflagit(g, firstpar, genwidth);
										int upflag2r = upflagit(flag2 >> 1, firstpar, genwidth);
										int upshiftr = upflagit(shiftflagmode >> 1, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));
										int marker = -(q + 1000);
										int impossibleval = generation * markerposes.size() + marker;

										if (impossible[shiftflagmode & 1][firstpar][f2n][upflagr][upflag2r + 1][upshiftr][marker & 3] == impossibleval)
										{
											std::cerr << "IMPOSSIBLE FACILITY USED" << std::endl;
											goto continueloop;
										}
									}
								}


								// This is the big main call
								val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, classicstop(q, g),
									flag2, true, 0, -40.0 + factor) - factor;

								if (_finite(val) && val > -40.0)
								{
									// shift mode not included, this is the "real" f2n, indicating what value
									// in the marker pair is used, not the strand phase (strand phase is flag2 xored
									// with the other stuff)
									int f2n = ((flag2 /*^ shiftflagmode*/)& 1);


									val = exp(val);
									int marker = -q - 1000;

									int mapval = 0;
									double outmapval = dous[j]->trackpossible<false, true>(tb, UnknownMarkerVal, 0, marker, g * 2, flag2, *(tb.shiftflagmode), trackpossibleparams(0, &mapval));
									double sidevals[2][2] = { 0 };
									double sidevalsums[2] = { 0 };

									if (DOINFPROBS)
									for (int side = 0; side < 2; side++)
									{
										for (auto markerval : { 1 * MarkerValue, 2 * MarkerValue })
										{
											double sideval = dous[j]->trackpossible<false, false>(tb, markerval, 0, marker, g * 2 + side, flag2 ^ side, *(tb.shiftflagmode), trackpossibleparams(0, nullptr));
											sidevals[side][markerval.value() - 1] += sideval;
											sidevalsums[side] += sideval;
										}
									}
//									double pairvalsum = accumulate(begin(pairvals), end(pairvals), 0.0);

/*									for (int i = 0; &pairvals[i] != end(pairvals); i++)
									{
									  reporter.addval(q, i, g, flag2, val * pairvals[i] / pairvalsum);
									}*/
									/*/*if (!outmapval && val > 1e-3)
									{
										//std::cerr << "ERROR TERROR " << -q - 1000 << " " << g * 2 << " " << flag2 << " " << *(tb.shiftflagmode) << "\t" << val << "\n";
									}
									if (mapval < 0 || mapval > 2)
									  {
									  std::cerr << "Incorrect mapval " << -q - 1000 << " " << dous[j]->n << std::endl;
									  }*/
									//reporter.addval(q, mapval, g, flag2, val);
									if (!full && HAPLOTYPING)
									{
										dous[j]->updatehaplo(tb, -q - 1000, g, flag2, val);

										if (DOINFPROBS)
										for (int side = 0; side < 2; side++)
										{
											for (auto markerval : { 1 * MarkerValue, 2 * MarkerValue })
											{
												//std::cout << "EXTREME VETTING IND " << dous[j]->n << " MARKER " << marker << ":" << markerval.value() << ", fl2" << flag2 << ", sfm " << *(tb.shiftflagmode) << ", VAL: " << val << " SIDEVAL " << sidevals[side][markerval.value() - 1] << ", SIDEVALSUM " << sidevals[side][markerval.value() - 1] << std::endl;
												double updateval = val * sidevals[side][markerval.value() - 1] / sidevalsums[side];
												dous[j]->trackpossible<GENOS, false>(tb, markerval, 0, marker, g * 2 + side, flag2 ^ side, *(tb.shiftflagmode), trackpossibleparams(updateval, nullptr));
											}
										}
									}
								}
							continueloop:;
							}
						}
					}

					// Coordinate estimated distances from all individuals.
					if (q <= -1000 && DOREMAPDISTANCES)
					{
						double colsums[NUMTYPES] = { 0 };
						double rowsums[NUMTYPES] = { 0 };
						double acc3 = 0;
						double acc4 = 0;
						calcdistancecolrowsums(mwvals, rowsums, colsums, acc3, acc4, mwfvals);

						double relinfo1 = 0;
						double relinfo2 = 0;
						double infosum[2] = { 0 };
						double recprob[2];
						double dist = markerposes[-q - 1000 + 1] - markerposes[-q - 1000];

						for (int k = 0; k < 2; k++)
						{
							recprob[k] = 0.5 * (1.0 - exp(actrec[k][-q - 1000 + 1] * (dist)));
							recprob[k] = max(1e-8, recprob[k]);
						}

						for (int g = 0; g < NUMTYPES; g++)
						{
							relinfo1 += colsums[g];
							relinfo2 += rowsums[g];
							for (int g2 = 0; g2 < NUMTYPES; g2++)
							{
								double infofactor = /*rowsums[g] * colsums[g2]*/ 1;
								double val = mwvals[g][g2] * infofactor;

								int mask = g ^ g2;
								for (int i = 0; i < TYPEBITS; i++)
								{
									double corr = 1;
									if (NUMGEN > 2)
									{
										if (i % (TYPEBITS >> 1))
										{
											// TODO: handle null
											corr /= dous[j]->pars[i / (TYPEBITS) >> 1]->children;
										}
									}
									bool switched = ((bool)(mask & (1 << i)));
									double expected = !switched ? 1.0 - recprob[TYPESEXES[i]] : recprob[TYPESEXES[i]];

									mwval[TYPESEXES[i] * 2 + switched] += (val / expected) * corr;
									infosum[TYPESEXES[i]] += corr;
								}
							}
						}

						double summw = (mwval[0] + mwval[1] + mwval[2] + mwval[3]);
						for (int z = 0; z < 2; z++)
						{
							infosum[z] /= NUMTYPES * NUMTYPES;
						}
						if (acc3 == 0)
						{
							acc3 = 1;
							infosum[0] = 0;
							infosum[1] = 0;
						}


						double delta = 0;
#pragma omp critical(markerweights)					       
						{
							for (int t = 0; t < 4; t++)
							{
								int tbase = (t >> 1) << 1;

								double dval = (mwval[t] / acc3 - infosum[t / 2]) /*/ summw*/
									* fabs(((t & 1) ? 0.0 : 1.0) - recprob[t / 2])
									/** (-1 + (t & 1) * 2)*/;

								markerweight[-q - 1000][t] += dval;
								if (t == 0) delta = dval;
							}
							markerweight[-q - 1000][4] = min(markerweight[-q - 1000][4], acc4);
						}
						if (delta < -0.1 && q > qend + 2) fprintf(out, "RECOMB:\t%d\t%d\t%lf\n", dous[j]->n, q, delta);
					}

					// TODO: NEGSHIFT DOESN'T TAKE RELMAP FLAG2 RESTRICTIONS INTO ACCOUNT
					// Consider doing haplotype reversal from a specific position and all the way down.
					if (HAPLOTYPING && !early && !full && dous[j]->gen >= 0)
					{
						int marker = -q - 1000;

						if (RELSKEWS && !RELSKEWSTATES && false)
#pragma omp critical(negshifts)
						{
							double prevval = dous[j]->haploweight[marker];

							// TODO: OVERRUN AT MARKER + 1 ?
							for (int k = 0; k < 2; k++)
							{
								double relval = fabs(k - dous[j]->relhaplo[marker]);
								double sum = 0;
								double term = dous[j]->haploweight[marker + 1] * (prevval * relval + (1 - prevval) * (1 - relval));
								double lo = term;
								sum += term;

								term = (1 - dous[j]->haploweight[marker + 1]) * ((1 - prevval) * relval + prevval * (1 - relval));
								sum += term;
								dous[j]->negshift[marker] += (0 == k ? 1 : -1) * log(sum);
							}
						}
						
						double rawvals[NUMTURNS][NUMSHIFTS];
						double rawervals[NUMTURNS][NUMSHIFTS];
						double sumnegval[TYPEBITS + 1] = { 0 };
						bool validg[NUMTURNS] = { 0 };

						for (int g = 0; g < NUMTURNS; g++)
						{
							for (int s = 0; s < NUMSHIFTS; s++)
							{
								rawvals[g][s] = -1;
								rawervals[g][s] = -1;
							}
						}

						for (int g = 0; g < NUMTURNS; g++)
						{
							if (g & (flag2ignore >> 1)) continue;

							aroundturner turn(g);
							for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
							{
								if (shiftflagmode & shiftignore) continue;
								// If we are above this limit, we are shifting shift mode
								// within the range and cannot use this heuristic of the
								// aggregated probability to know anything... anything at all!
								// (In other cases, our simulated crossover event amounts to 
								// far below a exp(20) change in probability.)
								int g2 = g;
								if (!g) g2 = (1 << 15) - 1;

								int oldshift = shiftflagmode;
								rawervals[g][oldshift] = exp(dous[j]->doanalyze<aroundturner>(tb, turn, chromstarts[i],
									chromstarts[i + 1] - 1, classicstop(q, -1), -1, true, 0, -5000 + factor) - factor);
								shiftflagmode = oldshift;

								/*								if (c > 1) continue;*/
								validg[g] = true;

								if (!isfinite(rawervals[g][oldshift]))
								{
									rawvals[g][oldshift] = 0;
									continue;
								}

								rawvals[g][oldshift] = rawervals[g][oldshift];

								for (int t = 0; t < TYPEBITS + 1; t++)
								{
									if (g2 & (1 << t))
									{
										sumnegval[t] += rawvals[g][shiftflagmode];
									}
								}

							}

						}

						for (int t = 0; t < TYPEBITS + 1; t++)
						{
							if (!sumnegval[t]) sumnegval[t] = 1;
						}

						// Changed for EBBA's degree project
						// Reminder: NUMTURNS is a bitflag with relevant individuals that should be changed
						// remaining individuals should not be changed (= negative numbers)
						// All individuals ancestors (assuming max 3 gens) should be included in clause
						// The weight is the sum of rawvals (?)
						// The current individual is dous[j]
						// structure of container: vector<vector<clause*>> toulInput;

						{
							//Markers stored in q (see line 2844 in original cnF2freq)
							//such that
							int mark = -q - 1000;
							int numbind = 1;

							//Do we have parents and grand parents?
							//Store their identifying numbers in an array.
							//Hard coded for max 3 gens.
							//std::fstream test("test2.txt", ios::out | ios::in | ios::trunc);//TEST//TEST
							vector<int> cands(7);
							vector<bool> exists(7, false);
							std::set<int> family;
							int temp = dous[j]->n;
							cands[6] = temp;
							exists[6] = true;
							family.insert(temp);



							// If incest, we preted the person did not sire anyone the second time they show up in the focus tree

							//test << "Mark: " << mark << " Individ: " << dous[j]->n;

							if (dous[j]->pars[0]) {
								temp = dous[j]->pars[0]->n;
								if (family.insert(temp).second) {//if family member is unique
									cands[0] = temp;
									exists[0] = true;
									numbind++;
								}
								//test << " Parent1: " << temp;

								if (dous[j]->pars[0]->pars[0]) {
									temp = dous[j]->pars[0]->pars[0]->n;
									if (family.insert(temp).second) {
										cands[1] = temp;
										exists[1] = true;
										numbind++;
									}
									//test << " Parent1's parent1: " << dous[j]->pars[0]->pars[0]->n;
								}
								if (dous[j]->pars[0]->pars[1]) {
									temp = dous[j]->pars[0]->pars[1]->n;
									if (family.insert(temp).second) {
										cands[2] = temp;
										exists[2] = true;
										numbind++;
									}
									//test << " Parent1's parent2: " << dous[j]->pars[0]->pars[1]->n;
								}
							}

							if (dous[j]->pars[1]) {
								temp = dous[j]->pars[1]->n;
								numbind++;
								if (family.insert(temp).second) {
									cands[3] = temp;
									exists[3] = true;
									numbind++;
								}
								//test << " Parent2: " << dous[j]->pars[1]->n;
								if (dous[j]->pars[1]->pars[0]) {
									temp = dous[j]->pars[1]->pars[0]->n;
									if (family.insert(temp).second) {
										cands[4] = temp;
										exists[4] = true;
										numbind++;
									}
									//test << " Parent2: " << dous[j]->pars[1]->pars[0]->n;
								}
								if (dous[j]->pars[1]->pars[1]) {
									temp = dous[j]->pars[1]->pars[1]->n;
									if (family.insert(temp).second) {
										cands[5] = temp;
										exists[5] = true;
										numbind++;
									}
									//test << " Parent2: " << dous[j]->pars[1]->pars[1]->n;
								}
							}
							//test << "Number of individuals:  " << numbind << " End of input into cands \n ";//remember, incesters only counted once

																											//Use structure clause to store weight and values.

																											//loop over NUMTURNS
																											// get relevant weights and individuals, insert to container

																											//for (int g = 0; g < NUMTURNS; g++) {
							double normsum = 0;
							for (int s = shifts; s < shiftend; s++) {
								if (s & shiftignore || rawvals[0][s] <= 0) continue;
								normsum += rawvals[0][s];
							}
							double normfactor = 1 / normsum;

							for (int g = 0; g < NUMTURNS; g++) {
								if (g & (flag2ignore >> 1)) continue;

								std::bitset<16> bits(g);
								vector<int> claus;
								int cind;
								int shiftcount = 0;
								for (int b = 0; b < 7; b++) {
									if (exists[b]) {
										cind = cands[b];
										//if (find(claus.begin(), claus.end(), cind | -cind) == claus.end()){// avoid inbreeding results.
										if (bits[b]) {
											claus.push_back(-cands[b]);
											shiftcount++;
										}
										else {
											claus.push_back(cands[b]);
										}
										//}
									}
								}
								double w = 1e-100;

								for (int s = shifts; s < shiftend; s++) {
									if (s & shiftignore || rawvals[g][s] <= 0) continue;
									w += rawvals[g][s];
								}
								w *= normfactor;
								//Now simply construct a clause type and send it to the right marker
								clause c;
								w = log(w) * 1000000000;
								c.weight = w;
								c.weight -= g;
								c.cinds = claus;
								//test << "Mark: " << mark << "ClausToString: " << c.toString() << " Current maxweight: " << maxweight << endl;//TEST
#pragma omp critical(negshifts)
								{
								  if (c.weight > maxweight) {
								    maxweight = c.weight;
								  }
								  toulInput[mark].push_back(c);
								}
							}
#pragma omp critical(negshifts)
							for (int b = 0; b < 7; b++) {
								if (exists[b]) {
									indnumbers.insert(cands[b]);
								}
							}
						}
					}

					reporter.report(outqueue[j]);

					if (!full)
					{
						int marker = -q - 1000;

						// Contribute haplotype data, but truncate it, i.e. a 50/50 contribution for either interpretation is not added.
						// Instead, we have a cap later on at the maximum change at any iteration.
						// critical section outside the loop to maintain symmetry in floating point ops
#pragma omp critical(update)
						{
#pragma ivdep
							for (int k = 0; k < (int)reltree.size(); k++)
							{
								int i = reltree[k]->n;
								moveinfprobs(i, k, marker);
								movehaplos(i, k, marker);
							}
						}
					}
				}

			}
			fprintf(stderr, "LAST: %d %s\n", dous[j]->n, dous[j]->name.c_str());
			fflush(stderr);
		}

		//Print all information to seperate files EBBA
		// stored by marker in toulIn  vector<vector<clause>>
		//Then run toulbar and save best solution in relevant negshift vectors
		//Remember: typedef boost::tuple<individ*, double, int> negshiftcand;

		int nbvar = indnumbers.size();
		long long minsumweight = std::numeric_limits<long long>::max();
		std::set<negshiftcand> bestcands;
		//for (int m=0; m < (int) toulInput.size(); m++ ){//TODO change so that it is valid for more than one chromosome
#pragma omp parallel for schedule(dynamic,1)
		for (int m = chromstarts[i]; m < chromstarts[i + 1]; m++) {
		  if (m % 10) continue;
			std::string tid = boost::lexical_cast<std::string>(omp_get_thread_num());
			std::string toulin(std::string("toul_in") + tid + ".wcnf");
			std::string toulout(std::string("toul_out") + tid + ".txt");
			std::string sol(std::string("sol") + tid);
			std::fstream infile(toulin, ios::out | ios::in | ios::trunc);			
			if (!infile) {
				perror("Toulbars input file failed to open to be written to because: ");
			}

			infile << "c In Weigthed Partial Max-SAT, the parameters line is 'p wcnf nbvar nbclauses top'\n";
			infile << "c Where p is the weight\n";
			infile << "c nbvar is the nu	mber of a variables appearing in the file (TYPEBITS +1)\n";
			infile << "c nbclauses is the exact number of clauses contained in the file\n";
			infile << "c see http://maxsat.ia.udl.cat/requirements/\n";

			int nbclauses = (int)toulInput[m].size();
			int nbc = nbclauses + nbvar * 2;
			//cout<<"nbvar: " <<nbvar<< "\n"; // problem solving
			//cout<<"nbclauses: " <<nbc<< "\n"; // problem solving
			infile << "p wcnf " << 999 << " " << nbc << "\n"; //" " <<std::numeric_limits<int>::max()<<"\n";

			for (auto cind : indnumbers) {//add clauses to get output variables sorted by size.
				infile << "1 " << cind << " 0\n";
				infile << "1 " << -cind << " 0\n";
			}

			for (int g = 0; g < nbclauses; g++) {
				clause& c = toulInput[m][g];
				c.weight = maxweight - c.weight + 1;
				if (c.weight < 0)
				  {
				    fprintf(stderr, "Negative weight marker %d, clause %d, weight %lld, maxweight %lld\n", m, g, c.weight, maxweight);
				  }
				infile << c.weighttostring() << c.clausetostring() << " 0\n";
				//infile<< toulInput[m][g].toString() << "\n";
				//cout<<"TEST " <<toulInput[m][g].toString()<< "\n"; // problem solving
			}
			infile.close();


			string str = "toulbar2 " + toulin + " -p=8 -m=1 -w=" + sol + " -s > " + toulout; //works as in it runs, not as in it actually does what we want
																			 //string str = "toulbar2 brock200_4.clq.wcnf -m=1 -w -s";//TEST


																			 // Convert string to const char * as system requires
			const char *command = str.c_str();
			system(command);

			//Read outfile and store best result in negshift
			std::fstream touloutput(sol, ios::in);
			//read from file to string of 0s and 1s
			int rawinput;
			vector<int> tf;
			while (touloutput >> rawinput) {
				tf.push_back(rawinput);
				if (rawinput)
				  {
				    fprintf(stderr, "Ind %d marker %d shift indicated\n", tf.size() - 1, m);
				  }
			}

			//Identify all violated clauses, elimination step means optimum cost data from toulbar not usable.
			long long sumweight = 0;
			int vc = 0;
			for (int g = 0; g < nbclauses; g++) {
				bool viol = true;
				for (int val : toulInput[m][g].cinds)
				{
					int ind = val < 0 ? -val : val;
					if (tf[ind - 1] == (val > 0))
					{
						viol = false;
					}
				}
				if (viol)
				  {
				    sumweight += toulInput[m][g].weight;
				    vc++;
				  }
			}

			fprintf(stderr, "Marker %d score %lld, %d/%d clauses, %d vars in tf\n", m, sumweight, vc, nbclauses, tf.size());
			
			

#pragma omp critical(negshifts)
			if (minsumweight > sumweight) {
				//vector containing all individuals numbers who should be shifted
				minsumweight = sumweight;
				//vector<int> neg;
				bestcands.clear();
				for (int g = 0; g<tf.size(); g++) {
					if (tf[g]) {
					  bestcands.emplace(getind(g + 1), sumweight, m);
						//neg.push_back(inds[g]);
					}
				}

				//if(neg.size()>1){
				//cout<< "There is a place where double shifts would be good!"<< endl;//(string) neg <<
				//}
			}


			//Close file
			touloutput.close();
			//output.close();
		}

		//Data structure to fill: vector<set<negshiftcand> > negshiftcands (0);
		if (bestcands.size() == 1) {
			cout << "Only one individual is switching!" << endl;
		}
		else if (bestcands.size()>1) {
			cout << "A switch!" << bestcands.size() << "individuals are changing!" << endl; //For Ebbas degree projects results
		}
		negshiftcands[i] = bestcands;
		bestcands.clear(); // uneccesary, just for clarity
						   //End of Ebbas code

		for (unsigned int j = 0; j < outqueue.size(); j++)
		{
			outqueue[j].push_back(0);
			if (out && printalot) fprintf(out, "%s:%d\n", dous[j]->name.c_str(), i + 1);
			if (out && printalot) fprintf(out, "%s\n", &outqueue[j].front());
		}

		if (out) fflush(out);

#ifdef F2MPI
		// TODO: Shouldn't this be turned on if F2MPI is enabled?
		if (false) reduce(world, markerweight, markerweight2, vectorplus<MWTYPE>(), 0);
#endif
		if (DOREMAPDISTANCES
#ifdef F2MPI
			&& world.rank() == 0
#endif
			&& !full)
		{
			markerweight = markerweight2;
			for (int q = chromstarts[i]; q < chromstarts[i + 1] - 2; q++)
			{
				double dist = markerposes[q + 1] - markerposes[q];
				if (dist > 0)
				{
					for (int t = 0; t < 2; t++)
					{
						double prob2 = -(markerweight[q][t * 2 + 1] - markerweight[q][t * 2]) / dous.size() / dist;
						if (prob2 > 0.3) prob2 = 0.3;
						if (prob2 < -0.3) prob2 = -0.3;

						actrec[t][q + 1] += prob2;
						if (actrec[t][q + 1] > -1e-5) actrec[t][q + 1] = -1e-4;
						if (actrec[t][q + 1] < -20) actrec[t][q + 1] = -10;

						if (q != chromstarts[i] && markerweight[q - 1][4] < dous.size() * 1e-3 && false)
						{
							double sum = actrec[t][q + 1] + actrec[t][q];
							sum /= 2;

							actrec[t][q] = sum;
							actrec[t][q + 1] = sum;
						}
					}
				}
				fprintf(out, "actrec: %0.3lf\t%0.3lf\n", actrec[0][q + 1] * 50, actrec[1][q + 1] * 50);
			}
		}


		if (!full && HAPLOTYPING)
		{			
			for (unsigned int i = 0; i < INDCOUNT; i++)
			{
				individ* ind = getind(i);
				if (!ind || !ind->haplocount.size()) continue;

				{
#ifdef F2MPI
					vector<float> templist;
					templist.resize(ind->haplocount.size());

					reduce(world, ind->haplocount, templist, vectorplus<vector<float> >(), 0);
					ind->haplocount = templist;
					world.barrier();

					reduce(world, ind->haplobase, templist, vectorplus<vector<float > >(), 0);
					ind->haplobase = templist;
					world.barrier();

					reduce(world, ind->negshift, templist, vectorplus<vector<float > >(), 0);
					ind->negshift = templist;
					world.barrier();
#endif


							fprintf(out, "FIRST PASS: %d\n", i);
							fflush(out);
				}

#ifdef F2MPI
				if (world.rank()) continue;
#endif	      	      		 

				if (/*ind->pars[0] || ind->pars[1] || */!ind->haplocount.size()) continue;

				// Perform the inversions indicated by the negshift data, at most a single one per individual
				// and chromosome, maybe 0.
				if (false) for (int c = 0; c < (int)chromstarts.size() - 1; c++)
				{
					int minstart = chromstarts[c + 1];
					double minval = -1e-10;
					bool prevlow = false;

					for (int p = chromstarts[c]; p < (int)chromstarts[c + 1]; p++)
					{
						if (ind->negshift[p] < minval)
						{
							minstart = p;
							minval = ind->negshift[p];
						}
						// TODO: Was a seed element ever needed here?
						/*						if (ind->negshift[p] < -1e-10)
						{
							if (!prevlow)
							{
								negshiftcand ourtuple(ind, 0, p);
								negshiftcands[c].insert(ourtuple);
							}
							prevlow = true;
							}*/
					}
					if (ind->lastinved[c] == minstart)
					{
						ind->lastinved[c] = -1;
						//continue;
					}

					ind->lastinved[c] = -1;

					if (minstart + 1 < chromstarts[c + 1])
					{
						negshiftcand ourtuple(ind, minval, minstart);
						inferiorrelated pred(ourtuple);
						for (set<negshiftcand>::iterator iter = negshiftcands[c].begin(); iter != negshiftcands[c].end(); )
						{
							if (pred(*iter))
							{
								negshiftcands[c].erase(iter++);
							}
							else
							{
								iter++;
							}

						}
						if (!pred.anymatch && rand() / (RAND_MAX / 5))
						{
							negshiftcands[c].insert(ourtuple);
						}
					}

				}
			}

#ifdef F2MPI
			if (!world.rank())
#endif
			{
			  bool hitnnn = false;
#pragma omp parallel for schedule(dynamic,1)
				for (unsigned int i = 0; i < INDCOUNT; i++)
				{

					individ* ind = getind(i, false);
					if (!ind || !ind->haplocount.size()) continue;

					fprintf(out, "SKEWNESS PASS: %d\n", i);
					fflush(out);

					int cno = 0;
					if (DOINFPROBS)
					for (unsigned int j = 0; j < ind->haplocount.size(); j++)
					{
						while (cno + 1 < chromstarts.size() && j >= chromstarts[cno + 1]) cno++;

						for (int side = 0; side < 2; side++)
						{
							double bestprob = 0;
							MarkerVal bestmarker = UnknownMarkerVal;
							double sum = 0;


							for (auto probpair : ind->infprobs[j][side])
							{
								sum += probpair.second;
								if ((ind->n == 433 && j >= 4086 && j <= 4087)) fprintf(stdout, "PROBPAIR A: %d %d %d %d %lf\n", ind->n, j, side, probpair.first.value(), probpair.second);
							}

							MarkerVal priorval = UnknownMarkerVal;
							if (ind->priormarkerdata.size() > j)
							{
								priorval = (&ind->priormarkerdata[j].first)[side];
							}							

							for (auto probpair : ind->infprobs[j][side])
							{
								double curprob = 0.5;
								auto curmarker = (&ind->markerdata[j].first)[side];

								if (curmarker != UnknownMarkerVal)
								{
									curprob = fabs((curmarker == probpair.first ? 1 : 0) - (&ind->markersure[j].first)[side]);
								}

								double d = (probpair.second - sum * curprob) / (curprob - curprob * curprob);
								d += log( 1 / curprob - 1);

								if (priorval != UnknownMarkerVal)
								{
									double priorprob = 1.0 - (&ind->priormarkersure[j].first)[side];

									MarkerVal nowval = probpair.first;
									if (nowval != priorval)
									{
										priorprob = 1.0 - priorprob;
									}							

									d += log(priorprob) - log(1 - priorprob);
								}
								if (d)
								{
									ind->infprobs[j][side][probpair.first] = caplogitchange(curprob + d * scalefactor, curprob, maxdiff / (ind->children + 1), hitnnn);
								}
							}

							for (auto probpair : ind->infprobs[j][side])
							{
								if (probpair.second > bestprob)
								{
									bestmarker = probpair.first;
									bestprob = probpair.second;
								}
								if ((ind->n == 433 && j >= 4086 && j <= 4087)) fprintf(stdout, "PROBPAIR B: %d %d %d %d %lf\n", ind->n, j, side, probpair.first.value(), probpair.second);
							}

							if (bestmarker != UnknownMarkerVal || bestprob > 0)
							{
								(&ind->markerdata[j].first)[side] = bestmarker;
								double intended = 1.0 - bestprob;
								(&ind->markersure[j].first)[side] = intended;
							}
							ind->infprobs[j][side].clear();
						}
						// oldinfprobslogic(ind, j, iter, cno, out);
					}

					{
						vector<bool> allhalf;
						vector<bool> anyinfo;
						vector<bool> cleared;
						vector<int> nudgeme;


						anyinfo.resize(chromstarts.size());
						allhalf.resize(chromstarts.size());
						cleared.resize(chromstarts.size());
						nudgeme.resize(chromstarts.size());
						for (int k = 0; k < chromstarts.size(); k++)
						{
							anyinfo[k] = false;
							allhalf[k] = true;
							cleared[k] = false;
							nudgeme[k] = -1;
						}

						cno = 0;
						double prevval = 0.5;
						for (unsigned int j = 0; j < ind->haplocount.size(); j++)
						{
							while (cno + 1 < chromstarts.size() && j >= chromstarts[cno + 1]) cno++;
							anyinfo[cno] = true;

							if (!(ind->haploweight[j] && ind->haploweight[j] != 1) && false)
							{
								double b1 = ind->haplobase[j];
								double b2 = ind->haplocount[j] - ind->haplobase[j];
								/*if (!allhalf[cno])
								{
								ind->haploweight[j] = 0.5 - (0.5 - ind->haploweight[j]) * 0.01;
								}*/

								if (!early && fabs(ind->negshift[j]) < 1e-5)
								{
									fprintf(out, "Clearing: %d %lf %lf %lf\n", ind->n, b1, b2, ind->negshift[j]);
									//lockhaplos(ind, cno);
									cleared[cno] = true;
									ind->haploweight[j] = 0.5 /*- (0.5 - ind->haploweight[j]) * 0.99*/;
									ind->haplocount[j] = 0;
								}
								else
									allhalf[cno] = false;
							}


							if ((ind->haplocount[j] || RELSKEWS) && ind->haploweight[j] && ind->haploweight[j] != 1 /*&& (!ind->founder || ind->children)*/)
							{
							  
							  double val;
							  if (ind->haplocount[j])
							    {
							      //val = exp(ind->haplobase[j] / ind->haplocount[j]);
								  val = ind->haplobase[j] / ind->haplocount[j];
							      val *= (1 - ind->haploweight[j]) / ind->haploweight[j];
							    }
							  else
							    {
							      val = 1;
							    }

								double baseterm = log(ind->haploweight[j] / (1 - ind->haploweight[j]));
								double relskewterm = 0;
								if (false && RELSKEWS && j)
								{
									// Modify haplotype based on relhaplo relationship with raw intended haploweight for j - 1
									double relval = ind->relhaplo[j - 1];
									double sum = 0;
									double term = ind->haploweight[j] * (prevval * relval + (1 - prevval) * (1 - relval));
									double lo = term;
									sum += term;

									term = (1 - ind->haploweight[j]) * ((1 - prevval) * relval + prevval * (1 - relval));
									sum += term;
									
									if (sum)
									{
										lo /= sum;
										relskewterm = log((lo + 1e-300) / ((1 - lo) + 1e-300)) - baseterm;
									}

									prevval = exp((log(val) * ind->haplocount[j] + relskewterm) + baseterm);
									prevval = prevval / (prevval + 1.0);
								}

								double scorea = 1.0 - ind->markersure[j].first;
								double scoreb = 1.0 - ind->markersure[j].second;
								if (ind->markerdata[j].first != ind->markerdata[j].second) scoreb = 1 - scoreb;

								double similarity = scorea * scoreb + (1 - scorea) * (1 - scoreb);
								if (similarity >= 1 - maxdiff)
								{
									ind->haplobase[j] = 0;
								}
								else
								{
									double count = ind->haplocount[j];
									ind->haplobase[j] -= count * ind->haploweight[j];
									count = count - similarity * count;
									ind->haplobase[j] += count * ind->haploweight[j];
									ind->haplobase[j] *= ind->haplocount[j] / count;
									if (ind->haplobase[j] < 0) ind->haplobase[j] = 0;
									if (ind->haplobase[j] >= count) ind->haplobase[j] = count;
								}

								
								double intended = ind->haploweight[j] + scalefactor * ((ind->haplobase[j] - ind->haploweight[j] * ind->haplocount[j]) / (ind->haploweight[j] - ind->haploweight[j] * ind->haploweight[j]) + log(1/ind->haploweight[j] - 1));

								if (!early && allhalf[cno] && fabs(intended - 0.5) > 0.1 &&
									ind->markerdata[j].first != UnknownMarkerVal && ind->markerdata[j].second != UnknownMarkerVal &&
									cleared[cno])
								{
									allhalf[cno] = false;
									fprintf(out, "Locking: %d %d %lf\n", ind->n, j, ind->negshift[j]);
									ind->haploweight[j] = (intended < 0.5) ? 0 : 1;
								}
								else
								{
								  if (/*ind->children &&*/ (ind->lastinved[cno] == -1 || true) /*&& !ind->pars[0] && !ind->pars[1]*/)
								    {
									// Cap the change if the net difference is small/miniscule
									  intended = caplogitchange(intended, ind->haploweight[j], maxdiff / (ind->children + 1), hitnnn);

									//								if ((ind->haploweight[j] - 0.5) * (intended - 0.5) < 0) intended = 0.5;
									
									
									  /*										if (!(intended < 0.5) && ind->haploweight[j] < 0.5)
										{
											cout << "CROSSOVER " << ind->name << " " << ind->n << " " << j << " " << intended << " " << ind->haploweight[j] << " " << limn << " " << limd1 << std::endl;
											}*/
									  ind->haploweight[j] = intended;


										// Nudging flag currently not respected
										if ((nudgeme[cno] == -1 || fabs(ind->haploweight[nudgeme[cno]] - 0.5) < fabs(ind->haploweight[j] - 0.5)) && ind->haploweight[j] > maxdiff && ind->haploweight[j] < 1 - maxdiff)
										{
											nudgeme[cno] = j;
										}
									}
								}

								/*							if (ind->haploweight[j] != 0.5)
								{
								allhalf[cno] = false;
								}*/
							}
						}
					}
					vector<pair<double, boost::tuple<individ*, individ*, int, int> > > allnegshifts;
					flat_map<individ*, double> bestshift;
					for (auto i = nsm.begin(); i != nsm.end(); i++)
					{
						for (auto j = i->second.begin(); j != i->second.end(); j++)
						{
							// 1-3 allows shifts, but not genotype switches
							// 2 only allows shifts for paren 2, e.g. assumption that paren 1 is part of several half sibships
							for (int k = 1; k <= 4; k++)
							{
								/*		      		      */
								if (k == 2/* || k == 4*/)
									allnegshifts.push_back(make_pair(j->second[k] - 1e-5, boost::make_tuple(i->first.first, i->first.second, k, j->first)));
								/*		      if (k == 4) printf("ANS: %lf %d %d %d %d\n",
								j->second[k], i->first.first->n, i->first.second->n, k, j->first);*/

								/*		      if (k == 4)
								allnegshifts.push_back(make_pair(j->second[k] - 1e-5, make_tuple(i->first.first, i->first.first, 1, j->first)));*/
							}
						}
					}

					sort(allnegshifts.begin(), allnegshifts.end());

					for (int k = allnegshifts.size() - 1; k >= 0; k--)
					{
						if (bestshift[allnegshifts[k].second.get<0>()] < allnegshifts[k].first &&
							bestshift[allnegshifts[k].second.get<1>()] < allnegshifts[k].first)
						{
							bestshift[allnegshifts[k].second.get<0>()] = allnegshifts[k].first;
							bestshift[allnegshifts[k].second.get<1>()] = allnegshifts[k].first;
							int c = 0;
							while (c < chromstarts.size() && (chromstarts[c] <= allnegshifts[k].second.get<3>()))
							{
								c++;
							}

							int phase = allnegshifts[k].second.get<2>();
							individ* inds[2] = { allnegshifts[k].second.get<0>(), allnegshifts[k].second.get<1>() };
							if (rand() > RAND_MAX / 10) continue;

							printf("Inv: %d %d %d %d %lf\n", inds[0]->n, inds[1]->n, phase, allnegshifts[k].second.get<3>(), allnegshifts[k].first);
							for (int m = allnegshifts[k].second.get<3>() + 1; m < chromstarts[c]; m++)
							{
								for (int z = 0; z < 2; z++)
								{
									if (phase & 4)
									{
										// These days, we are just emulating the shifts
										if (z == 0 && false)
										{
											swap(inds[0]->haploweight[m], inds[1]->haploweight[m]);
											swap(inds[0]->markerdata[m], inds[1]->markerdata[m]);
											swap(inds[0]->markersure[m], inds[1]->markersure[m]);
										}

										for (int k = 0; k < inds[z]->kids.size(); k++)
										{
											// ONLY inverting those that share both parents
											// This used to be that we would invert every child, but avoid inverting those that share both parents
											// twice, i.e. not inverting them at all
											if (z && inds[z]->kids[k]->pars[0] == inds[0]) {
												//printf("I2: %d %d\n", inds[z]->kids[k]->n, m);
												inds[z]->kids[k]->haploweight[m] = 1.0 - inds[z]->kids[k]->haploweight[m];
											}
										}
									}

									if (phase & (1 << z))
									{
										//			      printf("Inv %d at %d, was %lf\n", inds[z]->n, m, inds[z]->haploweight[m]);
										inds[z]->haploweight[m] = 1.0 - inds[z]->haploweight[m];
									}

								}
							}
						}
					}




				}
				if (hitnnn)
				  {
				    scalefactor /= 1.1;
				  }
				else
				  {
				    scalefactor *= 1.1;
				  }
				if (scalefactor < 0.01) scalefactor = 0.01;
				fprintf(stdout, "Scale factor now %lf\n", scalefactor);
				for (int c = 0; c < (int) chromstarts.size() - 1; c++)
				  {
				    for_each(negshiftcands[c].begin(), negshiftcands[c].end(), negshifter(c));

				  }
			}
		}
	}

	/*	FILE* binut = fopen("form_data/set009/ad_dataset_009_chr_01.bin", "wb");
	for (int pos = (int) markerposes[chromstarts[0]]; pos != (int) markerposes[chromstarts[1] - 1] + 1; pos++)
	{
	for (int half = 0; half < 2; half++)
	{
	for (int ind = 0; ind < dous.size(); ind++)
	{
	fwrite(&realgeno[ind][pos][half], sizeof(float), 1, binut);
	}
	}
	}*/

	generation++;
}

// Handle the fuss of trailing \r, \n characters when combining scanf and gets and stuff.
void clean(char* tlf)
{
	int i = strlen(tlf);
	while (i && tlf[i - 1] < 32)
	{
		tlf[i - 1] = 0;
		i--;
	}
}

void readhaplodata(int innum)
{
	char tlf[255];
	sprintf(tlf, "qtlmas14haplo_%d", innum);
	FILE* in = fopen(tlf, "r");
	bool inind = false;
	bool inact = false;
	individ* ind = 0;
	int mcc = 0;
	int mnum;

	while (fgets(tlf, 255, in))
	{
		float w;
		int a, b;
		float n;

		if (sscanf(tlf, "%f %d %d %f", &w, &a, &b, &n) == 4)
		{
			if (ind)
			{
				ind->haploweight[mnum++] = w;
			}
		}
		else
		{
			ind = 0;
			float a1, a2;

			if (sscanf(tlf, "%d", &a) == 1)
			{
				ind = getind(a);
				printf("Reading ind %d, was at marker %d\n", a, mnum);
				mnum = innum;
			}
			else
				if (sscanf(tlf, "actrec: %f %f", &a1, &a2) == 2)
				{
					if (!inact)
					{
						inact = true;
						mnum = innum;
					}

					mnum++;
					actrec[0][mnum] = a1 / 10000.0;
					actrec[1][mnum] = a2 / 10000.0;
					printf("Reading actrec %d : %f, %f\n", mnum, a1, a2);
				}
				else
					inact = false;
		}
	}
}

// Map individual names to inds
map<string, individ*> indmap;
individ* getind(string name, bool create = true)
{
	decltype(indmap.begin()) i;
	while ((i = indmap.find(name)) == indmap.end() && create)
	{
		int origindex = indmap.size();
		individ* ind = indmap[name] = getind(origindex, true);
		if (ind) ind->name = name;
	}

	return i == indmap.end() ? nullptr : i->second;
}

const individ* zeroguy = getind("0");

void readalphaped(FILE* in)
{
	char me[255], father[255], mother[255];
	while (fscanf(in, "%s %s %s", me, father, mother) == 3)
	{
		char line[255];
		fgets(line, 255, in);
		int gen = 0;
		if (sscanf(line, "%d", &gen) != 1) gen = 0;

		// We are all empty until we get some data, not implemented in getind yet due
		// to possible regressions.
		individ* ime = getind(me);
		individ* ifounderf = getind(father);
		individ* ifounderm = getind(mother);

		if (ime) ime->empty = true;
		if (ifounderf) ifounderf->empty = true;
		if (ifounderm) ifounderm->empty = true;

		if (gen >= 2 && !ifounderf->gen && !ifounderm->gen)
		{
			individ* realpars[2] = { getind(me + (string) "_aux_realf"), getind(me + (string) "_aux_realm") };
			for (int k = 0; k < 2; k++)
			{
				ime->pars[k] = realpars[k];
				realpars[k]->gen = 1;
				realpars[k]->pars[0] = ifounderf;
				realpars[k]->pars[1] = ifounderm;
				realpars[k]->empty = true;
			}
			ime->gen = gen;
		}
		else
		{
			ime->gen = gen;
			ime->pars[0] = ifounderf;
			ime->pars[1] = ifounderm;
		}

		if (gen >= 2)
		{
			dous.push_back(ime);
		}
	}
}

void readalphadata(FILE* in)
{
	individ* haplo = getind("haplo");
	for (int x = 0; x < markerposes.size(); x++)
	{
		haplo->markerdata[x] = make_pair(sexmarkerval, sexmarkerval);
		haplo->markersure[x] = make_pair(0, 0);
	}

	char me[255];
	while (fscanf(in, "%s", me) == 1)
	{
		individ* ime = getind(me);
		bool doublehaplo = (ime->pars[1] == getind("haplo"));
		ime->empty = false;
		int data;
		for (int x = 0; x < markerposes.size(); x++)
		{
			char datastr[255];
			fscanf(in, "%s", datastr);
			int data2;
			int numread = sscanf(datastr, "%d/%d", &data, &data2);
			ime->haploweight[x] = 0.5;
			if (numread == 1)
			{
				pair<MarkerVal, MarkerVal> marker;
				switch (data)
				{
				case 0:
					marker = make_pair(1 * MarkerValue, 1 * MarkerValue);
					break;
				case 1:
					marker = make_pair(1 * MarkerValue, 2 * MarkerValue);
					break;
				case 2:
					marker = make_pair(2 * MarkerValue, 2 * MarkerValue);
					break;
				default:
					marker = make_pair(UnknownMarkerVal, UnknownMarkerVal);
					break;
				}
				ime->markerdata[x] = marker;
				ime->markersure[x] = make_pair(0.02, 0.02);
			}
			else
			{
				printf("%s %d\n", datastr, numread);
				if (data == data2 && !data)
				{
					ime->markerdata[x] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
					ime->markersure[x] = make_pair(0, 0);
				}
				else
				{
					boost::math::binomial_distribution<double> dist1(data);
					boost::math::binomial_distribution<double> dist2(data2);
					double sure1 = 0;
					double sure2 = 0;
					double probsum = 0;

					for (int rl1 = 0; rl1 <= data; rl1++)
					{
						for (int rl2 = 0; rl2 <= data2; rl2++)
						{
							int l1 = rl1;
							int l2 = rl2;
							double overallprob = (data ? boost::math::pdf(dist1, l1) : 1) * (data2 ? boost::math::pdf(dist2, l2) : 1);

							double sureb1;
							double sureb2;
							while (true)
							{
								sureb1 = 0.5;
								sureb2 = 0.5;
								if (l1 + l2)
									sureb1 = l1 / (double)(l1 + l2);

								if (data + data2 - l1 - l2)
									sureb2 = (data2 - l2) / (double)(data + data2 - l1 - l2);

								if (sureb1 + 1e-9 > 1 - sureb2) break;
								l1 = data - l1;
								l2 = data2 - l2;
							}

							overallprob *= pow(sureb1, l1) * pow(1 - sureb1, l2) * pow(sureb2, (data2 - l2)) * pow(1 - sureb2, data - l1);

							sure1 += sureb1 * overallprob;
							sure2 += sureb2 * overallprob;
							probsum += overallprob;
						}
					}

					sure1 /= probsum;
					sure2 /= probsum;

					pair<MarkerVal, MarkerVal> marker = make_pair(2 * MarkerValue, 1 * MarkerValue);
					pair<double, double> markersure = make_pair(sure1, sure2);

					// Oh, error rate over 0.5, let's invert to make it more passable
					for (int k = 0; k < 2; k++)
					{
						if ((&(markersure.first))[k] > 0.5)
						{
							(&(markersure.first))[k] = 1 - (&(markersure.first))[k];
							(&(marker.first))[k] = (k + 1) * MarkerValue;
						}
					}

					ime->markerdata[x] = marker;
					ime->markersure[x] = markersure;
					printf("%d/%d turned into %d %d with %lf;%lf\n", data, data2, marker.first.value(), marker.second.value(), markersure.first, markersure.second);
				}
				if (doublehaplo)
				{
					ime->markerdata[x] = make_pair(ime->markerdata[x].first, sexmarkerval);
				}
			}
		}
		ime->priormarkerdata = ime->markerdata;
		ime->priormarkersure = ime->markersure;
	}
}

void readalphamap(FILE* in)
{
	double prevval = 1e30;
	double value;
	int nowmarker = 0;
	while (fscanf(in, "%lf", &value) == 1)
	{
		if (value < prevval)
		{
			chromstarts.push_back(nowmarker);
		}
		nowmarker++;
		markerposes.push_back(value);
		prevval = value;
	}
	chromstarts.push_back(markerposes.size());
}

void readmerlinmap(FILE* mapfile)
{
	// Assuming sex-averaged map
	// Assuming chromosomes and markers to be presented in order
	int chromo;
	int oldchromo = -1;
	double cmpos; // All positions will be shifted making first marker occurs at 0
	double cmbase = -1;
	int bppos;
	char markername[256];
	while (fscanf(mapfile, "%d %255s %lf %d", &chromo, markername, &cmpos, &bppos) == 4)
	{
		chromo--;
		if (chromo != oldchromo)
		{
			chromstarts.push_back(markerposes.size());
			oldchromo = chromo;
			cmbase = cmpos;
		}

		markerposes.push_back(cmpos - cmbase);
		for (int t = 0; t < 2; t++)
		{
			actrec[t].push_back(baserec[t]);
		}

	}

	chromstarts.push_back(markerposes.size());
}

void readmerlinped(FILE* pedfile)
{
	indmap["0"] = 0;
	char famname[256];
	char indname[256];
	char pname[256];
	char mname[256];
	int sex;
	double pheno;
	while (fscanf(pedfile, "%255s %255s %255s %255s %d %lf", famname, indname, pname, mname, &sex, &pheno) == 6)
	{
		printf("Reading main line for individual %s\n", indname);
		individ* ind = getind(indname);
		ind->pars[0] = getind(pname);
		ind->pars[1] = getind(mname);
		ind->gen = 0;
		ind->sex = --sex;

		if (ind->pars[0] || ind->pars[1])
		{
			dous.push_back(ind);
			ind->gen++;
		}

		ind->markerdata.resize(markerposes.size());
		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());

		ind->infprobs.resize(markerposes.size());
		ind->markersure.resize(markerposes.size());
		//		ind->semishift.resize(5000);
		ind->lastinved.resize(chromstarts.size());
		ind->lockstart.resize(chromstarts.size());

		for (int i = 0; i < chromstarts.size(); i++)
		{
			ind->lastinved[i] = -1;
			ind->lockstart[i] = 0;
		}

		int a, b;
		for (int k = 0; k < markerposes.size(); k++)
		{
			ind->haploweight[k] = 0.5;
			fscanf(pedfile, "%d %d", &a, &b);

			// Assume 1 percent allele error rate
			ind->markersure[k] = make_pair(a ? 1e-7 : 0.0, b ? 1e-7 : 0.0);

			ind->markerdata[k] = make_pair(a * MarkerValue, b * MarkerValue);

		}
	}
}

int family = 1;

void domerlinind(FILE* pedfile, individ* ind)
{
	int pn[2] = { 0 };

	for (int j = 0; j < 2; j++)
	{
		if (ind->pars[j]) pn[j] = ind->pars[j]->n;
	}

	fprintf(pedfile, "%d\t%d\t%d\t%d\t%d", family, ind->n, pn[0], pn[1], ind->sex + 1);
	for (int k = 0; k < chromstarts[1]; k++)
	{
		fprintf(pedfile, "\t%d\t%d", ind->markerdata[k].first.value(), ind->markerdata[k].second.value());
	}

	fprintf(pedfile, "\n");
}

template<class RuleType, class AttrType> void parseToEndWithError(istream& file, const RuleType& rule, AttrType& target)
{
	auto parseriter = boost::spirit::istream_iterator(file);
	boost::spirit::istream_iterator end;

	bool res = phrase_parse(parseriter, end, rule, x3::space - x3::eol, target);

	if (!res)
	{
	  std::string val;
	  file >> val;
		throw logic_error("Parsing failed. " + (std::string) __func__ + " " + val);
	}

	if (!file.eof())
	{
		throw logic_error("Not reaching end of file in parser. " + (std::string) __func__);
	}
}

template<class RuleType> void parseToEndWithError(istream& file, const RuleType& rule)
{
	auto parseriter = boost::spirit::istream_iterator(file);
	boost::spirit::istream_iterator end;

	bool res = phrase_parse(parseriter, end, rule, x3::space - x3::eol);

	if (!res)
	{
	  std::string val;
	  while (!file.eof())
	    {
	      file >> val;
	      std::cout << val;
	    }
	  std::cout << std::endl;
		throw logic_error("Parsing failed. " + (std::string) __func__ + " # " + val);

	}

	if (!file.eof())
	{
		throw logic_error("Not reaching end of file in parser. " + (std::string) __func__);
	}
}

#ifdef READHAPSSAMPLE
std::string filterExisting(const set<std::string>& names, std::string name)
{
	if (names.find(name) == names.end())
	{
		return "0";
	}

	return name;
}

vector<int> mapIndices;
vector<bool> hapmonomorphs;

auto word_ = x3::lexeme[+(x3::char_ - x3::space)];

typedef std::vector<std::tuple<std::string, std::string, std::string>> sampletype;

const auto marker_ = (x3::int_ > word_);
const auto hapsLine = (marker_ > x3::omit[x3::float_] > word_ > word_ > (+x3::int_));
const auto hapsLineIgnoreGenotypes = (marker_ > x3::long_long > word_ > word_> x3::omit[(+x3::int_)]);

struct samplereader
{
	sampletype samples;

	void read(istream& sampleFile)
	{
		sampleFile >> std::noskipws;
		using namespace x3;

		auto sampleHeader = omit[
			repeat(7)[word_] > eol >
				int_ > int_ > int_ > repeat(4)[word_] > eol];
		auto sampleLine = (omit[int_] > word_ > omit[int_] > word_ > word_ > omit[int_] > omit[int_]);

		try
		{
			parseToEndWithError(sampleFile, sampleHeader > (sampleLine % eol), samples);
			std::cout << samples.size() << " samples read." << std::endl;
		}
		catch (expectation_failure<boost::spirit::istream_iterator> const& x)
		{
			std::cerr << "expected: " << x.which();
			std::cerr << "got: \"" << x.where() << '"' << std::endl;

			throw x;
		}
	}
};

void readhaps(const sampletype& samples, istream& bimFile, vector<istream*>& hapsFile)
{
	using namespace x3;

	bimFile >> std::noskipws;
	*hapsFile[0] >> std::noskipws;
	
	std::vector<std::tuple<int, std::string, std::string, std::string, std::vector<int>>> snpData;
	map<std::pair<int, std::string>, pair<int, int> > geneMap;
		
	auto alleles_ = word_ > word_;
	auto bimLine = (marker_ > omit[float_] > int_ > omit[alleles_]);

	try
	{
		parseToEndWithError(bimFile, (bimLine[([&](auto& context)
		{
			using namespace boost::fusion;
			auto attr = _attr(context);
			int index = geneMap.size();
			geneMap[make_pair(at_c<0>(attr), at_c<1>(attr))] = make_pair(at_c<2>(attr), index);
		})])
			% eol);
		std::cout << geneMap.size() << " entries read in map." << std::endl;

		parseToEndWithError(*hapsFile[0], hapsLine % eol, snpData);
		std::cout << snpData.size() << " SNPs read." << std::endl;
	}
	catch (expectation_failure<boost::spirit::istream_iterator> const& x)
	{
		std::cerr << "expected: " << x.which();
		std::cerr << "got: \"" << x.where() << '"' << std::endl;

		throw x;
	}

	int lastchrom = -1;
	double lastpos = -1;
	double basepos = 0;

	// Walk over SNPs to create map info
	for (auto snp : snpData)
	{
		int chrom = get<0>(snp);
		int bppos;
		int index;

		std::tie(bppos, index) = geneMap[make_pair(chrom, get<1>(snp))];

		double pos = bppos * 1e-6;
		
		if (chrom != lastchrom)
		{
			chromstarts.push_back(markerposes.size());
			basepos = pos;
		}

		markerposes.push_back(pos - basepos);
		mapIndices.push_back(index);
		lastpos = pos;
		lastchrom = chrom;
		hapmonomorphs.push_back(get<2>(snp) == get<3>(snp));
	}

	chromstarts.push_back(markerposes.size());

	vector<individ*> sampleInds;
	set<std::string> sampleNames;

	// Collect all names for which we really have samples
	for (std::tuple<std::string, std::string, std::string> sample : samples)
	{
		sampleNames.insert(get<0>(sample));
	}

	for (std::tuple<std::string, std::string, std::string> sample : samples)
	{
		individ* me = getind(get<0>(sample));
		// We don't care about sex? Is this OK?
		me->sex = 0;
		me->pars[0] = getind(filterExisting(sampleNames, get<1>(sample)));
		me->pars[1] = getind(filterExisting(sampleNames, get<2>(sample)));

		// Hack the generation to make non-founders full citizens
		me->gen = 2 * (me->pars[0] || me->pars[1]);
		dous.push_back(me);

		sampleInds.push_back(me);
	}

	auto dohaploweight = [] (individ* ind) { return (ind->gen < 2); };

	for (int i = 0; i < snpData.size(); i++)
	{
		const vector<int>& markers = get<4>(snpData[i]);
		for (int j = 0; j < sampleInds.size(); j++)
		{
		  float sureVal = 0;
		  /*if (sampleInds[j]->gen == 2)*/ sureVal = 0;
			sampleInds[j]->markerdata[i] = make_pair((markers[j * 2] + 1) * MarkerValue, (markers[j * 2 + 1] + 1) * MarkerValue);
			if (dohaploweight(sampleInds[j])) sampleInds[j]->haploweight[i] = 1e-3;
			sampleInds[j]->markersure[i] = { sureVal, sureVal };
			if (RELSKEWS)
			  {
				if (i != markerposes.size() - 1)
				{
					//sampleInds[j]->relhaplo[i] = 0.51;//(markers[j * 2] == markers[j * 2 + 1]) ? 1.0 : 0.51;
					// Assume that our local bp/cM model is 1e-6 and the population-level rho used by ShapeIT is 0.0004
					// And that a good proxy for haplotype accuracy as a function of length is the global rho...
					// Most important point is to try to get negshift inversion boundaries located on actual marker map gaps.
				  sampleInds[j]->relhaplo[i] = 0.5 + 0.5 * exp(-(markerposes[i + 1] - markerposes[i]) /* * 1e6 * 0.0004 */ );
				}
			  }
		}
	}

	const double padding = 0.01;
	double unit = 1.0 / (hapsFile.size() + padding);
	for (int j = 0; j < sampleInds.size(); j++)
	  {
	    for (int i = 0; i < snpData.size(); i++)
	      {
		if (RELSKEWS)
		  {
		    sampleInds[j]->relhaplo[i] = unit;
		  }
		sampleInds[j]->markersure[i] = make_pair(padding * unit, padding * unit);
	      }
	  }

	for (int k = 1; k < hapsFile.size(); k++)
	{
		vector<int> phases;
		vector<int> origPhases;
		phases.resize(sampleInds.size());
		origPhases.resize(sampleInds.size());
		snpData.clear();

		*hapsFile[k] >> std::noskipws;
		parseToEndWithError(*hapsFile[k], hapsLine % eol, snpData);
		std::cout << snpData.size() << " SNPs read." << std::endl;
		for (int i = 0; i < snpData.size(); i++)
		{
			const vector<int>& markers = get<4>(snpData[i]);
			for (int j = 0; j < sampleInds.size(); j++)
			{
				int oldPhase = phases[j];
				int matchNum;
				int numMatches = 0;

				for (int p = 1; p <= 2; p++)
				{
					auto marker = (p == 1) ? make_pair((markers[j * 2] + 1) * MarkerValue, (markers[j * 2 + 1] + 1) * MarkerValue) : make_pair((markers[j * 2 + 1] + 1) * MarkerValue, (markers[j * 2] + 1) * MarkerValue);
					if (marker == sampleInds[j]->markerdata[i])
					{
						matchNum = p;
						numMatches++;
					}
				}

				// Not conclusive, assume same phase
				if (numMatches == 2 || numMatches == 0)
				{
					matchNum = oldPhase;
				}

				phases[j] = matchNum;
				if (phases[j] && !origPhases[j])
				  {
				    origPhases[j] = phases[j];
				  }
				if (dohaploweight(sampleInds[j]) && origPhases[j] && origPhases[j] != phases[j])
				  {
				    sampleInds[j]->haploweight[i] += unit;
				  }
				if (RELSKEWS && i)
				{
				// The definition of relhaplo is that marker i
				// defines the skewness shift to the NEXT marker.
				  sampleInds[j]->relhaplo[i - 1] += unit * (oldPhase == 0 || phases[j] == oldPhase);
				}
				if (!numMatches)
				{
					// TODO: Look at whether one of the two alleles matches
					sampleInds[j]->markersure[i] = make_pair(sampleInds[j]->markersure[i].first + unit, sampleInds[j]->markersure[i].second + unit);
				}
			}
		}
	}

	for (int j = 0; j < sampleInds.size(); j++)
	{
		sampleInds[j]->priormarkerdata = sampleInds[j]->markerdata;
		sampleInds[j]->priormarkersure = sampleInds[j]->markersure;
	}
}

void createhapfile(const sampletype& samples, istream& oldhapfile, ostream& newhapfile)
{
	vector<individ*> sampleInds;
	for (auto sample : samples)
	{
		sampleInds.push_back(getind(get<0>(sample), false));
	}

	std::vector<std::tuple<int, std::string, long long, std::string, std::string>> snpData;

	using namespace x3;
	
	try
	{
		parseToEndWithError(oldhapfile, hapsLineIgnoreGenotypes % eol, snpData);
		std::cout << snpData.size() << " SNPs read." << std::endl;
	}
	catch (expectation_failure<boost::spirit::istream_iterator> const& x)
	{
		std::cerr << "expected: " << x.which();
		std::cerr << "got: \"" << x.where() << '"' << std::endl;

		throw x;
	}

	auto translateMarker = [](MarkerVal x) -> std::string
	{
		if (x == UnknownMarkerVal) return "?";
		if (x == 1 * MarkerValue) return "0";
		if (x == 2 * MarkerValue) return "1";
		assert(false);
	};

	int i = 0;
	for (auto marker : snpData)
	{
	  newhapfile << get<0>(marker) << " " << get<1>(marker) << " " << get<2>(marker) << " " << get<3>(marker) << " " << get<4>(marker);
		for (auto ind : sampleInds)
		{
			pair<MarkerVal, MarkerVal> data = ind->markerdata[i];
			if (ind->haploweight[i] > 0.5 && (((ind->pars[0] && !ind->pars[0]->empty) || (ind->pars[1] && !ind->pars[1]->empty)) || ind->children))
			{
				swap(data.first, data.second);
			}

			newhapfile << " " << translateMarker(data.first) << " " << translateMarker(data.second);
		}
		newhapfile << "\n";
		i++;
	}
}

void readfambed(std::string famFileName, std::string bedFileName, bool readall = true)
{
	using namespace boost::interprocess;
	using namespace x3;

	auto word = lexeme[+(char_ - space)];
	std::ifstream file(famFileName);
	file >> std::noskipws;
	auto parseriter = boost::spirit::istream_iterator(file);
	boost::spirit::istream_iterator end;
	map<std::string, int> indNums;

	phrase_parse(parseriter, end, (omit[int_] > word > omit[word] > omit[word] > omit[int_] > omit[int_])
		[([&](auto& ctx)
	{
	  int index = indNums.size();
		indNums[_attr(ctx)] = index;
	})] % eol, space - eol);
	cout << indNums.size() << " individuals found." << std::endl;

	int blocksize = (indNums.size() + 3) / 4;
	cout << blocksize << " bytes per SNP block." << std::endl;

	file_mapping bedFile(bedFileName.c_str(), read_only);

	mapped_region snpRegion(bedFile, read_only, 3); // Skip header
	unsigned char* snpdata = (unsigned char*)snpRegion.get_address();
	size_t size = snpRegion.get_size();
	cout << size << " bytes, " << size/blocksize << " SNPs." << std::endl;

	vector<int> indArray;
	for (individ* ind : dous)
	{		
	        indArray.push_back(indNums[ind->name]);
		cout << ind->name << " " << indArray[indArray.size()-1] << std::endl;
	}

	for (int i = 0; i < markerposes.size(); i++)
	{
		unsigned char* thisSnp = &snpdata[mapIndices[i] * blocksize];
		for (int j = 0; j < dous.size(); j++)
		{
		  /*if (!(
			dous[j]->pars[0] && !dous[j]->pars[0]->empty &&
			dous[j]->pars[1] && !dous[j]->pars[1]->empty)) continue;*/
		  
			int index = indArray[j];
			int thisval = (thisSnp[index / 4] >> (2 * (index % 4))) & 3;
			pair<MarkerVal, MarkerVal> marker;

			bool isachange = false;
			switch (thisval)
			{
			case 0:
				marker = make_pair(1 * MarkerValue, 1 * MarkerValue);
				isachange = marker != dous[j]->priormarkerdata[i];
				break;
			case 1:
				marker = make_pair(UnknownMarkerVal, UnknownMarkerVal);
				cout << "MISSING BED " << dous[j]->name << " " << i << std::endl;
				break;
			case 2:			  
				marker = make_pair(1 * MarkerValue, 2 * MarkerValue);
				isachange = dous[j]->priormarkerdata[i].first == dous[j]->priormarkerdata[i].second;
				break;
			case 3:
				// ShapeIT will turn A A to 0 A, making all genotypes homozygotes for the second allele, rather than the first
				int val = 2 - hapmonomorphs[i];
				marker = make_pair(val * MarkerValue, val * MarkerValue);
				isachange = marker != dous[j]->priormarkerdata[i];
				break;
			}
			if (isachange)
			{
			  cout << "!!! " << dous[j]->name << " " << i << " " << marker.first.value() << marker.second.value() << " " << std::endl;
			}

			/*if (readall || marker.first == UnknownMarkerVal || isachange)
			  {
			    //cout << "/// " << dous[j]->name << " " << i << std::endl;
			    dous[j]->priormarkerdata[i] = marker;
			    if (marker.first == UnknownMarkerVal)
			      {
				/*if (RELSKEWS)
				  {
				    dous[j]->relhaplo[i] = 1;
				    }*/
				/*dous[j]->priormarkersure[i] = make_pair(0.f, 0.f);
			      }
			  }*/
			if (isachange || marker.first == UnknownMarkerVal)
			{
				dous[j]->priormarkersure[i] = make_pair(
					0.5 * (0.5 + dous[j]->priormarkersure[i].first),
					0.5 * (0.5 + dous[j]->priormarkersure[i].second));
				cout << "Increasing prior uncertainty individual " << dous[j]->n << ", marker " << i << std::endl;
			}
			if (dous[j]->priormarkerdata[i].first != UnknownMarkerVal)
			  {
			    dous[j]->priormarkersure[i] = make_pair(
								    max(1e-3, dous[j]->priormarkersure[i].first),
								    max(1e-3, dous[j]->priormarkersure[i].second));			  
			  }
		}
	}
}
#endif

void compareimputedoutput(istream& filteredOutput)
{
  int j = 0;
	while (!filteredOutput.eof())
	{
	  std::cout << j++ << std::endl;
		for (individ* ind : dous)
		{
			std::string name;
			do			  
			filteredOutput >> name;
			while (name == "--");
			  
			for (int i = 0; i < chromstarts[1]; i++)
			{
				double val[3];
				int maxval = 0;
				for (int k = 0; k < 3; k++)
				{
				  std::string temp;
				  filteredOutput >> temp;
				  if (sscanf(temp.c_str(), "%lf", &val[k]) != 1)
				    {
				      val[k] = -1;
				    }
					if (val[k] > val[maxval]) maxval = k;
				}

				int data = (ind->markerdata[i].first == 2 * MarkerValue) + (ind->markerdata[i].second == 2 * MarkerValue);
				int origmaxval = maxval;
				int oridata = data;

				// If reference allele is not aligned in shapit output:
				/*				  maxval = abs(1-maxval);
								  data = abs(1-data);*/
				  // That will only look at hetero vs. homo, not what hetero
				if (maxval != data &&
			      ind->pars[0] && !ind->pars[0]->empty &&
			      ind->pars[1] && !ind->pars[1]->empty &&
				    i != chromstarts[1] - 1 && val[origmaxval] >= 0 && ind->markerdata[i].first != UnknownMarkerVal)
				{
				  std::cout << ind->name << " " << j << ":" << i << " " << oridata << "\t";
					for (double oneVal : val)
					{
						std::cout << oneVal << "\t";
					}
					std::cout << std::endl;
				}
			}
		}
	}
}

void readhaplodata(FILE* in, int swap)
{
	char tlf[255];
	bool inind = false;
	bool inact = false;
	individ* ind = 0;
	int mnum = 0;
	int mcc = 0;
	int ncc = 0;

	int mixes[] = {
	  0, 0,
	  1215, 1739,
	  1241, 2254,
	  1264, 1886,
	  1314, 1905,
	  1332, 1485,
	  1351, 1430,
	  1391, 1758,
	  1457, 1576,
	  1425, 1860,
	  1482, 1823,
	  1527, 1564,
	  1538, 1855,
	  1552, 2259,
	  1590, 2057,
	  1694, 2101,
	  1772, 1931,
	  1773, 2096,
	  1805, 2020,
	  1849, 1900,
	  1876, 1935,
	  1877, 2076,
	  1936, 2029,
	  1953, 2288,
	  2012, 2233,
	  2034, 2293,
	  2118, 2237,
	  0, 0 };

	printf("Reading haplos\n");
	while (fgets(tlf, 255, in))
	{
		float w;
		int a, b;
		float n;

		if (sscanf(tlf, "%f %d %d %f", &w, &a, &b, &n) == 4)
		{
			if (ind)
			{
				int i = mnum;
				if ((!((a == ind->markerdata[i].first.value() && b == ind->markerdata[i].second.value()) ||
					(b == ind->markerdata[i].first.value() && a == ind->markerdata[i].second.value()))) && (a || b))
				{
					printf("Mismatch: %d %d\t%d %d %d %d\n", ind->n, mnum, a, b, ind->markerdata[i].first.value(), ind->markerdata[i].second.value());
					mcc++;
				}
				if (!(a || b)) ncc++;

				ind->haploweight[mnum++] = w;
			}
		}
		else
		{
			float a1, a2;

			if (sscanf(tlf, "%d", &a) == 1)
			{
				int i = 0;
				int n = a;
				/*      if (n <= 220 && swap)
	  {
		if (n > 20)
	  {
		n -= 21;
		n /= 10;
		n++;
	  }
		else
	  {
		n--;
		n *= 10;
		n += 21;
	  }
	  }*/
				a = n;
				while (mixes[i])
				{
					for (int j = 0; j < 2; j++)
					{
						if (a == mixes[i + j])
						{
							a = mixes[i + (!j)];
							break;
						}
					}
					i += 2;
				}

				if (mnum)
				{
					printf("%d Ind %d, %d mismatches (%d)\n", swap, ind->n, mcc, ncc);
				}

				ind = getind(a);
				printf("Reading ind %d, was at marker %d\n", a, mnum);
				if (a > 220) ind = 0;

				mnum = 0;
				mcc = 0;
				ncc = 0;
			}
			else
				if (sscanf(tlf, "actrec: %f %f", &a1, &a2) == 2)
				{
					if (!inact)
					{
						inact = true;
						mnum = 0;
					}

					mnum++;
					actrec[0][mnum] = a1 / 10000.0;
					actrec[1][mnum] = a2 / 10000.0;
					printf("Reading actrec %d : %f, %f\n", mnum, a1, a2);
				}
				else
					inact = false;
		}
	}
}

void deserialize(istream& stream)
{
	/*
	This is what we're unserializing
	fprintf(stdout, "%f\t%d\t%d\t\t%f\t%lf %lf %lf\n", ind->haploweight[j], ind->markerdata[j].first.value(), ind->markerdata[j].second.value(), ind->negshift[j],
												ind->markersure[j].first, ind->markersure[j].second, RELSKEWS ? ind->relhaplo[j] : 0);*/
  auto haploline = x3::double_ >> x3::int_ >> x3::int_ >> x3::double_ >> x3::double_ >> x3::double_;
	
	while (!stream.eof())
	{
		string line;
		std::getline(stream, line);

		pair<int, string> data;
	
		if (x3::phrase_parse(line.begin(), line.end(), x3::int_ >> word_ >> x3::eoi, x3::space, data))
		{
			int n;
			string name;
			std::tie(n, name) = data;

			individ* ind = getind(n, false);
			individ* indcheck = getind(name, false);

			if (ind && ind == indcheck)
			{
				int oldphase = 0;
				int switches = 0;
				for (int i = 0; i < markerposes.size(); i++)
				{
					std::getline(stream, line);
					std::tuple<double, int, int, double, double, double> output;
					if (!x3::phrase_parse(line.begin(), line.end(), haploline, x3::space, output))
					{
						std::cerr << "Reading haplotype for marker " << i << " for individual " << ind->name << " failed: " << line << std::endl;
					}
					else
					{
					  ind->haploweight[i] = std::get<0>(output);
						if (ind->haploweight[i] == 0.5) continue;

						int newphase = 1 + (ind->haploweight[i] > 0.5);
						if (oldphase && oldphase != newphase) switches++;

						oldphase = newphase;

						pair<MarkerVal, MarkerVal> pmv = make_pair(std::get<1>(output) * MarkerValue, std::get<2>(output) * MarkerValue);
						if (pmv != ind->markerdata[i])
						{
							std::cerr << "Genotype mismatch for marker " << i << " for individual " << ind->name << " (" << ind->markerdata[i].first.value() << "," << ind->markerdata[i].second.value() << ") to " <<
								" (" << pmv.first.value() << "," << pmv.second.value() << ")" << std::endl;
						}
						ind->markerdata[i] = pmv;
						ind->markersure[i] = make_pair(std::get<4>(output), std::get<5>(output));
					}
				}

				if (ind->children || (ind->pars[0] && !ind->pars[0]->empty) || (ind->pars[1] && !ind->pars[1]->empty)) std::cout << "Switches " << ind->n << " " << ind->name << "\t" << switches << std::endl;
			}
			else
			{
				std::cerr << "Supposed individual header not a header: " << line << std::endl;
			}
		}
	}
}

int main(int argc, char* argv[])
{
#ifdef _MSC_VER
	// Performant?
	// This turns off some aspects of IEEE precision, but we do not really need it, as we avoid denormalized values
	// anyway.
	_controlfp(_DN_FLUSH, _MCW_DN);
	_set_printf_count_output(1);
#endif
#ifdef F2MPI
	mpi::environment env(argc, argv);
	mpi::communicator world;
#endif
	discstep = 1;
	sexc = 2;
	if (NUMGEN == 3)
	{
		baserec[0] = -discstep / 50.0;
	}
	else
	{
		baserec[0] = -discstep / 50.0;
	}
	baserec[1] = baserec[0];
	for (int gen = 0; gen < 2; gen++)
	{
		genrec[gen] = baserec[0];
	}
	//selfgen = 1;
	genrec[2] = baserec[0] /** selfgen*/;
	 
#ifdef READHAPSSAMPLE
	po::options_description desc;
	po::variables_map inOptions;
	string impoutput, famfilename, bedfilename, deserializefilename, outputfilename, outputhapfilename, genfilename, pedfilename, mapfilename, samplefilename;
	int COUNT;

	desc.add_options()("samplefile", po::value<string>(&samplefilename), "ShapeIT-style .sample file")
		("bimfile", po::value<string>(), "BIM file")
		("hapfiles", po::value<vector<string> >()->multitoken(), "One or more HAP files, maximum realization followed by others.")
		("deserialize", po::value<string>(&deserializefilename), "Load existing Chaplink output as starting point, with reporting on number of inversions.")
		("impoutput", po::value<string>(&impoutput), "Imputed genotype output from previous run.")
		("famfile", po::value<string>(&famfilename), "Original PLINK fam file. Use with bedfile.")
		("bedfile", po::value<string>(&bedfilename), "Original PLINK bed file. Use with famfile.")
		("count", po::value<int>(&COUNT)->default_value(3), "Number of iterations")
		("output", po::value<string>(&outputfilename), "Output file name")
		("capmarker", po::value<int>()->notifier([&](int cap)
	{
		markerposes.resize(cap);
		chromstarts[1] = min(cap, (int)chromstarts[1]);
	}), "Limit to marker count.")
		("mapfile", po::value<string>(&mapfilename), "map file in original PlantImpute format, similar to AlphaImpute.")
	        ("pedfile", po::value<string>(&pedfilename), "ped file in original PlantImpute format, similar to AlphaImpute.")
		("genfile", po::value<string>(&genfilename), "Genotype file in original PlantImpute format, similar to AlphaImpute.")
		("createhapfile", po::value<string>(&outputhapfilename), "Output a hapfile based on input haplotypes.");

	auto parser = po::command_line_parser(argc, argv);
	parser.options(desc);
	po::store(parser.run(), inOptions);
	po::notify(inOptions);

	if (mapfilename != "")
	  {
		FILE* mapfile = fopen(mapfilename.c_str(), "rt");
		cerr << "Reading map file " << mapfilename << "\n";
		readalphamap(mapfile);
	  }

	if (pedfilename != "")
	  {
		FILE* pedfile = fopen(pedfilename.c_str(), "rt");
		cerr << "Reading pedigree file " << pedfilename << "\n";
		readalphaped(pedfile);
	  }

	if (genfilename != "")
	  {
		FILE* genofile = fopen(genfilename.c_str(), "rt");
		cerr << "Reading genotype file " << genfilename << "\n";
		readalphadata(genofile);
	  }

	// TODO: Make sets of required params.
	samplereader samples;
	vector<istream*> hapFiles;
	if (samplefilename != "")
	{
	  std::ifstream sampleFile(samplefilename);
	  samples.read(sampleFile);

	  std::ifstream bimFile(inOptions["bimfile"].as<string>());
	  vector<string> hapsfileOption = inOptions["hapfiles"].as<vector<string>>();

	  for (string filename : hapsfileOption)
	    {
	      hapFiles.push_back(new ifstream(filename));
	    }
	  
	  readhaps(samples.samples, bimFile, hapFiles);
	  std::cout << "readhapssample finished." << std::endl;
	}

	bool docompare = (impoutput != "");
	if (inOptions.count("famfile") + inOptions.count("bedfile") == 2)
	{
	  std::cout << "readfambed started." << std::endl;
	  readfambed(famfilename, bedfilename, docompare);
	}

	// Put generation 2 first, since those are more complex to analyze, avoiding a few threads
	// getting stuck towards the end.
	stable_sort(dous.begin(), dous.end(), [] (individ* a, individ* b) { return a->gen > b->gen; } );
    if (docompare)
	{
		std::ifstream filteredOutput(impoutput);
		compareimputedoutput(filteredOutput);

		return 0;
	}
#endif

	//	return 0;
	CORRECTIONINFERENCE = true;
	postmarkerdata();
	CORRECTIONINFERENCE = false;

	if (deserializefilename != "")
	{
		std::ifstream deserializationFile(deserializefilename);

		std::cout << "deserialize started." << std::endl;
		deserialize(deserializationFile);
		std::cout << "deserialize finished." << std::endl;
	}
	int chromnum;

	if (samplefilename != "" && outputhapfilename != "")
	  {
	    // Note that C++98 had strange EOF behavior (?)
	    hapFiles[0]->seekg(0);
	    std::ofstream outputhapfile(outputhapfilename);
	    createhapfile(samples.samples, *hapFiles[0], outputhapfile);
	    
	    return 0;
	  }

	FILE* out = stdout;
	if (outputfilename != "")
	{
		out = fopen(outputfilename.c_str(), "w");
	}

	if (HAPLOTYPING || true)
		for (int i = 0; i < COUNT; i++)
		{
			//		  	  	{
			early = (i < 1);
			if (!early) doit<false, genotypereporter>((i == COUNT - 1) ? out : stdout, /*i == COUNT - 1*/ true
#ifdef F2MPI
				, world
#endif
				);
			//		}

			//	doit<true>(out);
			//		if (i == COUNT - 1)
			{
				//		    fprintf(out, "\n");
				for (unsigned int j = 0; j < chromstarts[1]; j++)
				{
					for (int t = 0; t < 5; t++)
					{
						//			    fprintf(out, "%0.5lf\t", markerweight[j][t]);
						markerweight[j][t] = 0;
					}
					//			fprintf(out, "\n");		
				}
				//		    fprintf(out, "\n");
			}
			fflush(stdout);
			fflush(out);

			for (unsigned int i2 = 0; i2 < INDCOUNT; i2++)
			{
				individ* ind = getind(i2);
				if (!ind) continue;

				if (ind->haplocount.size())
				{
#ifdef F2MPI
					if (!world.rank())
#endif
						/*if (i == COUNT - 1)*/						fprintf(stdout, "%d %s\n", i2, ind->name.c_str());
					// Printing of haplotype data for each iteration
					for (unsigned int c = 0; c < chromstarts.size() - 1; c++)

					{
						for (unsigned int j = chromstarts[c]; j < chromstarts[c + 1]; j++)
						{

#ifdef F2MPI
							if (!world.rank())
#endif
								/*if (i == COUNT - 1)*/ fprintf(stdout, "%f\t%d\t%d\t\t%f\t%lf %lf %lf\n", ind->haploweight[j], ind->markerdata[j].first.value(), ind->markerdata[j].second.value(), ind->negshift[j],
												ind->markersure[j].first, ind->markersure[j].second, RELSKEWS ? ind->relhaplo[j] : 0);
							ind->negshift[j] = 0;
						}
						//if (!world.rank()) fprintf(out, "\n");
					}
				}
			}
			fflush(stdout);
			fflush(out);
		}

	return 0;
}

