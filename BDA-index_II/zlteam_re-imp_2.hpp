#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_II.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
using namespace std;
using namespace sdsl;

// Initial reimplmentation of LCP core algo to get us introduced to the implemenation format of this project
// and let us become familiar with the datatype that the authors have used.
// This function calculates the longest common prefix between the two input strings
// s1end scans from s1end to left all the way to 0
// s2start scans from starts right all the way to the length of s2.
INT lcp_zlteam ( string & s1, INT s1end, string & s2, INT s2start )
{
    // clear out the edge case where s1end as well as s2start is not in range.
    if (s1end < 0 || s2start >= s2.size())
        return 0;

    // Simply match one by one to find the longest common prefix.
    INT init_s2start = s2start;
    while (s1end >= 0 && s2start < s2.size()) {
        if (s1[s1end] == s2[s2start]) {
            --s1end;
            ++s2start;
        } else {
            return s2start - init_s2start;
        }
    }
    return 0;
}

INT red_minlexrot_zlteam( string &X, INT strLength){           //find lexicographically minimum rotation, return startPos which is the starting position of the LMR
    INT startPos = 0;
    string curMax = x;
    for (INT i = 0; i < n; i++){
        x = x.substr(1, n-1) + x[0];            //perform rotation
        if (curMax > x){                    //find lexicographically minimum rotation
            startPos = i+1;
            curMax = x;
        }
    }
    return startPos;
}



#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif
#include "bda.cc" // ensure we have access to LCParray function.

// This is the re-implemented bd_anchor compute function that is mentioned in paper section 4.1 and 
// partly 4.2. We have retained all the parameters to ensure that this function can be properly called
// by the original construct of the program. We wrote this function based on the thorems, images, as well
// as the original code.
// To further demonstrate our understanding we will write out the comments for each of the sections as we
// implement.
INT bd_anchors_zlteam(unsigned char * seq, INT pos, INT w, INT k, 
                      unordered_set<INT> &anchors, INT * SA, INT * LCP, INT * invSA, INT * rank)
{
    // Initial setup
    // Suffix Array computation with sorting algorithm from divsufsort.h, thus
    // the implementation is identical to the orignal paper Lorraine et al.
    INT n = strlen ( (char*) seq );
	#ifdef _USE_64
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif

	#ifdef _USE_32
  	if( divsufsort( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif

    // (w, k)-minimizer setup


    // utilizing (w, k)-minimizer instead of anchors (as it is a superset for reduced bd-anchors)
    // such that it can be used to compute reduced l-most bd anchors. 

    // Used LCP to compare anchored rotations for fragment of input string S
}