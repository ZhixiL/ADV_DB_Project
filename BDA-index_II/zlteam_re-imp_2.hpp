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
// 
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

