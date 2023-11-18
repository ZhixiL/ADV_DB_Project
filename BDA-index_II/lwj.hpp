#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_II.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
using namespace std;
using namespace sdsl;
//aabc
//abca
//bcaa
//caab
//aabc
INT red_minlexrot( string &X, INT strLength){           //find lexicographically minimum rotation, return startPos which is the starting position of the LMR
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
