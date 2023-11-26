#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_II.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
using namespace std;
using namespace sdsl;

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

INT bd_anchors_zlteam(unsigned char * seq, INT pos, INT w, INT k, 
                      unordered_set<INT> &anchors, INT * SA, INT * LCP, INT * invSA, INT * rank)
{
    // Initial setup 
    // --- IMPORTED FROM ORIGINAL CODE ---
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
    // a common way to invert suffix array that tells you the position 
    // of each suffix in the sorted list of suffixes
	for ( INT i = 0; i < n; i ++ ) 
        invSA [SA[i]] = i;
    // As the LCParray function is not part of this paper's contribution, we will not 
    // reimplement it. This portion essentially setup the empty input array of LCP to
    // containing the longest common prefixes. 
	if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
		exit( EXIT_FAILURE );
	}
    // --- END OF IMPORTING ORIGINAL CODE ---

    // (w, k)-minimizer setup
    // minimizer_rankings is using the same rank struct as provided and the deque to ensure random
    // access time. We didn't use other data type as that will greatly affect the 
    // thoretical computation time.
    deque<pair<INT,utils::Rank>> minimizer_rankings;
    // minimizer that acts as a superset to bd-anchors, as its easier to compute minimizer than 
    // bd-anchors due to shorter comparisons with window size w.
    // For data structure, we opted to use deque instead of the vector as it will be more efficient for
    // push_back even if it's slightly, and provides the same O(1) time for accessing elements.
    deque<utils::Rank> minimizers;

    // utilizing (w, k)-minimizer instead of anchors (as it is a superset for reduced bd-anchors)
    // such that it can be used to compute reduced l-most bd anchors. 

    //reused code from original paper to calculate rank arrays
    INT rank_count = 0;
    rank[SA[0]] = rank_count;

    for (INT j = 1; j < n; j++){
        rank_count += (LCP[j] < k)
        rank[SA[j]] = rank_count;
    }
    INT fragLength  = w - k - 1;

    //start by sorting the letters of F and assigning it a rank in {1 ... 2l} accordingly, compute reduced BD-Anchors and then use three LCP queries to identify lexi smallest rotation

    for (INT i = 0; i < fragLength; i++){
        while (!minimizer_rankings.empty() && rank[i] < minimizer_rankings.back().first){
            minimizer_rankings.pop_back()
        }
        utils::Rank bd{.start_pos = i, .rank_pos = i};
        minimizer_rankings.emplace_back(rank[i], std::move(bd));
    }

    for (INT j = 0; j <= n - w; j++){
        while (!minimizer_rankings.empty() && rank[i] < minimizer_rankings.back().first){
            minimizer_rankings.pop_back()
        }
        utils::Rank bd{.start_pos = i, .rank_pos = i};
        minimizer_rankings.emplace_back(rank[i], std::move(bd));
        while (!minimizer_rankings.empty() && minimizer_rankings.front().second.start_pos <= fragLengh - w + k)
            minimizer_rankings.pop_front();

        if (!minimizer_rankings.empty())
        {
            INT minVal = minimizer_rankings.front().first;
            for (const auto &elem : minimizer_rankings)
            {
                if (elem.first == minVal)
                    minimizers.push_back(elem.second);
                else if (elem.first > minVal)
                    break;
            }
        }
        fragLength++;

        //detect if theres more than one minimal rotation in each fragment, if so, we employ three LCP queries and three letter comparisons mentioned in the paper to determine the best candidate for minimizer
        if (minimizers.size() > 1){
            
            //first query: find h1 = LCP of F[i .. |F|] and F[j .. |F|]. If h1 < |F| - j + 1 (dist_end), we will compare F[i+h1] and F[j+h1] for the answer otherwise queue query 2
            INT bestCandidate = 0;
            for (INT i = 1; i < minimizers.size(); i++){
                INT dist_end = min(j + w - max(minimizer[i].rank_pos, minimizer[bestCandidate].rank_pos), w);
                INT h1 = 0;
                INT curPos = minimizers[i].rank_pos;
                INT bestPos = minimizers[bestCandidate].rank_pos;
                for (INT k = 0; k < n; k++){
                    if (seq[curPos + k] == seq[bestPos + k]){
                        h1++;
                    }
                    else{
                        break;
                    }
                }
                if (h1 < dist_end && invSA[bestPos] > invSA[curPos]){
                    bestCandidate = i;
                }
                else{
                    //second query: find h2 = LCP of F[i+h1..|F|] and F, if h2 < j - i then we compare F[i+h1+h2] and F[1 + h2] otherwise 
                    curPos = j;
                    bestPos += min(h1, dist_end);
                    dist_end = min(j + w - max(curPos, bestPos), w);
                    INT h2 = 0;

                    for (INT k = 0; k < n; k++){
                        if (seq[curPos + k] == seq[bestPos + k]){
                            h2++;
                        }
                        else{
                            break;
                        }
                    }
                    if (h2 < dist_end && invSA[bestPos] > invSA[curPos]){
                        bestCandidate = i;
                    }
                    else{
                        //third query: find h3 = LCP of F and F[j - i + 1 .. |F|]. If h3 < j - i we compare F[1 + h3] and F[j - i + 1 + h3]
                        curPos += min(h2, dist_end);
                        bestPos += min(h2, dist_end);
                        dist_end = min(minimizers[i].start_pos - max(curPos, bestPos), w);
                        INT h3 = 0;
                        for (INT k = 0; k < n; k++){
                            if (seq[curPos + k] == seq[bestPos + k]){
                                h3++;
                            }
                            else{
                                break;
                            }
                        }
                        if (h3 < dist_end && invSA[bestPos] > invSA[curPos]){
                            bestCandidate = i;
                        }
                }
            }
        }
            anchors.insert(minimizers[bestCandidate].start_pos + pos);
        }
        else{
            anchors.insert(minimizers[0].start_pos + pos);
        }
    minimizers.clear();
    }
    return 0;
}