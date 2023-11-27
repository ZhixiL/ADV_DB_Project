#ifndef ZLTEAM
#define ZLTEAM

#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_II.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
using namespace std;
using namespace sdsl;

/* ------ FOLLOWING ARE DEDICATED FOR STEP 1: BDA-COMPUTE ------ */ 

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

// This is the re-implemented bd_anchor compute function that is mentioned in paper 4.1, and partly 4.2
// on step 1. We have retained all the parameters to ensure that this function can be properly called
// by the original construct of the program. We wrote this function based on the thorems, images, as well
// as the original code. Certain logics we have opted to do it the original code's way if we found out that
// it may impact the potential computation complexity by a huge margin. Others like the initialization that
// Lorraine et al. importd from other people's work we will keep the same thing as they will not be the core
// objective of this project.
// To further demonstrate our understanding we will write out the comments for each of the sections as we
// implement.
INT bd_anchors_zlteam(unsigned char * seq, INT pos, INT w, INT k, unordered_set<INT> &anchors, INT * SA, INT * LCP, INT * invSA, INT * rank)
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
            anchors.insert(minimizers[bestCandidate].start_pos + pos);      //insert the result minimizer into BD-Anchor
        }
        else{
            anchors.insert(minimizers[0].start_pos + pos);
        }
    minimizers.clear();
    }
    return 0;
}


/* We ignore step 2 as the step 2 only involves in using external memories to build SA & LCP array */

/* ------ FOLLOWING ARE DEDICATED FOR STEP 3: Left & Right tree building preparation ------ */
/* Construct a right tree provided with RSA and RLCP, generated by the original input string. */
// One parameter change we made is remove the parameter g, because we believe it is irrelevant.
// *As the paper gave a very good explaination on how to implement this algorithm, we will quote with ""
// what paper described and our corresponding implementation. Our implementation will be very similar to
// the original implementation due to the code logic is not that complex, so to show our understanding
// we will attempt to make changes that will hopefully make the code more efficient, and comment as much
// as we can.
void build_RSA_RLCP_zlteam ( unordered_set<INT> &anchors, INT n, INT * RSA, INT * RLCP, INT ram_use, char * sa_fname, char * lcp_fname ) {
    /* Since the SA and LCP are residing in EM (disk) now, we read them and then store them to memory */
    // As this portion just read off the disk file to obtain SA and LCP, with no technical merit,
    // we will not re-implement but rather reuse the original construct.
	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);

    // Following code are implemented based on the paper description:
    deque<INT> SA_deq;
    // deque<INT> LCP_deq;
    INT curLCP = 0;
    // Attempted to optimize by reading the SA and LCP at once to save read operations, but this breaks
    // the purpose that SA/LCP will have to take spaces on memory again, on which the author did sequential
    // read to make sure they only take O(1) time, so this way is removed.
    // for (INT i = 0; i <=n; ++i) {
    //     SA_deq.push_back(SA->read());
    //     LCP_deq.push_back(LCP->read());
    // }
    
    // As our RLCP will be different from LCP, we need a value that keep track of the actual previous 
    // LCP between i and i-1 when we meet discontinuity, so we make use of this variable.
    INT recordLCP = n;
    // set the default value for RLCP's first element, because there's no previous index to compare with.
    RLCP[0] = 0; 
    // Now we carry out sampling with SA and LCP against hashtable (anchors) from bda_compute,
    // Such that we can build the RSA and RLCP array.
    for (INT i = 0; i <= n; ++i) {
        // sequentially read SA and LCP
        SA_deq.push_back(SA->read());
        curLCP = LCP->read();
        // make sure that we are only keep track at most 2 SA at a time such that
        // they only take O(1) memory space, on which we don't need to store the entire SA/LCP.
        if (SA_deq.size() > 2) SA_deq.pop_front();

        // "find if SA[i] exist in our anchors hash table"
        if (anchors.find(SA_deq.back()) != anchors.end())
        {
            // "if so, we set RSA[i] = SA[i]."
            RSA[i] = SA_deq.back();
            // "It is also well known that for any i1 < i2 the length of the LCP between S[SA[i1] .. n] 
            // and S[SA[i2] .. n] is the minimum value lying in LCP[i1 + 1], . . . , LCP[i2]." 
            if (i > 0) { // Since we want to compare past of SA and RSA, we can't use 0.
                // "Since we scan also the LCP array simultaneously, we maintain the value we 
                // need to store in RLCP[i]."
                if (RSA[i-1] == SA_deq.front()) { 
                    RLCP[i] = curLCP; // if prev SA is in anchor too, the following LCP will be the same.
                } else {
                    // But if they are not continuous, we need to use the recorded LCP 
                    // We used comparison instead of std::min as we believe this is faster.
                    RLCP[i] = recordLCP < curLCP ? recordLCP : curLCP;
                }
            }
            recordLCP = n; // as recordLCP may be modified later, we modify reset it.
        } else {
            // memorize the current longest common prefix length.
            recordLCP = recordLCP < curLCP ? recordLCP : curLCP; 
        }
    }

    // garbage collection.
    delete(LCP);
    delete(SA);
}

// For building the LSA and LLCP, we build it with reverse version of the original input, on which the only
// difference is how we access the elements in anchors, as the SA and LCP that we got is already reversed
// from author's other portion of the code in main().
void build_LSA_LLCP_zlteam ( unordered_set<INT> &anchors, INT n, INT * LSA, INT * LLCP, INT ram_use, char * sa_fname, char * lcp_fname )
{
    // Same initial setup as the right direction construction
	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);
    deque<INT> SA_deq;
    INT curLCP = 0;
    INT recordLCP = n;
    // set the default value for LLCP's first element, because there's no previous index to compare with.
    LLCP[0] = 0; 

    // Now we carry out sampling with SA and LCP against hashtable (anchors) from bda_compute,
    // Such that we can build the LSA and LLCP array.
    for (INT i = 0; i <= n; ++i) {
        // sequentially read SA and LCP
        SA_deq.push_back(SA->read());
        curLCP = LCP->read();
        // make sure that we are only keep track at most 2 SA at a time such that
        // they only take O(1) memory space, on which we don't need to store the entire SA/LCP.
        if (SA_deq.size() > 2) SA_deq.pop_front();

        /**
         * The only difference from the right tree construction is the query part on which we take the 
         * inverse of the current SA[i] such that we are utilizing anchors that is from the reverse
         * direction. We use what paper mentioned of (n-1)-SA[i] to obtain that inverse index.
         */
        if (anchors.find((n-1)-SA_deq.back()) != anchors.end())
        {
            // "if so, we set RSA[i] = SA[i]." replace with LSA.
            LSA[i] = SA_deq.back();
            // "It is also well known that for any i1 < i2 the length of the LCP between S[SA[i1] .. n] 
            // and S[SA[i2] .. n] is the minimum value lying in LCP[i1 + 1], . . . , LCP[i2]." 
            if (i > 0) { // Since we want to compare past of SA and LSA, we can't use 0.
                // "Since we scan also the LCP array simultaneously, we maintain the value we 
                // need to store in RLCP[i]." Replace with LLCP.
                if (LSA[i-1] == SA_deq.front()) { 
                    LLCP[i] = curLCP; // if prev SA is in anchor too, the following LCP will be the same.
                } else {
                    // But if they are not continuous, we need to use the recorded LCP 
                    // We used comparison instead of std::min as we believe this is faster.
                    LLCP[i] = recordLCP < curLCP ? recordLCP : curLCP;
                }
            }
            recordLCP = n; // as recordLCP may be modified later, we modify reset it.
        } else {
            // memorize the current longest common prefix length.
            recordLCP = recordLCP < curLCP ? recordLCP : curLCP; 
        }
    }

    // garbage collection.
    delete(LCP);
    delete(SA);
}


#endif