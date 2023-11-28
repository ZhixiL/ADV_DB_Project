#pragma once
#include "stream.h"
#include "uint40.h"
#include <math.h>
// #include "bda-index_II.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
using namespace std;
using namespace sdsl;

/**
 * We have grouped all the re-implementation we did in this single file. However, they are spreaded on
 * Index-I - [bda.cc, bda-index_I.cc, pattern-matching.cc]
 * Index-2 - [bda.cc, bda-index_II.cc, pattern-matching.cc]
 * Such that they can replace the original implementation.
 */

/* ------ FOLLOWING ARE DEDICATED FOR STEP 1: BDA-COMPUTE ------ */ 
// This is the function necessary for the pattern matching that occured during BDA-Compute on which 
// it is able to compute the longest common prefix.
INT lcp_zlteam ( string & x, INT x_begin, string & y, INT y_begin )
{
    // setup boundaries for how many itrs we can test for.
    INT bound = std::min(x.size()-x_begin, y.size()-y_begin);
    // handle corner cases when one of the param is outside of bound.
    if (bound <= 0) return 0;
    // now check from x_begin & y_begin up until the bound to see how many in common.
    for (i = 0; i < bound; i++) {
        if (x[x_begin+i] != y[y_begin+i]) {
            return i;
        }
    }
    // all prefixes in x matched with y, we have the longest matching.
    return bound-1;
}

// Similar to lcp, simple algorithm to test for the longest common suffix between the two strings.
INT lcs_zlteam(string & x, INT x_end, string & y, INT y_start)
{
    // setup bound for how many char we test at most
    INT bound = std::min(x_end + 1 || y.size() - y_start);
    // check for corner cases where we can't match at all;
    if (bound <= 0) return 0;
    // Now we iteratively checks for longest common suffix, 
    for (INT i = 0; i < bound; i++) {
        if (x[x_end-i] != y[y_start+i]) {
            return i
        }
    }
    // all suffixes in x matched with y, we have the longest matching.
    return bound-1;
}

/* ------ FOLLOWING ARE DEDICATED FOR STEP 1: BDA-COMPUTE ------ */ 

INT red_minlexrot_zlteam( string &s, INT *f, INT n, INT r){           //find lexicographically minimum rotation, return startPos which is the starting position of the LMR
    std::string ss = s + s;     //append the same string together to achieve something similar to a circular shift      
    int i = 0, ans = 0;
    
    while (i < n) {
        ans = i;
        int j = i + 1, k = i;
         
        while (j < n + i) {     
            if (ss[k] > ss[j]) {        //break the loop if smaller character is found
                break;
            }
            if (ss[k] < ss[j]) {        //otherwise update k based on the comparison of characters and increment j to move on to the next character
                k = i;
            } else {
                k++;
            }
            j++;
        }
        
        while (i <= k) {        //update i for the next iteration
            i += j - k;
        }
    }
    
    return ans;
}

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif
//#include "bda.cc" // ensure we have access to LCParray function.

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
        rank_count += (LCP[j] < k);
        rank[SA[j]] = rank_count;
    }
    INT fragLength  = w - k - 1;
    //iterate through the first fragment with fragment size being w-k-1 (w,k) minimizer length
    for (INT i = 0; i < fragLength; i++){
	//remove elements from the back if the current rank is smaller so its in ascending order
        while (!minimizer_rankings.empty() && rank[i] < minimizer_rankings.back().first){
            minimizer_rankings.pop_back();
        }
	//create a potential bd and add it to the minimizer rankings
        utils::Rank bd{.start_pos = i, .rank_pos = i};
        minimizer_rankings.emplace_back(rank[i], std::move(bd));
    }

    //iterate through the remaining fragments and perform the same actions
    for (INT j = 0; j <= n - w; j++){
        while (!minimizer_rankings.empty() && rank[fragLength] < minimizer_rankings.back().first){
            minimizer_rankings.pop_back();
        }
        utils::Rank bd{.start_pos = fragLength, .rank_pos = fragLength};
        minimizer_rankings.emplace_back(rank[fragLength], std::move(bd));
        
	//remove elements from the front of minimizer rankings if its outside of the current fragment
 	while (!minimizer_rankings.empty() && minimizer_rankings.front().second.start_pos <= fragLength - w + k)
            minimizer_rankings.pop_front();
	
	//populate the minimizer deque with potential bd anchors with minimum rank
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
		//determine the distance to the end of the fragment, select the smallest distance between current position to the end of the fragment from all minimizers 
                INT dist_end = min(j + w - max(minimizers[i].rank_pos, minimizers[bestCandidate].rank_pos), w);
                INT h1 = 0;
		//get the start position of the current minimizer (curPos) and the smallest minimizer (bestPos)
                INT curPos = minimizers[i].rank_pos;            
                INT bestPos = minimizers[bestCandidate].rank_pos;
		//simply find LCP of F[i..|F|] and F[j..|F|]
                for (INT k = 0; k < n; k++){
                    if (seq[curPos + k] == seq[bestPos + k]){
                        h1++;
                    }
                    else{
                        break;
                    }
                }
		//first comparison makes sure the length of LCP isnt larger than the distance to the end of the fragment, and second comparison makes sure that the current minimizer has larger inverse SA value
                if (h1 < dist_end && invSA[bestPos] > invSA[curPos]){
                    bestCandidate = i;
                }
                else{
                    //second query: find h2 = LCP of F[i+h1..|F|] and F, if h2 < j - i then we compare F[i+h1+h2] and F[1 + h2] otherwise 
                    curPos = j;
		    //update minimizer rank position 
                    bestPos += min(h1, dist_end);
		    //calculate the new distance to end between current pos and minimum pos
                    dist_end = min(j + w - max(curPos, bestPos), w);
                    INT h2 = 0;

		    //find length of LCP which is h2
                    for (INT k = 0; k < n; k++){
                        if (seq[curPos + k] == seq[bestPos + k]){
                            h2++;
                        }
                        else{
                            break;
                        }
                    }
		    //similar comparison as previous
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
void build_RSA_RLCP_zlteam (unordered_set<INT> &anchors, INT n, INT * RSA, INT * RLCP, INT ram_use, char * sa_fname, char * lcp_fname) {
    /* Since the SA and LCP are residing in EM (disk) now, we read them and then store them to memory */
    // As this portion just read off the disk file to obtain SA and LCP, with no technical merit,
    // we will not re-implement but rather reuse the original construct.
	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);

    // Following code are implemented based on the paper description:
    deque<INT> SA_deq;
    INT curLCP = 0;
    // Attempted to optimize by reading the SA and LCP at once to save read operations, but this breaks
    // the purpose that SA/LCP will have to take spaces on memory again, on which the author did sequential
    // read to make sure they only take O(1) time, so this way is removed.
	// deque<INT> LCP_deq;
    // for (INT i = 0; i <=n; ++i) {
    //     SA_deq.push_back(SA->read());
    //     LCP_deq.push_back(LCP->read());	
    // }
    
    // As our RLCP will be different from LCP, we need a value that keep track of the actual previous 
    // LCP between i and i-1 when we meet discontinuity, so we make use of this variable.
    INT recordLCP = n;
	INT RPos = 0; // trackers where we at for RSA & RLCP as they will be shorter than SA/LCP
    // set the default value for RLCP's first element, because there's no previous index to compare with.
    RLCP[0] = 0; 
    // Now we carry out sampling with SA and LCP against hashtable (anchors) from bda_compute,
    // Such that we can build the RSA and RLCP array.
    for (INT i = 0; i <= n; i++) {
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
            RSA[RPos] = SA_deq.back();
            // "It is also well known that for any i1 < i2 the length of the LCP between S[SA[i1] .. n] 
            // and S[SA[i2] .. n] is the minimum value lying in LCP[i1 + 1], . . . , LCP[i2]." 
            if (RPos > 0) { // Since we want to compare past of SA and RSA, we can't use 0.
                // "Since we scan also the LCP array simultaneously, we maintain the value we 
                // need to store in RLCP[i]."
                if (RSA[RPos-1] == SA_deq.front()) { 
                    RLCP[RPos] = curLCP; // if prev SA is in anchor too, the following LCP will be the same.
                } else {
                    // But if they are not continuous, we need to use the recorded LCP 
                    // We used comparison instead of std::min as we believe this is faster.
                    RLCP[RPos] = recordLCP < curLCP ? recordLCP : curLCP;
                }
            }
            recordLCP = n; // as recordLCP may be modified later, we modify reset it.
			RPos++;
        } else {
            // memorize the current longest common prefix length.
            recordLCP = recordLCP < curLCP ? recordLCP : curLCP; 
        }
    }

    // garbage collection.
    delete(SA);
	delete(LCP);
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
	INT LPos = 0;

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
            LSA[LPos] = SA_deq.back();
            // "It is also well known that for any i1 < i2 the length of the LCP between S[SA[i1] .. n] 
            // and S[SA[i2] .. n] is the minimum value lying in LCP[i1 + 1], . . . , LCP[i2]." 
            if (i > 0) { // Since we want to compare past of SA and LSA, we can't use 0.
                // "Since we scan also the LCP array simultaneously, we maintain the value we 
                // need to store in RLCP[i]." Replace with LLCP.
                if (LSA[LPos-1] == SA_deq.front()) { 
                    LLCP[LPos] = curLCP; // if prev SA is in anchor too, the following LCP will be the same.
                } else {
                    // But if they are not continuous, we need to use the recorded LCP 
                    // We used comparison instead of std::min as we believe this is faster.
                    LLCP[LPos] = recordLCP < curLCP ? recordLCP : curLCP;
                }
            }
            recordLCP = n; // as recordLCP may be modified later, we modify reset it.
			LPos++;
        } else {
            // memorize the current longest common prefix length.
            recordLCP = recordLCP < curLCP ? recordLCP : curLCP; 
        }
    }

    // garbage collection.
    delete(LCP);
    delete(SA);
}

