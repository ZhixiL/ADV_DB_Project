#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_II.h"
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
//#include "zlteam_re-imp_2.hpp"

using namespace std;
using namespace sdsl;

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif

/* Kasai et al algorithm for O(n)-time LCP construction */
INT LCParray ( unsigned char * text, INT n, INT * SA, INT * ISA, INT * LCP )
{
        INT i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }
        return ( 1 );
}

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


/* Computes the bd-anchors of a string of length n in O(n) time */
INT bd_anchors(  unsigned char * seq, INT pos, INT ell, INT k, unordered_set<INT> &anchors, INT * SA, INT * LCP, INT * invSA, INT * rank  )
{
 	return bd_anchors_zlteam(seq, pos, ell, k, anchors, SA, LCP, invSA, rank);

	// INT w = ell;
	// INT n = strlen ( (char*) seq );
	
		
	// /* Compute suffix array for block */
	// #ifdef _USE_64
  	// if( divsufsort64( seq, SA,  n ) != 0 )
  	// {
  	// 	fprintf(stderr, " Error: SA computation failed.\n" );
    //       	exit( EXIT_FAILURE );
  	// }
	// #endif

	// #ifdef _USE_32
  	// if( divsufsort( seq, SA,  n ) != 0 )
  	// {
  	// 	fprintf(stderr, " Error: SA computation failed.\n" );
    //       	exit( EXIT_FAILURE );
  	// }
	// #endif
	
	
	

	// for ( INT i = 0; i < n; i ++ )
	// {
	//         invSA [SA[i]] = i;
	        
	// }

	
	// //std::chrono::steady_clock::time_point  start_lcp = std::chrono::steady_clock::now();
	
	// /* Compute the LCP array for block */
	// if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
	// {
	// 	fprintf(stderr, " Error: LCP computation failed.\n" );
	// 	exit( EXIT_FAILURE );
	// }


	// INT rank_count =  0;	
	// rank[SA[0]] = rank_count;
		
	// /* Compute the ranks for block */
	// for(INT j =1; j<n; j++)
	// {
	// 	if( LCP[j] >= k )
	// 	{
	// 		rank[SA[j]] = rank_count;
	// 	}
	// 	else 
	// 	{
	// 		rank_count++;
	// 		rank[SA[j]] = rank_count;
	// 	}
	// }	
	
	
	// deque<pair<INT,utils::Rank>> min_rank;
	// vector<utils::Rank> minimizers;
		
	// INT j;
   	// for (INT j = 0; j < w - k - 1; j++) 
   	// {
 	// 	while ( !min_rank.empty() && rank[j] < min_rank.back().first )
 	// 		min_rank.pop_back();
 
    //    		utils::Rank potential_bd;
	// 	potential_bd.start_pos = j;
	// 	potential_bd.rank_pos = j;
				
	// 	min_rank.push_back(std::make_pair(rank[j], potential_bd));
		
    // 	}
    	
	// /* Compute reduced bd-anchors for every window of size ell */
	
	// INT i = w - k - 1;
	// for( j = 0; j<=n-w; j++ )
	// {
		
	// 	while (!min_rank.empty() && min_rank.back().first > rank[i])
	// 		min_rank.pop_back();
					
	// 	utils::Rank potential_bd;
	// 	potential_bd.start_pos = i;
	// 	potential_bd.rank_pos = i;
				
	// 	min_rank.push_back(std::make_pair(rank[i], potential_bd));
		
	
	// 	while( min_rank.front().second.start_pos <= i - w + k)
	// 	{
	// 		min_rank.pop_front();
	// 	}	
		

	// 	INT min_ = min_rank.at(0).first;
	// 	for(INT i = 0; i<min_rank.size(); i++)
	// 	{
	// 		if( min_rank.at(i).first == min_ )
	// 		{
	// 			minimizers.push_back( min_rank.at(i).second );
	// 		}
	// 		else if( min_rank.at(i).first >  min_ )
	// 			break;
	// 	}
		
	// 	i++;
	
	// 	/* Filter draws if there are more than one minimum rank, otherwise 
	// 	only one potential bd-anchor in window */			
	// 	if( minimizers.size() > 1 )
	// 	{ 	
			
	// 		INT minimum = 0;
			
	// 		for(INT i = 1; i<minimizers.size(); i++)
	// 		{
		
	// 			INT dist_to_end = w;
				
	// 			INT rank_pos = minimizers.at(i).rank_pos;
	// 			INT min_rank_pos = minimizers.at(minimum).rank_pos;
				
	// 			if( ( (j+ w ) - rank_pos ) < dist_to_end )
	// 				dist_to_end = ( (j+ w ) - rank_pos );
				
	// 			if( ( (j+ w ) -  min_rank_pos ) < dist_to_end )
	// 				dist_to_end = ( (j+ w ) -  min_rank_pos );
				
	// 			INT min_inv = min( invSA[min_rank_pos], invSA[rank_pos])+1 ;
				
					
	// 			INT max_inv = max( invSA[min_rank_pos], invSA[rank_pos]) ;
				
	// 			INT lcp1 = 0; 
				
	// 			while ( seq[min_rank_pos+lcp1] == seq[rank_pos+lcp1] )
	// 				lcp1++;
			
				
	// 			if( lcp1 < dist_to_end )
	// 			{
					
	// 				if( invSA[ min_rank_pos ] > invSA[ rank_pos ] )
	// 				{
	// 					minimum = i;
	// 				}
	// 			}
	// 			else
	// 			{
					
					
	// 				min_rank_pos =  min_rank_pos + min(lcp1,dist_to_end) ;
	// 				rank_pos = j;
				
	// 				dist_to_end = w;
	// 				if( ( (j+ w ) - rank_pos ) < dist_to_end )
	// 					dist_to_end = ( (j+ w ) - rank_pos );
				
	// 				if( ( (j+ w ) -  min_rank_pos ) < dist_to_end )
	// 					dist_to_end = ( (j+ w ) -  min_rank_pos );
				
	// 				INT min_inv = min( invSA[min_rank_pos], invSA[rank_pos])+1 ;
					
						
	// 				INT lcp2 = 0; 
				
	// 				while ( seq[min_rank_pos+lcp2] == seq[rank_pos+lcp2] )
	// 					lcp2++;
				
	// 				if( lcp2 < dist_to_end )
	// 				{
						
	// 					if( invSA[ min_rank_pos ] > invSA[ rank_pos ] )
	// 					{
	// 						minimum = i;
	// 					}
	// 				}
	// 				else
	// 				{
						
						
						
	// 				 	min_rank_pos = min_rank_pos + min(lcp2,dist_to_end);
	// 					rank_pos = rank_pos + min(lcp2,dist_to_end);
	// 					dist_to_end = w;
	// 					if( ( (minimizers.at(i).start_pos) - rank_pos ) < dist_to_end )
	// 						dist_to_end = ( (minimizers.at(i).start_pos) - rank_pos );
						
	// 					if( ( (minimizers.at(i).start_pos) -  min_rank_pos ) < dist_to_end )
	// 						dist_to_end = ( (minimizers.at(i).start_pos) -  min_rank_pos );
						
	// 					INT min_inv = min( invSA[min_rank_pos], invSA[rank_pos])+1 ;
	// 					INT max_inv = max( invSA[min_rank_pos], invSA[rank_pos]) ;
						
	// 					INT lcp3 = 0; 
						
	// 					while ( seq[min_rank_pos+lcp3] == seq[rank_pos+lcp3] )
	// 						lcp3++;
						
	// 					if( lcp3 < dist_to_end )
	// 					{
	// 						if( invSA[ min_rank_pos ] > invSA[ rank_pos ] )
	// 						{
	// 							minimum = i;
	// 						}
	// 					}
	// 				}
					
					
				
	// 			}
	// 		}
	// 		anchors.insert( minimizers.at(minimum).start_pos+pos );
		

	// 	}
	// 	else 
	// 	{
	// 		anchors.insert( minimizers.at(0).start_pos+pos );
	// 	}
			
		
	// 	minimizers.clear();
	// }
						
	return 0;
}
