#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "grid.h"
#include "bda-index_I.h"
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>

using namespace std;
using namespace sdsl;

/* ------ FOLLOWING ARE DEDICATED FOR STEP 0: Helpers ------ */ 
// This is the function necessary for the pattern matching that occured during BDA-Compute on which 
// it is able to compute the longest common prefix.
INT lcp_zlteam ( string & x, INT x_begin, string & y, INT y_begin )
{
    // setup boundaries for how many itrs we can test for.
    INT bound = std::min(x.size()-x_begin, y.size()-y_begin);
    // handle corner cases when one of the param is outside of bound.
    if (bound <= 0) return 0;
    // now check from x_begin & y_begin up until the bound to see how many in common.
    for (INT i = 0; i < bound; i++) {
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
    INT bound = std::min(x_end + 1, (INT)y.size() - y_start);
    // check for corner cases where we can't match at all;
    if (bound <= 0) return 0;
    // Now we iteratively checks for longest common suffix, 
    for (INT i = 0; i < bound; i++) {
        if (x[x_end-i] != y[y_start+i]) {
            return i;
        }
    }
    // all suffixes in x matched with y, we have the longest matching.
    return bound-1;
}

INT red_minlexrot_zlteam( string &s, INT *f, INT n, INT r){		//same as minlexrot, but with reduce parameter R so we know what to not consider
	std::string ss = s + s;     //append the same string together to achieve something similar to a circular shift      
    int i = 0, ans = 0;
    
    while (i < n) {
        ans = i;
        int j = i + 1, k = i;
         
        while (j < n + i) {     
            if (k < n - r && ss[k] > ss[j]) {        //break the loop if smaller character is found, consider parameter R for comparison to ignore the trailing R rotations
                break;
            }
            if (k < n - r && ss[k] < ss[j]) {        //otherwise update k based on the comparison of characters and increment j to move on to the next character
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

INT minlexrot_zlteam( string &s, INT *f, INT n){           //find lexicographically minimum rotation, return startPos which is the starting position of the LMR
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


INT red_minlexrot( string &X, INT *f, INT n, INT r )
{  
	return red_minlexrot_zlteam(X, f, n, r);

	// INT n_d = n<<1;
  	// for(INT i = 0; i < n_d; ++i)	f[i] = (INT) -1;

  	// INT k = 0;
  	// for (INT j = 1; j < n_d; ++j)
  	// {
    //             unsigned char sj = X[j%n];
    //             INT i = f[j - k - 1];
    //             while (i != (INT)-1 && sj != X[(k + i + 1)%n])
    //             {
    //                     if (sj < X[(k + i + 1)%n] && j - i - 1 < n - r )        k = j - i - 1;
    //                     i = f[i];
    //             }
				
    //             if (i == (INT) - 1 && sj != X[(k + i + 1)%n])
    //             {
    //                     if (sj < X[(k+i+1)%n] && j - i - 1 < n - r )    k = j;
    //                     f[j - k] = -1;
    //             }
    //             else
    //                     f[j - k] = i + 1;
   	// }
   	// return k;
}


/* Booth's O(n)-time algorithm -- slightly adapted for efficiency */
INT minlexrot( string &X, INT *f, INT n)
{  
	return minlexrot_zlteam(X, f, n);
	// INT n_d = n<<1;
  	// for(INT i = 0; i < n_d; ++i)	f[i] = (INT) -1;

  	// INT k = 0;
  	// for (INT j = 1; j < n_d; ++j)
  	// {
    //             unsigned char sj = X[j%n];
    //             INT i = f[j - k - 1];
    //             while (i != (INT)-1 && sj != X[(k + i + 1)%n])
    //             {
    //                     if (sj < X[(k + i + 1)%n])        k = j - i - 1;
    //                     i = f[i];
    //             }
				
    //             if (i == (INT) - 1 && sj != X[(k + i + 1)%n])
    //             {
    //                     if (sj < X[(k+i+1)%n])    k = j;
    //                     f[j - k] = -1;
    //             }
    //             else
    //                     f[j - k] = i + 1;
   	// }
   	// return k;
}





/* Computes the length of lcp of two suffixes of two strings */
INT lcp ( string & x, INT M, string & y, INT l )
{
	// Replace with the rewritten funciton.
	return lcp_zlteam(x, M, y, l);
	// INT xx = x.size();
	// if ( M >= xx ) return 0;
	// INT yy = y.size();
	// if ( l >= yy ) return 0;

	// INT i = 0;
	// while ( ( M + i < xx ) && ( l + i < yy ) )
	// {
	// 	if ( x[M+i] != y[l+i] )	break;
	// 	i++;
	// }
	// return i;
}

/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<INT,INT> pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n )
{

	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;

	pair<INT,INT> interval;

	while ( d + 1 < f )
	{
		
		INT i = (d + f)/2;
		
		/* lcp(i,f) */
		INT lcpif;
		
		if( f == n )
			lcpif = 0;
		else lcpif = LCP[rmq ( i + 1, f ) ];
			
		/* lcp(d,i) */
		INT lcpdi;
		
		if( i == n )
			lcpdi = 0;
		else lcpdi = LCP[rmq ( d + 1, i ) ];
	
		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			l = l + lcp ( a, SA[i] + l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;

					/* lcp(j,e) */
					INT lcpje;
					
					
					if( e == n )
						lcpje = 0;
					else lcpje = LCP[rmq ( j + 1, e ) ];
					
					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				
				
				
				if( e == n )
					lcpde = 0;
				else lcpde = LCP[rmq ( d + 1, e ) ];
				
				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;

					/* lcp(e,j) */
					INT lcpej;
					
					
					if( j == n )
						lcpej = 0;
					else lcpej = LCP[rmq ( e + 1, j ) ];
					
					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				
				
				
				if( f == n )
					lcpef = 0;
				else lcpef = LCP[rmq ( e + 1, f ) ];
				
				if ( lcpef >= m )	f = std::min (f+1,n);

				interval.first = d + 1;
				interval.second = f - 1;
				return interval;


			}
			else if ( ( l == N - SA[i] ) || ( ( SA[i] + l < N ) && ( l != m ) && ( a[SA[i]+l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}
		}
	}
	

	interval.first = d + 1;
	interval.second = f - 1;
	return interval;
}

/* Computes the length of lcs of two suffixes of two strings */
INT lcs ( string & x, INT M, string & y, INT l )
{
	return lcs_zlteam(x, M, y, l);
	// if ( M < 0 ) return 0;
	// INT yy = y.size();
	// if ( l >= yy ) return 0;

	// INT i = 0;
	// while ( ( M - i >= 0 ) && ( l + i < yy ) )
	// {
	// 	if ( x[M-i] != y[l+i] )	break;
	// 	i++;
	// }
	// return i;
}


/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<INT,INT> rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n )
{
	
	
	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;

	pair<INT,INT> interval;

	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		INT revSA = N - 1 - SA[i];
		//std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		
		if( f == n )
			lcpif = 0;
		else lcpif = LCP[rmq ( i + 1, f ) ];
		
		/* lcp(d,i) */
		INT lcpdi;
		//it = rmq.find(make_pair(d+1, i));
		
		if( i == n )
			lcpdi = 0;
		else lcpdi = LCP[rmq ( d + 1, i ) ];
		
	
		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			
			// avoid the function call if revSA-1<0 or l>=w.size() by changing lcs?
			l = l + lcs ( a, revSA - l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		    	{
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;

					/* lcp(j,e) */
					INT lcpje;
					
					
					if( e == n )
						lcpje = 0;
					else lcpje = LCP[rmq ( j + 1, e ) ];
					
					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				
				
				if( e == n )
					lcpde = 0;
				else lcpde = LCP[rmq ( d + 1, e ) ];
				
				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );
			
				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;

					/* lcp(e,j) */
					INT lcpej;
					
					
					if( j == n )
						lcpej = 0;
					else lcpej = LCP[rmq ( e + 1, j ) ];
					
					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
					
				if( f == n )
					lcpef = 0;
				else lcpef = LCP[rmq ( e + 1, f ) ];
				
				if ( lcpef >= m )	f = std::min (f+1,n);

				interval.first = d + 1;
				interval.second = f - 1;
				return interval;


			}
			else if ( ( l == N - SA[i] ) || ( ( revSA - l >= 0 ) && ( l != m ) && ( a[revSA - l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}

		}
	}
	

	interval.first = d + 1;
	interval.second = f - 1;
	
	return interval;
}
