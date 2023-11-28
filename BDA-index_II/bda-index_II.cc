

#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_II.h"
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>
// #include "zlteam_re-imp_2.hpp"

using namespace std;

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif

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

// #endif


using namespace sdsl;

double vm, vm0, rss, rss0;

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}



/* Sorting comparison */
bool sort_sa(const pair<INT,INT> &a,const pair<INT,INT> &b)
{
       return a.first<b.first;
}

/* Constructs the right compacted trie given the anchors and the SA of the whole string in O(n) time */
void right_compacted_trie ( unordered_set<INT> &anchors, INT n, INT * RSA, INT * RLCP, INT g, INT ram_use, char * sa_fname, char * lcp_fname )
{
	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);
	
	INT ii = 0;
	INT minLCP = n;
	INT prevSA = 0;
	INT currSA = 0;
	INT currLCP = 0;
	
	for( INT i = 0; i <= n; i++ ) // in lex order
	{
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		prevSA = currSA;
		currSA = SA->read();
		currLCP = LCP->read();
		//cout<<currSA<<" "<<currLCP<<endl;
		auto it = anchors.find( currSA );
		if( it != anchors.end() )
		{
			
			RSA[ii] = currSA;		// store this suffix
			
			if ( ii == 0 )	RLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( prevSA == RSA[ii-1] )     // if the immediately prior suffix was added
					RLCP[ii] = currLCP;	// then the LCP value is the correct one for RLCP
				else
					RLCP[ii] = std::min(minLCP, currLCP);	//otherwise, we should take the minimum in the range
					
				
			}
  			//cout<<"RSA[i]: "<< RSA[ii] <<" RLCP[i]: "<<RLCP[ii]<<"\n"; getchar();
			minLCP = n; // set this to something high to get the _FIRST_ next minimum value, because a new range STARTS
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( currLCP < minLCP )
				minLCP = currLCP;
		}
	}
	
	delete(SA);
	delete(LCP);
}

/* Constructs the left compacted trie given the anchors and the SA of the whole string in O(n) time */
void left_compacted_trie ( unordered_set<INT> &anchors, INT n, INT * LSA, INT * LLCP, INT g, INT ram_use, char * sa_fname, char * lcp_fname )
{

	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);
	
		
	INT ii = 0;
	INT minLCP = n;
	INT prevSA = 0;
	INT currSA = 0;
	INT currLCP = 0;

	for( INT i = 0; i <= n; i++ ) // in lex order
	{
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		prevSA = currSA;
		currSA = SA->read();
		currLCP = LCP->read();
		auto it = anchors.find( ( n - 1 ) - currSA );
		//cout<<currSA<<" "<<currLCP<<endl;
		if( it != anchors.end() )
		{
			
			LSA[ii] = currSA;		// store this suffix

			if ( ii == 0 )	LLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( prevSA == LSA[ii-1] ) // if the immediately prior suffix was added
					LLCP[ii] = currLCP;	//then the LCP value is the correct one for RLCP
				else
					LLCP[ii] = std::min(minLCP, currLCP);	//otherwise, we should take the minimum in the range
				
			}

			//cout<<"LSA[i]: "<< RSA[ii] <<" LLCP[i]: "<< RLCP[ii]<<"\n"; getchar();
			minLCP = n; //set this to something high to get the FIRST next minimum value
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( currLCP < minLCP )
				minLCP = currLCP;
		}
	}
	
	delete(SA);
	delete(LCP);
}

int main(int argc, char **argv)
{
	unordered_set<char> alphabet;

	if( argc < 7 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./index <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>\n";
 		exit(-1);
 	}

 	ifstream is;
 	is.open (argv[1], ios::in | ios::binary);

 	std::string str2(argv[2]);

 	ifstream is2;
 	is2.open (argv[3], ios::in | ios::binary);

 	INT ell;
 	std::stringstream(str2)>>ell;
 	
 	char * output_filename; 
 	output_filename = (char *) malloc(strlen(argv[4])+1);    
    	strcpy(output_filename,argv[4]);

 	std::string str5(argv[5]);
 	INT ram_use;
 	std::stringstream(str5)>>ram_use;
 	
 	std::string str3(argv[6]);
 	
 	INT block;
 	std::stringstream(str3)>>block;
 	
 	char * index_filename; 
 	index_filename = (char *) malloc(strlen(argv[7])+1);    
    	strcpy(index_filename,argv[7]);
    	
    	ifstream is3;
 	is3.open (argv[7], ios::in | ios::binary);

 	ifstream in_file(argv[1], ios::binary);
   	in_file.seekg(0, ios::end);
   	INT text_file_size = in_file.tellg();
   	
   	
   	
  	char c = 0;
  	INT text_size = 0;
	for (INT i = 0; i < text_file_size; i++)
	{	
		is.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c == '\n')
			continue;
			
		else
		{
			alphabet.insert( (unsigned char) c );
			text_size++;
		}
		
	}
	is.close();
	
	INT k  = ceil(4*log2(ell)/log2(alphabet.size()));//
	if( ell - k - 1 < 0 )
		k = 2;
	
	if( text_size < block )
	{
	
		fprintf( stderr, " Error: Block size cannot be larger than sequence length!\n");
		return ( 1 );
	}
	
	if( block < ell )
	{
	
		fprintf( stderr, " Error: Block size cannot be smaller than window size!\n");
		return ( 1 );
	}
	
	if( text_size < ell )
	{
	
		fprintf( stderr, " Error: Window size (ell) cannot be larger than sequence length!\n");
		return ( 1 );
	}

    	std::chrono::steady_clock::time_point  start_bd = std::chrono::steady_clock::now();
 	unordered_set<INT> text_anchors;
	
	char * index_name =  ( char * ) malloc ( ( strlen( argv[7] ) + 1 ) * sizeof ( char ) );
	strcpy(index_name, argv[7]);
	const char * bd_anchors_suffix = ".bd";
	char * bd = strcat(index_name, bd_anchors_suffix);
 	
 
 	ifstream is_bd_anchors;
 	is_bd_anchors.open (index_name, ios::in | ios::binary);
 	
	ifstream in_bd_anchors(index_name, ios::binary);
 	in_bd_anchors.seekg (0, in_bd_anchors.end);
   	INT file_size = in_bd_anchors.tellg();
   	string bd_anchor = "";
   	INT bd_anchor_int = 0;
   	
   	INT * SA;
	INT * LCP;
	INT * invSA;
	INT * rank;
	
	rank = ( INT * ) malloc( ( block  ) *  sizeof( INT ) );
	SA = ( INT * ) malloc( ( block ) * sizeof( INT ) );

	
	if( ( SA == NULL) )
	{
	        fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
	        return ( 0 );
	}
	
	/* Compute the inverse SA array for block */
	invSA = ( INT * ) calloc( block , sizeof( INT ) );
	if( ( invSA == NULL) )
	{
	        fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
	        return ( 0 );
	}
	
	LCP = ( INT * ) calloc  ( block, sizeof( INT ) );
	if( ( LCP == NULL) )
	{
	        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
	        return ( 0 );
	}

   	if( file_size > 0 )
	{
		// Read in from .bd file
	    	c = 0;
		for (INT i = 0; i < file_size; i++)
		{	
			is_bd_anchors.read(reinterpret_cast<char*>(&c), 1);
		
			if( (unsigned char) c == '\n' )
			{
				bd_anchor_int = stoi( bd_anchor );
				text_anchors.insert( bd_anchor_int );
				bd_anchor = "";
			}
			else bd_anchor += (unsigned char) c;
			
		}
		is_bd_anchors.close();
	}
	else
	{
  		
		unsigned char * text_block = ( unsigned char * ) malloc (  ( block + 1 ) * sizeof ( unsigned char ) );
		unsigned char * suffix_block = ( unsigned char * ) malloc (  ( ell  ) * sizeof ( unsigned char ) );
		
		ifstream is_block;
	 	is_block.open (argv[1], ios::in | ios::binary);
	  	  
	 	c = 0;
	 	INT count = 0;
	 	INT pos = 0;
		for (INT i = 1; i < text_file_size; i++)
		{	
			is_block.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c != '\n' && (unsigned char) c != ' ' )
			{
				text_block[count] = (unsigned char) c ;
				count++;
				
				if( count == block || i == text_file_size - 1 )
				{
					text_block[count] = '\0';
					
					bd_anchors( text_block, pos, ell, k, text_anchors, SA, LCP, invSA, rank );
					
					memcpy( &suffix_block[0], &text_block[ block - ell + 1], ell -1 );
					memcpy( &text_block[0], &suffix_block[0], ell -1 );
					
					pos = pos + ( block - ell + 1 );
					count = ell - 1;
				}
			}
			
		}

			
		is_block.close();
		
		free( text_block );	
		free( suffix_block );
		
		ofstream bd_output;
		bd_output.open(index_name);
		
		for (auto &anchor : text_anchors)	
			bd_output<<anchor<<endl;
			
		bd_output.close();
	}
	
	free( rank );
	free( SA );
	free( invSA );
	free( LCP );
	
	INT g = text_anchors.size();
	INT n = text_size;
		
	std::chrono::steady_clock::time_point  end_bd = std::chrono::steady_clock::now();
	std::cout <<"bd construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_bd - start_bd).count() << "[ms]" << std::endl;
	cout<<"The text is of length "<< n << ", its alphabet size is "<< alphabet.size()<<", and it has "<<g<<" bd-anchors of order "<<ell<<endl;
	cout<<"The density is "<<(double) g / text_size<<endl;
	
	string text_string = "";
	ifstream is_full;
 	is_full.open (argv[1], ios::in | ios::binary);
 	
	c = 0;
	for (INT i = 0; i < text_size; i++)
	{	
		is_full.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c == '\n' )
			continue;
			
		else text_string.push_back( (unsigned char) c );
	}
	is_full.close();
	
	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	/* Constructing right and left compacted tries */
	INT * RSA;
	INT * RLCP;

	RSA = ( INT * ) malloc( ( g+1 ) * sizeof( INT ) );
	if( ( RSA == NULL) )
	{
	 	fprintf(stderr, " Error: Cannot allocate memory for RSA.\n" );
		return ( 0 );
	}

	RLCP = ( INT * ) malloc( ( g+1 ) * sizeof( INT ) );
	if( ( RLCP == NULL) )
	{
	 	fprintf(stderr, " Error: Cannot allocate memory for RLCP.\n" );
		return ( 0 );
	}
	
	const char * sa_suffix = ".RSA";
	char * index_name_0 =  ( char * ) malloc (  strlen( argv[7] ) * sizeof ( char ) );
	strcpy(index_name_0, argv[7]);
	strcat(index_name_0, sa_suffix);
	
	ifstream is_RSA;
 	is_RSA.open (index_name_0, ios::in | ios::binary);
 	
	ifstream in_RSA(index_name_0, ios::binary);

   	char sa_fname[strlen((const char *) index_filename)+20] ;
	sprintf(sa_fname, "%s_SA.sa5", index_filename);
	
	ifstream in_SA(sa_fname, ios::binary);
	in_RSA.seekg (0, in_RSA.end);
	INT file_size_sa = in_RSA.tellg();
	
	if( !(in_SA)  )
	{

	  	char commandesa[ strlen(sa_fname) + 1000 ];
	  	char * fullpathstart = dirname(realpath(argv[0], NULL));
	  	char command1[ strlen(sa_fname) + 1000 ];
	  	strcpy(command1, fullpathstart);
	  	strcat(command1, "/psascan/construct_sa %s -m %ldMi -o %s");
	  	sprintf(commandesa, command1, argv[1], ram_use, sa_fname);
	  	int outsa=system(commandesa);
	  	
	}
	
	if( file_size_sa > 0 )
	{
	   	string sa = "";
	   	INT sa_int = 0;
		c = 0;
		INT p = 0;
		
		for (INT i = 0; i < file_size; i++)
		{	
			is_RSA.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c == '\n' )
			{
				sa_int = stoi( sa);
				RSA[p] = sa_int;
				sa = "";
				p++;
			}
			else sa += (unsigned char) c;
			
		}
		is_RSA.close();
	
	
	}
	
	const char * lcp_suffix = ".RLCP";
	index_name =  ( char * ) malloc (  strlen( argv[7] ) * sizeof ( char ) );
	strcpy(index_name, argv[7]);
	strcat(index_name, lcp_suffix);
	
	ifstream is_RLCP;
 	is_RLCP.open (index_name, ios::in | ios::binary);
 	
	ifstream in_RLCP(index_name, ios::binary);
	in_RLCP.seekg (0, in_RLCP.end);
	file_size = in_RLCP.tellg();
	
   	
   	char lcp_fname[strlen((const char*) index_filename)+20] ;
	sprintf(lcp_fname, "%s_LCP.lcp5", index_filename);
	
	ifstream in_LCP(lcp_fname, ios::binary);
	

	if ( !(in_LCP) || file_size <= 0 || file_size_sa <= 0 )
	{
		if( !(in_LCP ) )
		{
			char commande[strlen(sa_fname) + strlen(lcp_fname) + 1000];
			char * fullpathstart = dirname(realpath(argv[0], NULL));
			char command2[strlen(sa_fname) + strlen(lcp_fname) + 1000];
			strcpy(command2, fullpathstart);
			strcat(command2, "/sparsePhi/src/construct_lcp_parallel -m %ldG -o %s -s %s %s");
			sprintf(commande, command2, ram_use, lcp_fname, sa_fname, argv[1]);
			int out=system(commande);
		}
		
	  	// right_compacted_trie ( text_anchors, n, RSA, RLCP, g, ram_use, sa_fname, lcp_fname );
	  	build_RSA_RLCP_zlteam(text_anchors, n, RSA, RLCP, ram_use, sa_fname, lcp_fname);


	  	ofstream rsa_output;
		rsa_output.open(index_name_0);
		
		for(INT i = 0; i<g; i++)	
			rsa_output<<RSA[i]<<endl;
			
		rsa_output.close();
		
		ofstream rlcp_output;
		rlcp_output.open(index_name);
		
		for(INT i = 0; i<g; i++)	
			rlcp_output<<RLCP[i]<<endl;
			
		rlcp_output.close();
		
		cout<<"Right Compacted trie constructed "<<endl;	
	}
	
	if( file_size > 0 )
	{
	   	string lcp = "";
	   	INT lcp_int = 0;
		c = 0;
		INT p = 0;
		
		for (INT i = 0; i < file_size; i++)
		{	
			is_RLCP.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c == '\n' )
			{
				lcp_int = stoi(lcp);
				RLCP[p] = lcp_int;
				lcp = "";
				p++;
			}
			else lcp += (unsigned char) c;
			
		}
		is_RLCP.close();
	
	
	}
	
	INT * LSA;
  	INT * LLCP;

  	LSA = ( INT * ) malloc( ( g+1 ) * sizeof( INT ) );
  	if( ( LSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LSA.\n" );
        	return ( 0 );
  	}
  	LLCP = ( INT * ) malloc( ( g+1 ) * sizeof( INT ) );
  	if( ( LLCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LLCP.\n" );
        	return ( 0 );
  	}
	
	char * output_reverse;
	const char * reversed_text = "_reverse";
	
	output_reverse = (char *) malloc(strlen(argv[1])+9);
	strcpy( output_reverse, argv[1]);
	strcat(output_reverse, reversed_text);
		
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	std::ofstream output_r;
  	output_r.open (output_reverse);
  	reverse(text_string.begin(), text_string.end());
  	output_r << text_string;
    	output_r.close();
 
	const char * sa_reverse_suffix = ".LSA";
	index_name_0 =  ( char * ) malloc (  strlen( argv[7] ) * sizeof ( char ) );
	strcpy(index_name_0, argv[7]);
	strcat(index_name_0, sa_reverse_suffix);
  	
  	
  	ifstream is_LSA;
 	is_LSA.open (index_name_0, ios::in | ios::binary);
 	
  	ifstream in_LSA(index_name_0, ios::binary);
  	
 
  	in_LSA.seekg (0, in_LSA.end);
  	file_size_sa = in_LSA.tellg();
  	
   	char sa_fname_reverse[strlen((const char *) index_filename)+20] ;
	sprintf(sa_fname_reverse, "%s%s_SA.sa5", index_filename, reversed_text);
	
	ifstream in_SA_reverse(sa_fname_reverse, ios::binary);
	
	
	if ( !(in_SA_reverse)  )
	{
	  	char commandesa_reverse[ strlen(sa_fname_reverse) + 1000 ];
	  	char * fullpathstart_reverse = dirname(realpath(argv[0], NULL));
	  	char command1_reverse[ strlen(sa_fname_reverse) + 1000 ];
	  	strcpy(command1_reverse, fullpathstart_reverse);
	  	strcat(command1_reverse, "/psascan/construct_sa %s -m %ldMi -o %s");
	  	sprintf(commandesa_reverse, command1_reverse, output_reverse, ram_use, sa_fname_reverse);
	  	int outsa_reverse=system(commandesa_reverse);
	}
	
	if( file_size_sa  > 0 )
	{
	   	string sa = "";
	   	INT sa_int = 0;
		c = 0;
		INT p = 0;
		
		for (INT i = 0; i < file_size; i++)
		{	
			is_LSA.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c == '\n' )
			{
				sa_int = stoi( sa);
				LSA[p] = sa_int;
				sa = "";
				p++;
			}
			else sa += (unsigned char) c;
			
		}
		is_LSA.close();	
	}
	

	
	const char * lcp_reverse_suffix = ".LLCP";
	index_name =  ( char * ) malloc (  strlen( argv[7] ) * sizeof ( char ) );
	strcpy(index_name, argv[7]);
	strcat(index_name, lcp_reverse_suffix);
	
	ifstream is_LLCP;
 	is_LLCP.open (index_name, ios::in | ios::binary);
	
	ifstream in_LLCP(index_name, ios::binary);
	in_LLCP.seekg (0, in_LLCP.end);
	file_size = in_LLCP.tellg();
   	
 
	char lcp_fname_reverse[strlen((const char*) index_filename)+20] ;
        sprintf(lcp_fname_reverse, "%s%s_LCP.lcp5", index_filename, reversed_text);
        
        ifstream in_LCP_reverse(lcp_fname_reverse, ios::binary);
	
	
        if( !(in_LCP_reverse) || file_size <= 0 || file_size_sa <= 0 )
	{
	
		if( !(in_LCP_reverse ) )
		{
			char commande_reverse[strlen(sa_fname_reverse) + strlen(lcp_fname_reverse) + 1000];
			char * fullpathstart_reverse = dirname(realpath(argv[0], NULL));
			char command2_reverse[strlen(sa_fname_reverse) + strlen(lcp_fname_reverse) + 1000];
			strcpy(command2_reverse, fullpathstart_reverse);
			strcat(command2_reverse, "/sparsePhi/src/construct_lcp_parallel -m %ldG -o %s -s %s %s");
			sprintf(commande_reverse, command2_reverse, ram_use, lcp_fname_reverse, sa_fname_reverse, output_reverse);
			int out_reverse=system(commande_reverse);
		}
		
		// left_compacted_trie ( text_anchors, n, LSA, LLCP, g, ram_use, sa_fname_reverse, lcp_fname_reverse );
		build_LSA_LLCP_zlteam(text_anchors, n, LSA, LLCP, ram_use, sa_fname_reverse, lcp_fname_reverse);
  		
  		
  		ofstream lsa_output;
		lsa_output.open(index_name_0);
		
		for(INT i = 0; i<g; i++)	
			lsa_output<<LSA[i]<<endl;
			
		lsa_output.close();
		
		ofstream llcp_output;
		llcp_output.open(index_name);
		
		for(INT i = 0; i<g; i++)	
			llcp_output<<LLCP[i]<<endl;
			
		llcp_output.close();
	
		cout<<"Left Compacted trie constructed"<<endl;
	}
	
	
	
	if( file_size > 0 )
	{
	   	string lcp = "";
	   	INT lcp_int = 0;
		c = 0;
		INT p = 0;
		
		for (INT i = 0; i < file_size; i++)
		{	
			is_LLCP.read(reinterpret_cast<char*>(&c), 1);
			
			if( (unsigned char) c == '\n' )
			{
				lcp_int = stoi(lcp);
				LLCP[p] = lcp_int;
				lcp = "";
				p++;

			}
			else lcp += (unsigned char) c;
			
		}
		is_LLCP.close();
	}
	

  	/* After constructing the tries these DSs over the whole string are not needed anymore, our data structure must be of size O(g) */
  	text_anchors.clear();

	
  	/* The following RMQ data structures are used for spelling pattern over the LSA and RSA */
  
	const char * rmq_left_suffix = ".lrmq";
	index_name =  ( char * ) malloc (  strlen( argv[7] ) * sizeof ( char ) );
	strcpy(index_name, argv[7]);
	strcat(index_name, rmq_left_suffix);
  		
	ifstream in_rmq_left(index_name, ios::binary);
   
  	rmq_succinct_sct<> lrmq;
  	
  	
  	if( in_rmq_left )
  	{
		cout<<"load rmq "<<endl;  	
  		load_from_file(lrmq, index_name); 
	}
  	else
  	{
	  	int_vector<> llcp( g , 0 ); // create a vector of length n and initialize it with 0s

		for ( INT i = 0; i < g; i ++ )
		{
			llcp[i] = LLCP[i];
		}

		util::assign(lrmq, rmq_succinct_sct<>(&llcp));
		util::clear(llcp);
		
		cout<<"Left RMQ DS constructed "<<endl;	
			
		store_to_file(lrmq, index_name);
	}
	
	
	const char * rmq_right_suffix = ".rrmq";
	index_name =  ( char * ) malloc (  strlen( argv[7] ) * sizeof ( char ) );
	strcpy(index_name, argv[7]);
	strcat(index_name, rmq_right_suffix);
	
	ifstream in_rmq_right(index_name, ios::binary);
  	  
  	int_vector<> rlcp( g , 0 ); // create a vector of length n and initialize it with 0s
  	rmq_succinct_sct<> rrmq;
  	
  	if( in_rmq_right )
  	{

  		load_from_file(rrmq, index_name); 
  	
  	}
  	else
  	{
		int_vector<> rlcp( g , 0 ); // create a vector of length n and initialize it with 0s

		for ( INT i = 0; i < g; i ++ )
		{
			rlcp[i] = RLCP[i];
		}
		
		util::assign(rrmq, rmq_succinct_sct<>(&rlcp));
		util::clear(rlcp);
		 

	  	cout<<"Right RMQ DS constructed "<<endl;
	  	
	  	store_to_file(rrmq, index_name);
	}
	 
    	cout<<"The whole index is constructed"<<endl;
	
	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"Index took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index- start_index).count() << "[ms]" << std::endl;
	  
	std::chrono::steady_clock::time_point  begin_pt = std::chrono::steady_clock::now();
    	reverse(text_string.begin(), text_string.end()); 				//I re-reverse to take the original string
  
  	vector<vector<unsigned char> > all_patterns;
    	vector<unsigned char> pattern;
    	c = 0;
    	while (is2.read(reinterpret_cast<char*>(&c), 1))
    	{
        	if(c == '\n')
        	{
  			if(pattern.empty())	break;
  			all_patterns.push_back(pattern);
  			pattern.clear();
        	}
        	else	pattern.push_back((unsigned char)c);
    	}
    	is2.close();
    	pattern.clear();

	vector<string> new_all_pat;
	for(auto &it_pat : all_patterns)	new_all_pat.push_back(string(it_pat.begin(), it_pat.end()));
	all_patterns.clear();
	
	INT hits=0;
    
	INT *f = new INT[ell<<1];
  	std::chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();
  	ofstream pattern_output;
	pattern_output.open(output_filename);
	for(auto &pattern : new_all_pat)
   	{
  		if ( pattern.size() < ell )
  		{
  			pattern_output<<"Pattern skipped: its length is less than ell!\n";
  			continue;
  		}
		//num_of_patterns++;
		
		string first_window = pattern.substr(0, ell).c_str();
  		INT j = red_minlexrot( first_window, f, ell, k );
  		
		if ( pattern.size() - j >= j ) //if the right part is bigger than the left part, then search the right part to get a smaller interval on RSA (on average)
		{ 
  			string right_pattern = pattern.substr(j, pattern.size()-j);
			pair<INT,INT> right_interval = pattern_matching ( right_pattern, text_string, RSA, RLCP, rrmq, g );
  														

			if(right_interval.first > right_interval.second)	continue;
		
			for(INT i = right_interval.first; i <= right_interval.second; i++ ) //this can be a large interval and only one occurrence is valid.
			{
				INT index = RSA[i];
				INT jj = j;		//this is the index of the anchor in the pattern
				index--; 	jj--;	//jump the index of the anchor and start looking on the left
				while ( ( jj >= 0 ) && ( index >= 0 ) && ( text_string[index] == pattern[jj] ) )
				{
					index--; jj--;
				}
				if ( jj < 0 ) //we have matched the pattern completely
				{
					pattern_output<< pattern <<" found at position "<< index + 1 << " of the text"<<endl;
				}					
			}
		}
		else //otherwise, search the left part to get a smaller interval on LSA (on average)
		{ 
			string left_pattern = pattern.substr(0, j+1);
			reverse(left_pattern.begin(), left_pattern.end());
			pair<INT,INT> left_interval = rev_pattern_matching ( left_pattern, text_string, LSA, LLCP, lrmq, g );
  														

			if(left_interval.first > left_interval.second)	continue;
		
			for(INT i = left_interval.first; i <= left_interval.second; i++ ) //this can be a large interval and only one occurrence is valid.
			{
				INT index = n-1-LSA[i];
				INT jj = j;		//this is the index of the anchor in the pattern
				index++; 	jj++;	//jump the index of the anchor and start looking on the right
				while ( ( jj < pattern.size() ) && ( index < n ) && ( text_string[index] == pattern[jj] ) )
				{
					index++; jj++;
				}
				if ( jj == pattern.size() ) //we have matched the pattern completely
				{ 
					if ( index == n - 1 )	pattern_output<< pattern <<" found at position "<< index - pattern.size() + 1 << " of the text"<<endl;					
					else			pattern_output<< pattern <<" found at position "<<  index - pattern.size() << " of the text"<<endl;
				}
			}
		}
				
   	}
 	
	std::chrono::steady_clock::time_point  end_pt = std::chrono::steady_clock::now();
	std::cout <<"Pattern matching took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pt - begin_pt).count() << "[ms]" << std::endl;
  	free ( RSA );
  	free ( RLCP );
  	free ( LSA );
  	free ( LLCP );

	return 0;
}

