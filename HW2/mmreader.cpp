#include "mmreader.hpp"

using namespace std;

bool CoordinateCompare( const Coordinate &c1, const Coordinate &c2 )
{
    if( c1.x != c2.x )
        return ( c1.x < c2.x );
    else
        return ( c1.y < c2.y );
}

// Class definition

int MatrixMarketReader::MMReadHeader( FILE* mm_file )
{
    int status = MMReadBanner( mm_file );
    if( status != 0 )
    {
        printf( "Error Reading Banner in Matrix-Market File !\n" );
        return 1;
    }

    if( !mm_is_coordinate( Typecode ) )
    {
        printf( "Handling only coordinate format\n" ); return( 1 );
    }

    if( mm_is_complex( Typecode ) )
    {
        printf( "Error: cannot handle complex format\n" );
        return ( 1 );
    }

    if( mm_is_symmetric( Typecode ) )
        isSymmetric = 1;

    status = MMReadMtxCrdSize( mm_file );
    if( status != 0 )
    {
        printf( "Error reading Matrix Market crd_size %d\n", status );
        return( 1 );
    }

    return 0;
}

int MatrixMarketReader::MMReadHeader( const std::string &filename )
{
    FILE *mm_file = ::fopen( filename.c_str( ), "r" );
    if( mm_file == NULL )
    {
        printf( "Cannot Open Matrix-Market File !\n" );
        return 1;
    }

    if ( MMReadHeader( mm_file ) )
    {
        printf ("Matrix not supported !\n");
        return 2;
    }

    // If symmetric MM stored file, double the reported size
    if( mm_is_symmetric( Typecode ) )
        nNZ <<= 1;

    ::fclose( mm_file );

    std::clog << "Matrix: " << filename << " [nRow: " << GetNumRows( ) << "] [nCol: " << GetNumCols( ) << "] [nNZ: " << GetNumNonZeroes( ) << "]" << std::endl;

    return 0;
}

int MatrixMarketReader::MMReadFormat( const std::string &filename)
{
    FILE *mm_file = ::fopen( filename.c_str( ), "r" );
    if( mm_file == NULL )
    {
        printf( "Cannot Open Matrix-Market File !\n" );
        return 1;
    }

    if ( MMReadHeader( mm_file ) )
    {
        printf ("Matrix not supported !\n");
        return 2;
    }

    if( mm_is_symmetric( Typecode ) )
        unsym_coords = new Coordinate[ 2 * nNZ ];
    else
        unsym_coords = new Coordinate[ nNZ ];

    MMGenerateCOOFromFile( mm_file);
    ::fclose( mm_file );

    return 0;
}

void FillCoordData( char Typecode[ ],
                    Coordinate *unsym_coords,
                    int32_t &unsym_actual_nnz,
                    int32_t ir,
                    int32_t ic,
                    float val )
{
    if( mm_is_symmetric( Typecode ) )
    {
        unsym_coords[ unsym_actual_nnz ].x = ir - 1;
        unsym_coords[ unsym_actual_nnz ].y = ic - 1;
        unsym_coords[ unsym_actual_nnz++ ].val = val;

        if( unsym_coords[ unsym_actual_nnz - 1 ].x != unsym_coords[ unsym_actual_nnz - 1 ].y )
        {
            unsym_coords[ unsym_actual_nnz ].x = unsym_coords[ unsym_actual_nnz - 1 ].y;
            unsym_coords[ unsym_actual_nnz ].y = unsym_coords[ unsym_actual_nnz - 1 ].x;
            unsym_coords[ unsym_actual_nnz ].val = unsym_coords[ unsym_actual_nnz - 1 ].val;
            unsym_actual_nnz++;
        }
    }
    else
    {
        unsym_coords[ unsym_actual_nnz ].x = ir - 1;
        unsym_coords[ unsym_actual_nnz ].y = ic - 1;
        unsym_coords[ unsym_actual_nnz++ ].val = val;
    }
}

void MatrixMarketReader::MMGenerateCOOFromFile(FILE *infile)
{
    int32_t unsym_actual_nnz = 0;
    float val;
    int32_t ir, ic;

    //silence warnings from fscanf (-Wunused-result)
    int rv = 0;

    for ( int32_t i = 0; i < nNZ; i++)
    {
        if( mm_is_real( Typecode ) )
        {
            rv = fscanf(infile, "%" SCNuLEAST32, &ir);
            if(rv == EOF)
                printf("failed to read.\n");
            rv = fscanf(infile, "%" SCNuLEAST32, &ic);
            if(rv == EOF)
                printf("failed to read.\n");
            rv = fscanf(infile, "%f\n", (float*)(&val));
            if(rv == EOF)
                printf("failed to read.\n");
            if(val == 0)
                continue;
            else
                FillCoordData( Typecode, unsym_coords, unsym_actual_nnz, ir, ic, val );
        }
        else if( mm_is_integer( Typecode ) )
        {
            rv = fscanf(infile, "%" SCNuLEAST32, &ir);
            if(rv == EOF)
                printf("failed to read.\n");
            rv = fscanf(infile, "%" SCNuLEAST32, &ic);
            if(rv == EOF)
                printf("failed to read.\n");
            rv = fscanf(infile, "%f\n", (float*)( &val ) );
            if(rv == EOF)
                printf("failed to read.\n");
            if(val == 0)
                continue;
            else
                FillCoordData( Typecode, unsym_coords, unsym_actual_nnz, ir, ic, val );

        }
        else if( mm_is_pattern( Typecode ) )
        {
            rv = fscanf(infile, "%" SCNuLEAST32, &ir);
            if(rv == EOF)
                printf("failed to read.\n");
            rv = fscanf(infile, "%" SCNuLEAST32, &ic);
            if(rv == EOF)
                printf("failed to read.\n");
            val = 1.0f;

            if(val == 0)
                continue;
            else
                FillCoordData( Typecode, unsym_coords, unsym_actual_nnz, ir, ic, val );
        }
    }
    nNZ = unsym_actual_nnz;
}

int MatrixMarketReader::MMReadBanner( FILE *infile )
{
    char line[ MM_MAX_LINE_LENGTH ];
    char banner[ MM_MAX_TOKEN_LENGTH ];
    char mtx[ MM_MAX_TOKEN_LENGTH ];
    char crd[ MM_MAX_TOKEN_LENGTH ];
    char data_type[ MM_MAX_TOKEN_LENGTH ];
    char storage_scheme[ MM_MAX_TOKEN_LENGTH ];
    char *p;

    mm_clear_typecode( Typecode );

    if( fgets( line, MM_MAX_LINE_LENGTH, infile ) == NULL )
        return MM_PREMATURE_EOF;

    if( sscanf( line, "%s %s %s %s %s", banner, mtx, crd, data_type,
        storage_scheme ) != 5 )
        return MM_PREMATURE_EOF;

    for( p = mtx; *p != '\0'; *p = tolower( *p ), p++ );  /* convert to lower case */
    for( p = crd; *p != '\0'; *p = tolower( *p ), p++ );
    for( p = data_type; *p != '\0'; *p = tolower( *p ), p++ );
    for( p = storage_scheme; *p != '\0'; *p = tolower( *p ), p++ );

    /* check for banner */
    if( strncmp( banner, MatrixMarketBanner, strlen( MatrixMarketBanner ) ) != 0 )
        return MM_NO_HEADER;

    /* first field should be "mtx" */
    if( strcmp( mtx, MM_MTX_STR ) != 0 )
        return  MM_UNSUPPORTED_TYPE;
    mm_set_matrix( Typecode );

    /* second field describes whether this is a sparse matrix (in coordinate
            storgae) or a dense array */
    if( strcmp( crd, MM_SPARSE_STR ) == 0 )
        mm_set_sparse( Typecode );
    else if( strcmp( crd, MM_DENSE_STR ) == 0 )
        mm_set_dense( Typecode );
    else
        return MM_UNSUPPORTED_TYPE;


    /* third field */

    if( strcmp( data_type, MM_REAL_STR ) == 0 )
        mm_set_real( Typecode );
    else
        if( strcmp( data_type, MM_COMPLEX_STR ) == 0 )
            mm_set_complex( Typecode );
        else
            if( strcmp( data_type, MM_PATTERN_STR ) == 0 )
                mm_set_pattern( Typecode );
            else
                if( strcmp( data_type, MM_INT_STR ) == 0 )
                    mm_set_integer( Typecode );
                else
                    return MM_UNSUPPORTED_TYPE;


    /* fourth field */

    if( strcmp( storage_scheme, MM_GENERAL_STR ) == 0 )
        mm_set_general( Typecode );
    else
        if( strcmp( storage_scheme, MM_SYMM_STR ) == 0 )
            mm_set_symmetric( Typecode );
        else
            if( strcmp( storage_scheme, MM_HERM_STR ) == 0 )
                mm_set_hermitian( Typecode );
            else
                if( strcmp( storage_scheme, MM_SKEW_STR ) == 0 )
                    mm_set_skew( Typecode );
                else
                    return MM_UNSUPPORTED_TYPE;

    return 0;

}

int MatrixMarketReader::MMReadMtxCrdSize( FILE *infile )
{
    char line[ MM_MAX_LINE_LENGTH ];
    int num_items_read;

    /* now continue scanning until you reach the end-of-comments */
    do
    {
        if( fgets( line, MM_MAX_LINE_LENGTH, infile ) == NULL )
            return MM_PREMATURE_EOF;
    } while( line[ 0 ] == '%' );

    /* line[] is either blank or has M,N, nz */
    std::stringstream s(line);
    nRows = 0;
    nCols = 0;
    nNZ   = 0;
    s >> nRows >> nCols >> nNZ;
    if (nRows && nCols && nNZ)
        return 0;
    else
        do
        {
            num_items_read = 0;
            num_items_read += fscanf( infile, "%" SCNuLEAST32, &nRows );
            if (num_items_read == EOF) return MM_PREMATURE_EOF;
            num_items_read += fscanf(infile,  "%" SCNuLEAST32, &nCols);
            if (num_items_read == EOF) return MM_PREMATURE_EOF;
            num_items_read += fscanf(infile,  "%" SCNuLEAST32, &nNZ);
            if( num_items_read == EOF ) return MM_PREMATURE_EOF;
        } while( num_items_read != 3 );

    return 0;
}
