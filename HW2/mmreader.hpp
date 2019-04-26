/*
Portions of this file include code provided by The National Institute of
Standards and Technology (NIST).  The code includes
macro definitions from mmio.h and is subject to the following disclaimer.

Software Disclaimer

NIST-developed software is provided by NIST as a public service. You may use,
copy and distribute copies of the software in any medium, provided that you
keep intact this entire notice. You may improve, modify and create derivative
works of the software or any portion of the software, and you may copy and
distribute such modifications or works. Modified works should carry a notice
stating that you changed the software and should note the date and nature of
any such change. Please explicitly acknowledge the National Institute of
Standards and Technology as the source of the software.

NIST-developed software is expressly provided "AS IS" NIST MAKES NO WARRANTY
OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY OPERATION OF LAW,
INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA ACCURACY. NIST
NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE
UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES
NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR
THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY,
RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using and
distributing the software and you assume all risks associated with its use,
including but not limited to the risks and costs of program errors, compliance
with applicable laws, damage to or loss of data, programs or equipment, and
the unavailability or interruption of operation. This software is not intended
to be used in any situation where a failure could cause risk of injury or
damage to property. The software developed by NIST employees is not subject
to copyright protection within the United States.
*/
/* Modified date: 2019-03-19. Modification for a custom sparse matrix structure. */

#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <typeinfo>
#include <algorithm>
#include <cinttypes>

#define mm_is_matrix(typecode)	((typecode)[0]=='M')
#define mm_is_sparse(typecode)	((typecode)[1]=='C')
#define mm_is_coordinate(typecode)((typecode)[1]=='C')
#define mm_is_dense(typecode)	((typecode)[1]=='A')
#define mm_is_array(typecode)	((typecode)[1]=='A')

#define mm_is_complex(typecode)	((typecode)[2]=='C')
#define mm_is_real(typecode)		((typecode)[2]=='R')
#define mm_is_pattern(typecode)	((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode)	((typecode)[3]=='G')
#define mm_is_skew(typecode)	((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

/********************* MM_typecode modify fucntions ***************************/

#define mm_set_matrix(typecode)	((typecode)[0]='M')
#define mm_set_coordinate(typecode)	((typecode)[1]='C')
#define mm_set_array(typecode)	((typecode)[1]='A')
#define mm_set_dense(typecode)	mm_set_array(typecode)
#define mm_set_sparse(typecode)	mm_set_coordinate(typecode)

#define mm_set_complex(typecode)((typecode)[2]='C')
#define mm_set_real(typecode)	((typecode)[2]='R')
#define mm_set_pattern(typecode)((typecode)[2]='P')
#define mm_set_integer(typecode)((typecode)[2]='I')

#define mm_set_symmetric(typecode)((typecode)[3]='S')
#define mm_set_general(typecode)((typecode)[3]='G')
#define mm_set_skew(typecode)	((typecode)[3]='K')
#define mm_set_hermitian(typecode)((typecode)[3]='H')

#define mm_clear_typecode(typecode) ((typecode)[0]=(typecode)[1]= \
                                    (typecode)[2]=' ',(typecode)[3]='G')
#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)

/********************* Matrix Market error codes ***************************/

#define MM_COULD_NOT_READ_FILE  11
#define MM_PREMATURE_EOF        12
#define MM_NOT_MTX              13
#define MM_NO_HEADER            14
#define MM_UNSUPPORTED_TYPE     15
#define MM_LINE_TOO_LONG        16
#define MM_COULD_NOT_WRITE_FILE 17

#define MM_MTX_STR      "matrix"
#define MM_ARRAY_STR    "array"
#define MM_DENSE_STR    "array"
#define MM_COORDINATE_STR "coordinate"
#define MM_SPARSE_STR   "coordinate"
#define MM_COMPLEX_STR  "complex"
#define MM_REAL_STR     "real"
#define MM_INT_STR      "integer"
#define MM_GENERAL_STR  "general"
#define MM_SYMM_STR     "symmetric"
#define MM_HERM_STR     "hermitian"
#define MM_SKEW_STR     "skew-symmetric"
#define MM_PATTERN_STR  "pattern"

#define MM_MAX_LINE_LENGTH 1025
#define MM_MAX_TOKEN_LENGTH 64

#define MatrixMarketBanner "%%MatrixMarket"

#define MAX_RAND_VAL 5.0

struct Coordinate
{
    int32_t x;
    int32_t y;
    float val;
};

bool CoordinateCompare( const Coordinate &c1, const Coordinate &c2 );

/** @struct Sparse Matrix */
struct sparse_mtx
{
    uint32_t nrow; /**< Number of rows */
    uint32_t ncol; /**< Number of columns */
    uint32_t nnze; /**< Number of non-zero elements. */
    int32_t *row; /**< Row index of this matrix. */
    int32_t *col; /**< Column index of this matrix. */
    float *val; /**< Value of this matrix. */
};

/** @struct Dense Matrix */
struct dense_mtx
{
    uint32_t nrow; /**< Number of rows */
    uint32_t ncol; /**< Number of columns */
    float *val; /**< Value of this matrix. */
};

class MatrixMarketReader
{
    char Typecode[ 4 ];
    int32_t nNZ;
    int32_t nRows;
    int32_t nCols;
    int isSymmetric;
    int isDoubleMem;
    Coordinate *unsym_coords;

public:
    MatrixMarketReader( ): nNZ( 0 ), nRows( 0 ), nCols( 0 ), isSymmetric( 0 ), isDoubleMem( 0 )
    {
        int i;
        for(i = 0; i < 4; i++)
            Typecode[ i ] = '\0';

        unsym_coords = NULL;
    }

    int MMReadHeader( FILE* infile );
    int MMReadHeader( const std::string& filename );
    int MMReadFormat( const std::string& _filename);
    int MMReadBanner( FILE* infile );
    int MMReadMtxCrdSize( FILE* infile );
    void MMGenerateCOOFromFile( FILE* infile);

    int32_t GetNumRows( )
    {
        return nRows;
    }

    int32_t GetNumCols( )
    {
        return nCols;
    }

    int32_t GetNumNonZeroes( )
    {
        return nNZ;
    }

    int GetSymmetric( )
    {
        return isSymmetric;
    }

    Coordinate *GetUnsymCoordinates( )
    {
        return unsym_coords;
    }

    ~MatrixMarketReader( )
    {
        delete[ ] unsym_coords;
    }
};
