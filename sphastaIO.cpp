#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

#include "sphastaIO.h"

#include "rdtsc.h"
#define clockRate 2670000000.0


unsigned long long start, end;
double opentime_total = 0.0;
int LAST_FILE_ID;

void startTimer(unsigned long long* start) {
	*start =  rdtsc();
}

void endTimer(unsigned long long* end) {
	*end = rdtsc();
}

void computeTime(unsigned long long* start, unsigned long long* end) {
	double time = (double)((*end-*start)/clockRate);
	opentime_total += time;
}


#define VERSION_INFO_HEADER_SIZE 8192
#define DB_HEADER_SIZE 1024
#define TWO_MEGABYTE 2097152
#define ENDIAN_TEST_NUMBER 12180 // Troy's Zip Code!!
#define MAX_PHASTA_FILES 64
#define MAX_PHASTA_FILE_NAME_LENGTH 1024
#define MAX_FIELDS_NUMBER 48
#define MAX_FIELDS_NAME_LENGTH 128
#define DefaultMHSize (4*1024*1024)
int MasterHeaderSize = DefaultMHSize;
int diff_endian = 0;
long long counter = 0;

enum PhastaIO_Errors
{
	MAX_PHASTA_FILES_EXCEEDED = -1,
	UNABLE_TO_OPEN_FILE = -2,
	NOT_A_MPI_FILE = -3,
	GPID_EXCEEDED = -4,
	DATA_TYPE_ILLEGAL = -5,
};

int partID_counter;

typedef struct
{
	bool Wrong_Endian;                            /* default to false */

	char filename[MAX_PHASTA_FILE_NAME_LENGTH];   /* defafults to 1024 */
	int nppp;
	int nPPF;
	int nFiles;
	int nFields;
	unsigned long long my_offset;

	char * master_header;
	double * double_chunk;
	int * int_chunk;
	double * read_double_chunk;
	int * read_int_chunk;
	unsigned long long **my_offset_table;
	unsigned long long **my_read_table;
	int field_count;
	int part_count;
	int read_field_count;
	int read_part_count;
	int GPid;
	int start_id;
	unsigned long long next_start_address;

	int myrank;
	int numprocs;
	int local_myrank;
	int local_numprocs;
} phastaio_file_t;

//default: Paraview disabled

typedef struct
{
	int fileID;
	int nppf, nfields;
	int GPid;
	int read_field_count;
	char * masterHeader;
	unsigned long long **offset_table;
	unsigned long long my_offset;

}serial_file;

serial_file *SerialFile;
phastaio_file_t *PhastaIOActiveFiles[MAX_PHASTA_FILES];
int PhastaIONextActiveIndex = 0; /* indicates next index to allocate */

#define swap_char(A,B) { ucTmp = A; A = B ; B = ucTmp; }

std::map< int , char* > LastHeaderKey;
std::vector< FILE* > fileArray;
std::vector< int > byte_order;
std::vector< int > header_type;
int DataSize=0;
int LastHeaderNotFound = 0;
int Wrong_Endian = 0 ;
int Strict_Error = 0 ;
int binary_format = 0;

// the caller has the responsibility to delete the returned string
char* StringStripper( const char  istring[] )
{
	int length = strlen( istring );
	char* dest = new char [ length + 1 ];
	strcpy( dest, istring );
	dest[ length ] = '\0';

	if ( char* p = strpbrk( dest, " ") )
	{
		*p = '\0';
	}

	return dest;
}

int cscompare( const char teststring[],
		const char targetstring[] )
{

	char* s1 = const_cast<char*>(teststring);
	char* s2 = const_cast<char*>(targetstring);

	while( *s1 == ' ') { s1++; }
	while( *s2 == ' ') { s2++; }
	while( ( *s1 )
			&& ( *s2 )
			&& ( *s2 != '?')
			&& ( tolower( *s1 )==tolower( *s2 ) ) )
	{
		s1++;
		s2++;
		while( *s1 == ' ') { s1++; }
		while( *s2 == ' ') { s2++; }
	}
	if ( !( *s1 ) || ( *s1 == '?') )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void isBinary( const char iotype[] )
{
	char* fname = StringStripper( iotype );
	if ( cscompare( fname, "binary" ) )
	{
		binary_format = 1;
	}
	else
	{
		binary_format = 0;
	}
	delete [] fname;

}

size_t typeSize( const char typestring[] )
{
	char* ts1 = StringStripper( typestring );

	if ( cscompare( "integer", ts1 ) )
	{
		delete [] ts1;
		return sizeof(int);
	}
	else if ( cscompare( "double", ts1 ) )
	{
		delete [] ts1;
		return sizeof( double );
	}
	else if ( cscompare( "float", ts1 ) )
	{
		delete [] ts1;
		return sizeof( float );
	}
	else
	{
		delete [] ts1;
		fprintf(stderr,"unknown type : %s\n",ts1);
		return 0;
	}
}

int readHeader( FILE*       fileObject,
		const char  phrase[],
		int*        params,
		int         expect )
{
	char* text_header;
	char* token;
	char Line[1024];
	char junk;
	int FOUND = 0 ;
	int real_length;
	int skip_size, integer_value;
	int rewind_count=0;

	if( !fgets( Line, 1024, fileObject ) && feof( fileObject ) )
	{
		rewind( fileObject );
		clearerr( fileObject );
		rewind_count++;
		fgets( Line, 1024, fileObject );
	}

	while( !FOUND  && ( rewind_count < 2 ) )
	{
		if ( ( Line[0] != '\n' ) && ( real_length = strcspn( Line, "#" )) )
		{
			text_header = new char [ real_length + 1 ];
			strncpy( text_header, Line, real_length );
			text_header[ real_length ] =static_cast<char>(NULL);
			token = strtok ( text_header, ":" );
			if( cscompare( phrase , token ) )
			{
				FOUND = 1 ;
				token = strtok( NULL, " ,;<>" );
				skip_size = atoi( token );
				int i;
				for( i=0; i < expect && ( token = strtok( NULL," ,;<>") ); i++)
				{
					params[i] = atoi( token );
				}
				if ( i < expect )
				{
					fprintf(stderr,"Expected # of ints not found for: %s\n",phrase );
				}
			}
			else if ( cscompare(token,"byteorder magic number") )
			{
				if ( binary_format )
				{
					fread((void*)&integer_value,sizeof(int),1,fileObject);
					fread( &junk, sizeof(char), 1 , fileObject );
					if ( 362436 != integer_value )
					{
						Wrong_Endian = 1;
					}
				}
				else
				{
					fscanf(fileObject, "%d\n", &integer_value );
				}
			}
			else
			{
				/* some other header, so just skip over */
				token = strtok( NULL, " ,;<>" );
				skip_size = atoi( token );
				if ( binary_format)
				{
					fseek( fileObject, skip_size, SEEK_CUR );
				}
				else
				{
					for( int gama=0; gama < skip_size; gama++ )
					{
						fgets( Line, 1024, fileObject );
					}
				}
			}
			delete [] text_header;
		}

		if ( !FOUND )
		{
			if( !fgets( Line, 1024, fileObject ) && feof( fileObject ) )
			{
				rewind( fileObject );
				clearerr( fileObject );
				rewind_count++;
				fgets( Line, 1024, fileObject );
			}
		}
	}

	if ( !FOUND )
	{
		fprintf(stderr, "Error: Cound not find: %s\n", phrase);
		return 1;
	}
	return 0;
}

void SwapArrayByteOrder( void* array,
		int   nbytes,
		int   nItems )
{
	/* This swaps the byte order for the array of nItems each
		 of size nbytes , This will be called only locally  */
	int i,j;
	unsigned char ucTmp;
	unsigned char* ucDst = (unsigned char*)array;

	for(i=0; i < nItems; i++)
	{
		for(j=0; j < (nbytes/2); j++)
		{
			swap_char( ucDst[j] , ucDst[(nbytes - 1) - j] );
		}
		ucDst += nbytes;
	}
}

void queryphmpiio(const char filename[],int *nfields, int *nppf)
{

	FILE * fileHandle;
	char* fname = StringStripper( filename );

	fileHandle = fopen (fname,"rb");
	if (fileHandle == NULL ) {
		printf("\n File %s doesn't exist! Please check!\n",fname);
    exit(1);
  }
	else
	{
		SerialFile =(serial_file *)calloc( 1,  sizeof( serial_file) );

		SerialFile->masterHeader = (char *)malloc(MasterHeaderSize);
		fread(SerialFile->masterHeader,1,MasterHeaderSize,fileHandle);

		char read_out_tag[MAX_FIELDS_NAME_LENGTH];
		char * token;
		int magic_number;
		memcpy( read_out_tag,
				SerialFile->masterHeader,
				MAX_FIELDS_NAME_LENGTH-1 );

		if ( cscompare ("MPI_IO_Tag",read_out_tag) )
		{
			// Test endianess ...
			memcpy ( &magic_number,
					SerialFile->masterHeader+sizeof("MPI_IO_Tag :"),
					sizeof(int) );

			if ( magic_number != ENDIAN_TEST_NUMBER )
			{
				diff_endian = 1;
			}

			char version[MAX_FIELDS_NAME_LENGTH/4];
			int mhsize;

			memcpy(version,
					SerialFile->masterHeader + MAX_FIELDS_NAME_LENGTH/2,
					MAX_FIELDS_NAME_LENGTH/4 - 1); //TODO: why -1?

			if( cscompare ("version",version) )
			{
				// if there is "version" tag in the file, then it is newer format
				// read master header size from here, otherwise use default
				// TODO: if version is "1", we know mhsize is at 3/4 place...

				token = strtok(version, ":");
				token = strtok(NULL, " ,;<>" );
				int iversion = atoi(token);

				if( iversion == 1) {
					memcpy( &mhsize,
							SerialFile->masterHeader + MAX_FIELDS_NAME_LENGTH/4*3 + sizeof("mhsize : ")-1,
							sizeof(int));
					if ( diff_endian)
						SwapArrayByteOrder(&mhsize, sizeof(int), 1);
					free(SerialFile->masterHeader);
					SerialFile->masterHeader = (char *)malloc(mhsize);
					fseek(fileHandle, 0, SEEK_SET);
					fread(SerialFile->masterHeader,1,mhsize,fileHandle);
				}
				//TODO: check if this is a valid int??
				MasterHeaderSize = mhsize;
			}
			else { // else it's version 0's format w/o version tag, implicating MHSize=4M
				MasterHeaderSize = DefaultMHSize;
				//printf("-----> version = 0; mhsize = %d\n", MasterHeaderSize);
			}

			// END OF CHANGE FOR VERSION
			//
			memcpy( read_out_tag,
					SerialFile->masterHeader+MAX_FIELDS_NAME_LENGTH+1,
					MAX_FIELDS_NAME_LENGTH );

			// Read in # fields ...
			token = strtok ( read_out_tag, ":" );
			token = strtok( NULL," ,;<>" );
			*nfields = atoi( token );
			SerialFile->nfields=*nfields;

			memcpy( read_out_tag,
					SerialFile->masterHeader+
					*nfields * MAX_FIELDS_NAME_LENGTH +
					MAX_FIELDS_NAME_LENGTH * 2,
					MAX_FIELDS_NAME_LENGTH);

			token = strtok ( read_out_tag, ":" );
			token = strtok( NULL," ,;<>" );
			*nppf = atoi( token );
			SerialFile->nppf=*nppf;
		}
		else
		{
			printf("The file you opened is not new format, please check!\n");
		}

		fclose(fileHandle);
	}
	delete [] fname;

}

void finalizephmpiio( int *fileDescriptor )
{
  //printf("total open time is %lf\n", opentime_total);
	// free master header, offset table [][], and serial file struc
	free( SerialFile->masterHeader);
	int j;
	for ( j = 0; j < SerialFile->nfields; j++ )
	{
		free( SerialFile->offset_table[j] );
	}
	free( SerialFile->offset_table);
	free( SerialFile );
}

char* StrReverse(char* str)
{
	char *temp, *ptr;
	int len, i;

	temp=str;
	for(len=0; *temp !='\0';temp++, len++);

	ptr=(char*)malloc(sizeof(char)*(len+1));

	for(i=len-1; i>=0; i--)
		ptr[len-i-1]=str[i];

	ptr[len]='\0';
	return ptr;
}

void openfile( const char filename[],
		const char mode[],
		int*  fileDescriptor )
{
	//printf("in open(): counter = %ld\n", counter++);

	FILE* file=NULL ;
	*fileDescriptor = 0;
	char* fname = StringStripper( filename );
	char* imode = StringStripper( mode );

	int string_length = strlen( fname );
	char* buffer = (char*) malloc ( string_length+1 );
	strcpy ( buffer, fname );
	buffer[ string_length ] = '\0';

	char* tempbuf = StrReverse(buffer);
	free(buffer);
	buffer = tempbuf;

	//printf("buffer is %s\n",buffer);

	char* st2 = strtok ( buffer, "." );
	//st2 = strtok (NULL, ".");

	//printf("st2 is %s\n",st2);

	string_length = strlen(st2);
	char* buffer2 = (char*)malloc(string_length+1);
	strcpy(buffer2,st2);
	buffer2[string_length]='\0';

	char* tempbuf2 = StrReverse(buffer2);
	free(buffer2);
	buffer2 = tempbuf2;
	//printf("buffer2 is %s\n",buffer2);

	SerialFile->fileID = atoi(buffer2);
	if ( char* p = strpbrk(buffer, "@") )
		*p = '\0';

  startTimer(&start);
	if ( cscompare( "read", imode ) ) file = fopen(fname, "rb" );
	else if( cscompare( "write", imode ) ) file = fopen(fname, "wb" );
	else if( cscompare( "append", imode ) ) file = fopen(fname, "ab" );
  endTimer(&end);
  computeTime(&start, &end);


	if ( !file ){
		fprintf(stderr,"unable to open file : %s\n",fname ) ;
	} else {
		fileArray.push_back( file );
		byte_order.push_back( false );
		header_type.push_back( sizeof(int) );
		*fileDescriptor = fileArray.size();
	}

	////////////////////////////////////////////////
	//unsigned long long **header_table;
	SerialFile->offset_table = ( unsigned long long ** )calloc(SerialFile->nfields,
			sizeof(unsigned long long *));

	int j;
	for ( j = 0; j < SerialFile->nfields; j++ )
	{
		SerialFile->offset_table[j]=( unsigned long long * ) calloc( SerialFile->nppf ,
				sizeof( unsigned long long));
	}

	// Read in the offset table ...
	for ( j = 0; j < SerialFile->nfields; j++ )
	{

		memcpy( SerialFile->offset_table[j],
				SerialFile->masterHeader +
				VERSION_INFO_HEADER_SIZE +
				j * SerialFile->nppf * sizeof(unsigned long long),
				SerialFile->nppf * sizeof(unsigned long long) );

		if(diff_endian) {
			SwapArrayByteOrder(  SerialFile->offset_table[j],
			sizeof(unsigned long long int),
			SerialFile->nppf);
			}
		// Swap byte order if endianess is different ...
		/*if ( PhastaIOActiveFiles[i]->Wrong_Endian )
			{
			SwapArrayByteOrder_( PhastaIOActiveFiles[i]->my_read_table[j],
			sizeof(long long int),
			PhastaIOActiveFiles[i]->nppp );
			}
			*/
	}

	////////////////////////////////////////////////
	delete [] fname;
	delete [] imode;
	//free(fname);
	//free(imode);
	free(buffer);
	free(buffer2);
}

void closefile( int* fileDescriptor,
		const char mode[] )
{
	char* imode = StringStripper( mode );

	if( cscompare( "write", imode )
			|| cscompare( "append", imode ) ) {
		fflush( fileArray[ *fileDescriptor - 1 ] );
	}

	fclose( fileArray[ *fileDescriptor - 1 ] );
	delete [] imode;
}

void readheader( int* fileDescriptor,
		const char keyphrase[],
		void* valueArray,
		int*  nItems,
		const char  datatype[],
		const char  iotype[] )
{
	int filePtr = *fileDescriptor - 1;
	FILE* fileObject;
	int* valueListInt;

	if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
		fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
		fprintf(stderr,"openfile_ function has to be called before \n") ;
		fprintf(stderr,"acessing the file\n ") ;
		fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
		return;
	}

	LastHeaderKey[ filePtr ] = const_cast< char* >( keyphrase );
	LastHeaderNotFound = false;

	fileObject = fileArray[ filePtr ] ;
	Wrong_Endian = byte_order[ filePtr ];

	isBinary( iotype );
	typeSize( datatype );   //redundant call, just avoid a compiler warning.

	// right now we are making the assumption that we will only write integers
	// on the header line.

	valueListInt = static_cast< int* >( valueArray );

	/////////////////////////////////////////////////////////
	int j;
	bool FOUND = false ;
	unsigned int skip_size;
	char * token;
	char readouttag[MAX_FIELDS_NUMBER][MAX_FIELDS_NAME_LENGTH];

	int string_length = strlen( keyphrase );
	char* buffer = (char*) malloc ( string_length+1 );
	strcpy ( buffer, keyphrase );
	buffer[ string_length ] = '\0';

	char* st2 = strtok ( buffer, "@" );
	st2 = strtok (NULL, "@");
	SerialFile->GPid = atoi(st2);
	if ( char* p = strpbrk(buffer, "@") )
		*p = '\0';

	//printf("field is %s and nfields is %d\n",keyphrase,SerialFile->nfields);

	for ( j = 0; j<SerialFile->nfields; j++ )
	{
		memcpy( readouttag[j],
				SerialFile->masterHeader + j*MAX_FIELDS_NAME_LENGTH+MAX_FIELDS_NAME_LENGTH*2+1,
				MAX_FIELDS_NAME_LENGTH-1 );
	}

	for ( j = 0; j<SerialFile->nfields; j++ )
	{
		token = strtok ( readouttag[j], ":" );

		if ( cscompare( buffer, token ) )
		{
			SerialFile->read_field_count = j;
			FOUND = true;
			break;
		}
	}
	if (!FOUND)
	{
		printf("Not found %s \n",keyphrase);
		return;
	}

	int read_part_count =  SerialFile->GPid - ( SerialFile->fileID - 1 ) * SerialFile->nppf - 1;
	SerialFile->my_offset = SerialFile->offset_table[SerialFile->read_field_count][read_part_count];


	//printf("GP id is %d and fileID is %d and nppf is %d; ",SerialFile->GPid,SerialFile->fileID,SerialFile->nppf);
	//printf("read field count is %d and read part count is %d; ",SerialFile->read_field_count,read_part_count);

	char read_out_header[MAX_FIELDS_NAME_LENGTH];
	fseek(fileObject, SerialFile->my_offset+1, SEEK_SET);
	fread( read_out_header, 1, MAX_FIELDS_NAME_LENGTH-1, fileObject );


	token = strtok ( read_out_header, ":" );

	if( cscompare( keyphrase , token ) )
	{
		FOUND = true ;
		token = strtok( NULL, " ,;<>" );
		skip_size = atoi( token );
		for( j=0; j < *nItems && ( token = strtok( NULL," ,;<>") ); j++ )
			valueListInt[j] = atoi( token );
		//printf("$$Keyphrase is %s Value list [0] is %d \n",keyphrase,valueListInt[0] );
		if ( j < *nItems )
		{
			fprintf( stderr, "Expected # of ints not found for: %s\n", keyphrase );
		}
	}

	/////////////////////////////////////////////////////////

	byte_order[ filePtr ] = Wrong_Endian ;

	//if ( ierr ) LastHeaderNotFound = true;

	free(buffer);

	return;
}


void readdatablock( int*  fileDescriptor,
		const char keyphrase[],
		void* valueArray,
		int*  nItems,
		const char  datatype[],
		const char  iotype[] )
{

	int filePtr = *fileDescriptor - 1;
	FILE* fileObject;

	if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
		fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
		fprintf(stderr,"openfile_ function has to be called before \n") ;
		fprintf(stderr,"acessing the file\n ") ;
		fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
		return;
	}

	// error check..
	// since we require that a consistant header always preceed the data block
	// let us check to see that it is actually the case.

	if ( ! cscompare( LastHeaderKey[ filePtr ], keyphrase ) ) {
		fprintf(stderr, "Header not consistant with data block\n");
		fprintf(stderr, "Header: %s\n", LastHeaderKey[ filePtr ] );
		fprintf(stderr, "DataBlock: %s\n ", keyphrase );
		fprintf(stderr, "Please recheck read sequence \n");
		if( Strict_Error ) {
			fprintf(stderr, "fatal error: cannot continue, returning out of call\n");
			return;
		}
	}

	if ( LastHeaderNotFound ) return;

	fileObject = fileArray[ filePtr ];
	Wrong_Endian = byte_order[ filePtr ];
	//printf("in readdatablock(): wrong_endian = %d\n", Wrong_Endian);

	size_t type_size = typeSize( datatype );
	int nUnits = *nItems;
	isBinary( iotype );

	if ( binary_format ) {
		fseek(fileObject, SerialFile->my_offset+DB_HEADER_SIZE, SEEK_SET);

		fread( valueArray, type_size, nUnits, fileObject );
		//fread( &junk, sizeof(char), 1 , fileObject );
		//if ( Wrong_Endian ) SwapArrayByteOrder_( valueArray, type_size, nUnits );
		if ( diff_endian )
      SwapArrayByteOrder( valueArray, type_size, nUnits ); // fj
	} else {

		char* ts1 = StringStripper( datatype );
		if ( cscompare( "integer", ts1 ) ) {
			for( int n=0; n < nUnits ; n++ )
				fscanf(fileObject, "%d\n",(int*)((int*)valueArray+n) );
		} else if ( cscompare( "double", ts1 ) ) {
			for( int n=0; n < nUnits ; n++ )
				fscanf(fileObject, "%lf\n",(double*)((double*)valueArray+n) );
		}
		delete [] ts1;
	}
	return;
}
