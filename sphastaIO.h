/*  Primary interface for the Phasta Binary read and write routines these*/
/*  functions are 'C' callable.( All arguments have been kept as pointers to*/
/*  facilitate calling from Fortran )*/
/*  Anil Kumar Karanam  Spring 2003*/
#ifndef _SPHASTAIO_H_
#define _SPHASTAIO_H_

#include <FCMangle.h>

#define queryphmpiio FortranCInterface_GLOBAL_(queryphmpiio,QUERYPHMPIIO)
#define initphmpiio FortranCInterface_GLOBAL_(initphmpiio,INITPHMPIIO)
#define finalizephmpiio FortranCInterface_GLOBAL_(finalizephmpiio,FINALIZEPHMPIIO)

#define openfile FortranCInterface_GLOBAL_(openfile, OPENFILE)
#define closefile FortranCInterface_GLOBAL_(closefile, CLOSEFILE)
#define readheader FortranCInterface_GLOBAL_(readheader, READHEADER)
#define readdatablock FortranCInterface_GLOBAL_(readdatablock, READDATABLOCK)
#define writeheader FortranCInterface_GLOBAL_(writeheader, WRITEHEADER)
#define writedatablock FortranCInterface_GLOBAL_(writedatablock, WRITEDATABLOCK)
#define writestring FortranCInterface_GLOBAL_(writestring, WRITESTRING)
#define togglestrictmode FortranCInterface_GLOBAL_(togglestrictmode, TOGGLESTRICTMODE)
#define SwapArrayByteOrder FortranCInterface_GLOBAL_(swaparraybyteorder, SWAPARRAYBYTEORDER)
#define isLittleEndian FortranCInterface_GLOBAL_(islittleendian, ISLITTLEENDIAN)


#if defined (__cplusplus)
extern "C" {
#endif

  void mem_alloc( void* p, char* type, int size );

  void
  queryphmpiio( const char filename[],
		 int *nfields,
		 int *nppf );

  int
  initphmpiio( int *nfields,
		int *nppf,
		int *nfiles,
		int *filehandle,
		const char mode[] );

  void
  finalizephmpiio( int *fileDescriptor );

    void
    SwapArrayByteOrder( void* array,
                         int   nbytes,
                         int   nItems ) ;
    void

    openfile( const char filename[],
               const char mode[],
               int* fileDescriptor );

    void
    closefile( int* fileDescriptor,
                const char mode[] );

    void
    readheader( int*   fileDescriptor,
                 const char  keyphrase[],
                 void*  valueArray,
                 int*   nItems,
                 const char   datatype[],
                 const char   iotype[] );

    void
    writeheader( const int*  fileDescriptor,
                  const char keyphrase[],
                  const void* valueArray,
                  const int*  nItems,
                  const int*  ndataItems,
                  const char  datatype[],
                  const char  iotype[] );

    void
    readdatablock( int*  fileDescriptor,
                    const char keyphrase[],
                    void* valueArray,
                    int*  nItems,
                    const char  datatype[],
                    const char  iotype[] );


    void
    writedatablock( const int*   fileDescriptor,
                     const char  keyphrase[],
                     const void*  valueArray,
                     const int*   nItems,
                     const char   datatype[],
                     const char   iotype[]  );

    void
    writestring( int* fileDescriptor,
                  const char inString[] );

    void
    togglestrictmode( );

  int
  isLittleEndian( ) ;

#ifdef __cplusplus
} // end of extern "C".
#endif
#endif
