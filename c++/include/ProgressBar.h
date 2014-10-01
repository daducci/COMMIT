#ifndef __PROGRESSBAR_H__
#define __PROGRESSBAR_H__

#include <string>
#include <cmath>
#include <cstdio>


class ProgressBar
{
	private:
        string          progress, indent;
        int             counter, size;
        float           step, ratio, nextRatio;

	public:
        void            reset();
        void            reset( unsigned int _size );
        void            inc();
        void            close();

        ProgressBar( unsigned int _size, unsigned int _indent );
        ~ProgressBar();
};


ProgressBar::ProgressBar( unsigned int _size, unsigned int _indent = 0 )
{
    size = _size;
    step = _size / 100.0;
    indent = string( _indent, ' ' );
    reset();
}

ProgressBar::~ProgressBar()
{
};


/* ******************** */
/* Set the bar to empty */
/* ******************** */
void ProgressBar::reset()
{
    counter   = 0;
    ratio     = 0;
    nextRatio = 1;
    progress  = string( 100, ' ' );
    printf( "%s[%s] %5.1f%%\r", indent.c_str(), ratio, progress.c_str() );
    fflush( stdout );
}

void ProgressBar::reset( unsigned int _size )
{
    size = size;
    reset();
}


/* ***************************************************** */
/* Increment the counte and, if needed, the progress bar */
/* ***************************************************** */
void ProgressBar::inc()
{
    if ( counter>size ) return;
    counter++;

    ratio = 100.0*(float)counter/(float)size;
    if ( ratio >= nextRatio || counter == size )
    {
        nextRatio = floor(ratio);
        progress  = string(nextRatio,'=') + string(100-nextRatio,' ');
        printf( "%s[%s] %5.1f%%\r", indent.c_str(), ratio, progress.c_str() );
        fflush( stdout );
    }
}


/* **************************************** */
/* Fill the progress bar and go to new line */
/* **************************************** */
void ProgressBar::close()
{
    counter  = size;
    ratio    = 100.0;
    progress = string( 100, '=' );
    printf( "%s[%s] %5.1f%%\n", indent.c_str(), ratio, progress.c_str() );
    fflush( stdout );
}


#endif
