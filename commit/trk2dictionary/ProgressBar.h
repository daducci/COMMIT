#ifndef __PROGRESSBAR_H__
#define __PROGRESSBAR_H__

#include <string>
#include <cmath>
#include <cstdio>

#define EXTRA_ERASE_CHARS 50

class ProgressBar
{
    private:
        int           width, i, N;
        std::string   prefix, msg;

    public:
        void    reset( unsigned int _N );
        void    setPrefix( std::string _prefix );
        void    inc();
        void    close();

        ProgressBar( unsigned int _N, unsigned int _width );
        ~ProgressBar();
};


ProgressBar::ProgressBar( unsigned int _N, unsigned int _width = 25 )
{
    width = _width;
    reset( _N );
}

ProgressBar::~ProgressBar()
{
};


void ProgressBar::reset( unsigned int _N )
{
    i     = 1;
    N     = _N;
    msg   = "";
}

void ProgressBar::setPrefix( std::string _prefix )
{
    prefix = _prefix;
}

/* ****************************************************** */
/* Increment the counter and, if needed, the progress bar */
/* ****************************************************** */
void ProgressBar::inc()
{
    if ( i < 1 || i > N ) return;

    if ( i % width == 0 || i==N )
    {
        int p = floor( float(i)/N*width );
        msg = prefix + "[" + std::string(p,'=') + std::string(width-p,' ') + "]";
        printf( "%s\r", msg.c_str() );
        fflush( stdout );
    }
    i++;
}


/* **************************************** */
/* Fill the progress bar and go to new line */
/* **************************************** */
void ProgressBar::close()
{
    msg = std::string(prefix.length()+width+2+EXTRA_ERASE_CHARS,' ');
    printf("%s\r",msg.c_str());
    fflush( stdout );
    i    = 0;
    N    = 0;
    msg  = "";
}


#endif
