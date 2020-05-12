#define GL_GLEXT_PROTOTYPES 1
#ifdef __APPLE__
    #include <OpenGL/gl.h>
    #include <OpenGL/glext.h>
    #include <GLUT/glut.h>
#else
    #include <GL/gl.h>
    #include <GL/glext.h>
    #include <GL/glut.h>
#endif

#include "OPENGL_utils.h"
using namespace OPENGL_utils;

/* global variables */
GLfloat			id[16], rot[16], rot1[16], rot2[16], rot3[16];
Vec3Df			translation;
Vec3Di			start;
GLint			moving;
GLfloat			zoom;

float ScreenX, ScreenY;

void drawString( const char *string )
{
    static int y = glutGet( GLUT_WINDOW_HEIGHT ) - 50;
    if ( string=="" )
        y = glutGet( GLUT_WINDOW_HEIGHT ) - 50;
    else
    {
        glRasterPos2i(10, y);
        for (const char* c=string; *c != '\0'; c++) 
            glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *c);
        y -= 18;
    }
}

void PrintConfig()
{
    if ( !showConfig )
        return;

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();             
    glLoadIdentity();
    glMatrixMode( GL_MODELVIEW ) ;
    glPushMatrix() ;
    glLoadIdentity() ;
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    glOrtho( 0, w, 0, h, -1, 1 );
    glDisable( GL_DEPTH_TEST ); 

    char s[1024];
    glColor3f(1, 1, 0);
    drawString( "" ); // reset initial position

    drawString( "MAP" );
    sprintf( s, "   - value(%d,%d,%d) = %.2f", VOXEL.x, VOXEL.y, VOXEL.z, MAP(VOXEL.x, VOXEL.y, VOXEL.z) );
    drawString( s );
    sprintf( s, "   - range = [ %.1f ... %.1f ]", MAP_min_view, MAP_max_view );
    drawString( s );
    sprintf( s, "   - opacity = %.1f", MAP_opacity );
    drawString( s );

    drawString( "SIGNAL" );
    sprintf( s, "   - shell = %d/%d  (b=%.1f)", GLYPHS_shell+1, SCHEME_shells_b.size(), SCHEME_shells_b[GLYPHS_shell] );
    drawString( s );
    sprintf( s, "   - use affine = %s", GLYPHS_use_affine?"true":"false" );
    drawString( s );
    sprintf( s, "   - flip = [ %d, %d, %d ]", GLYPHS_flip[0], GLYPHS_flip[1], GLYPHS_flip[2] );
    drawString( s );
    sprintf( s, "   - b0 thr = %.1f", GLYPHS_b0_thr );
    drawString( s );

    if ( PEAKS_n>0 )
    {
        drawString( "PEAKS" );
        sprintf( s, "   - use affine = %s", PEAKS_use_affine?"true":"false" );
        drawString( s );
        sprintf( s, "   - flip = [ %d, %d, %d ]", PEAKS_flip[0], PEAKS_flip[1], PEAKS_flip[2] );
        drawString( s );
        sprintf( s, "   - thr = %.1f", PEAKS_thr );
        drawString( s );
        sprintf( s, "   - normalize = %s", PEAKS_doNormalize?"true":"false" );
        drawString( s );
    }

    if ( TRK_nTractsPlotted>0 )
    {
        drawString( "FIBERS" );
        sprintf( s, "   - shift = [ %.1f %.1f %.1f ]  (voxels)", TRK_offset.x, TRK_offset.y, TRK_offset.z );
        drawString( s );
        sprintf( s, "   - slab thickness = %.1f  (voxels)", TRK_crop );
        drawString( s );
    }

    glEnable (GL_DEPTH_TEST);     
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}


// KEYBOARD callback
// -----------------
void GLUT__keyboard( unsigned char key, GLint x=0, GLint y=0 )
{
    bool doRedraw = true;

    switch( key )
    {
        case 'l': showConfig = 1 - showConfig; break;

        case '1': showPlane[0] = 1 - showPlane[0]; break;
        case '2': showPlane[1] = 1 - showPlane[1]; break;
        case '3': showPlane[2] = 1 - showPlane[2]; break;
        case '4':
            showPlane[0] = 1;
            showPlane[1] = 0;
            showPlane[2] = 0;
            translation.x	= translation.y = 0;
            OPENGL_utils::identity(rot1);
            OPENGL_utils::rotateX(rot1, 90.0, rot2);
            OPENGL_utils::rotateZ(rot2, 90.0, rot);
            break;
        case '5':
            showPlane[0] = 0;
            showPlane[1] = 1;
            showPlane[2] = 0;
            translation.x	= translation.y = 0;
            OPENGL_utils::identity(rot1);
            OPENGL_utils::rotateX(rot1, 90.0, rot);
            break;
        case '6':
            showPlane[0] = 0;
            showPlane[1] = 0;
            showPlane[2] = 1;
            translation.x	= translation.y = 0;
            OPENGL_utils::identity( rot );
            break;

        case '0': showAxes = 1 - showAxes; break;
        case '-': zoom += 10.0; break;
        case '+': zoom -= 10.0; break;
        case 'm': MAP_max_view = fmaxf(0.0,MAP_max_view-MAP_max*0.05); break;
        case 'M': MAP_max_view = fminf(MAP_max,MAP_max_view+MAP_max*0.05); break;
        case 'o': MAP_opacity = fmaxf(0.0,MAP_opacity-0.1); break;
        case 'O': MAP_opacity = fminf(1.0,MAP_opacity+0.1); break;
        case 'w': LINE_width = fmaxf( 1,LINE_width-1); break;
        case 'W': LINE_width = fminf(10,LINE_width+1); break;
        case 'r':
            showPlane[0] = showPlane[1] = showPlane[2] = 1;
            translation.x	= translation.y = 0;
            zoom			= 0;
            OPENGL_utils::identity( rot );
            break;

        case 's': GLYPHS_show = 1 - GLYPHS_show; break;
        case 'S': GLYPHS_shell = (GLYPHS_shell+1) % SCHEME_shells_idx.size(); break;
        case 'a': GLYPHS_use_affine = 1 - GLYPHS_use_affine; break;
        case 'x': GLYPHS_flip[0] = 1 - GLYPHS_flip[0]; for(int d=0; d < SCHEME_dirs.size() ;d++) SCHEME_dirs[d].x *= -1; break;
        case 'y': GLYPHS_flip[1] = 1 - GLYPHS_flip[1]; for(int d=0; d < SCHEME_dirs.size() ;d++) SCHEME_dirs[d].y *= -1; break;
        case 'z': GLYPHS_flip[2] = 1 - GLYPHS_flip[2]; for(int d=0; d < SCHEME_dirs.size() ;d++) SCHEME_dirs[d].z *= -1; break;
        case 'b': GLYPHS_b0_thr = fmaxf(0.0,GLYPHS_b0_thr-10.0); break;
        case 'B': GLYPHS_b0_thr = fminf(MAP_max,GLYPHS_b0_thr+10.0); break;

        case 'p': if ( PEAKS_n>0 ) PEAKS_show  = 1 - PEAKS_show; break;
        case 'A': PEAKS_use_affine = 1 - PEAKS_use_affine; break;
        case 'X': PEAKS_flip[0] = 1 - PEAKS_flip[0]; break;
        case 'Y': PEAKS_flip[1] = 1 - PEAKS_flip[1]; break;
        case 'Z': PEAKS_flip[2] = 1 - PEAKS_flip[2]; break;
        case 't': PEAKS_thr = fmaxf(PEAKS_thr - 0.1, 0.0); break;
        case 'T': PEAKS_thr = fminf(PEAKS_thr + 0.1, 1.0); break;
        case 'n': PEAKS_doNormalize = 1 - PEAKS_doNormalize; break;

        case 'f': if ( TRK_nTractsPlotted>0 ) TRK_show = 1 - TRK_show; break;
        case 'c': TRK_crop = fmaxf( 0.0,TRK_crop-0.5); break;
        case 'C': TRK_crop = fminf(max(dim.x,max(dim.y,dim.z)),TRK_crop+0.5); break;
        case ' ': TRK_crop_mode = 1 - TRK_crop_mode; break;

        case 'q':
        case 27 : exit(0); break;

        default: doRedraw = false;
    }

    if ( doRedraw )
        glutPostRedisplay();
}


// MENU callback
// -------------
void GLUT__menu( int id ) 
{
    switch( id )
    {
        case   0: GLUT__keyboard('q'); break;

        case 101: GLUT__keyboard('s'); break;
        case 102: GLUT__keyboard('S'); break;
        case 103: GLUT__keyboard('a'); break;
        case 104: GLUT__keyboard('x'); break;
        case 105: GLUT__keyboard('y'); break;
        case 106: GLUT__keyboard('z'); break;
        case 107: GLUT__keyboard('b'); break;
        case 108: GLUT__keyboard('B'); break;

        case 201: GLUT__keyboard('p'); break;
        case 202: GLUT__keyboard('A'); break;
        case 203: GLUT__keyboard('X'); break;
        case 204: GLUT__keyboard('Y'); break;
        case 205: GLUT__keyboard('Z'); break;
        case 206: GLUT__keyboard('t'); break;
        case 207: GLUT__keyboard('T'); break;
        case 208: GLUT__keyboard('n'); break;

        case 301: GLUT__keyboard('f'); break;
        case 302: GLUT__keyboard('c'); break;
        case 303: GLUT__keyboard('C'); break;
        case 304: GLUT__keyboard(' '); break;

        case 401: GLUT__keyboard('1'); break;
        case 402: GLUT__keyboard('2'); break;
        case 403: GLUT__keyboard('3'); break;
        case 404: GLUT__keyboard('4'); break;
        case 405: GLUT__keyboard('5'); break;
        case 406: GLUT__keyboard('6'); break;
        case 407: GLUT__keyboard('0'); break;
        case 408: GLUT__keyboard('-'); break;
        case 409: GLUT__keyboard('+'); break;
        case 410: GLUT__keyboard('m'); break;
        case 411: GLUT__keyboard('M'); break;
        case 412: GLUT__keyboard('o'); break;
        case 413: GLUT__keyboard('O'); break;
        case 414: GLUT__keyboard('w'); break;
        case 415: GLUT__keyboard('W'); break;
        case 416: GLUT__keyboard('r'); break;
        case 417: GLUT__keyboard('l'); break;
    }
}


// Create the dropdown MENU
// ------------------------
void GLUT__createMenu()
{
    int submenu_SIGNAL_id, submenu_PEAKS_id, submenu_FIBERS_id, submenu_VIEW_id;

    submenu_SIGNAL_id = glutCreateMenu( GLUT__menu );
    glutAddMenuEntry("[s] Show/hide",         101);
    glutAddMenuEntry("[S] Change shell",      102);
    glutAddMenuEntry("[a] Use affine",        103);
    glutAddMenuEntry("[x] Flip X axis",       104);
    glutAddMenuEntry("[y] Flip Y axis",       105);
    glutAddMenuEntry("[z] Flip Z axis",       106);
    glutAddMenuEntry("[b] Decrease b0 thr",   107);
    glutAddMenuEntry("[B] Increase b0 thr",   108);

    if ( PEAKS_n>0 )
    {
        submenu_PEAKS_id = glutCreateMenu( GLUT__menu );
        glutAddMenuEntry("[p] Show/hide",         201);
        glutAddMenuEntry("[A] Use affine",        202);
        glutAddMenuEntry("[X] Flip X axis",       203);
        glutAddMenuEntry("[Y] Flip Y axis",       204);
        glutAddMenuEntry("[Z] Flip Z axis",       205);
        glutAddMenuEntry("[t] Decrease threshold",206);
        glutAddMenuEntry("[T] Increase threshold",207);
        glutAddMenuEntry("[n] Normalize length",  208);
    }

    if ( TRK_nTractsPlotted>0 )
    {
        submenu_FIBERS_id = glutCreateMenu( GLUT__menu );
        glutAddMenuEntry("[f] Show/hide",         301);
        glutAddMenuEntry("[c] Decrease crop size",302);
        glutAddMenuEntry("[C] Increase crop size",303);
        glutAddMenuEntry("[ ] Change crop mode",  304);
    }

    submenu_VIEW_id = glutCreateMenu( GLUT__menu );
    glutAddMenuEntry("[1] Show/hide YZ plane", 401);
    glutAddMenuEntry("[2] Show/hide XZ plane", 402);
    glutAddMenuEntry("[3] Show/hide XY plane", 403);
    glutAddMenuEntry("[4] Reset to YZ plane",  404);
    glutAddMenuEntry("[5] Reset to XZ plane",  405);
    glutAddMenuEntry("[6] Reset to XY plane",  406);
    glutAddMenuEntry("[0] Show/hide axes",     407);
    glutAddMenuEntry("[-] Decrease zoom",      408);
    glutAddMenuEntry("[+] Increase zoom",      409);
    glutAddMenuEntry("[m] Decrease max value", 410);
    glutAddMenuEntry("[M] Increase max value", 411);
    glutAddMenuEntry("[o] Decrease opacity",   412);
    glutAddMenuEntry("[O] Increase opacity",   413);
    glutAddMenuEntry("[t] Decrease line width",414);
    glutAddMenuEntry("[T] Increase line width",415);
    glutAddMenuEntry("[r] Reset view",         416);
    glutAddMenuEntry("[l] Show/hide log",      417);

    int menu_id = glutCreateMenu( GLUT__menu );
    glutAddSubMenu("Signal", submenu_SIGNAL_id);
    if ( PEAKS_n>0 )
        glutAddSubMenu("Peaks", submenu_PEAKS_id);
    if ( TRK_nTractsPlotted>0 )
        glutAddSubMenu("Fibers", submenu_FIBERS_id);
    glutAddSubMenu("View options", submenu_VIEW_id);
    glutAddMenuEntry("Quit", 0);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}


// RESHAPE callback
// ----------------
void GLUT__reshape( GLint w, GLint h )
{
    ScreenX = w;
    ScreenY = h;

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 45.0f, (GLfloat)w / (GLfloat)h, 1.0f, 5000.0f );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt(
        0.0, 0.0, 2.0 * max(pixdim.x*dim.x,pixdim.y*dim.y) * (GLfloat)ScreenY/(GLfloat)ScreenX, // eye point
        0.0, 0.0, 0.0, // reference point
        0.0, 1.0, 0.0  // up vector
    );
}


// SPECIALKEY callback
// -------------------
void GLUT__specialkey( GLint key, GLint x, GLint y )
{
    bool doRedraw = true;
    GLint modif = glutGetModifiers();
    GLint ALT   = modif & GLUT_ACTIVE_ALT;
    GLint CTRL  = modif & GLUT_ACTIVE_CTRL;

    switch( key )
    {
        case GLUT_KEY_LEFT:
            if ( ALT )
                TRK_offset.x -= 0.5;
            else if ( CTRL )
                translation.x -= 2.0;
            else
                VOXEL.x--;
            break;
        case GLUT_KEY_RIGHT:
            if ( ALT )
                TRK_offset.x += 0.5;
            else if ( CTRL )
                translation.x += 2.0;
            else
                VOXEL.x++;
            break;
        case GLUT_KEY_DOWN:
            if ( ALT )
                TRK_offset.y -= 0.5;
            else if ( CTRL )
                translation.y -= 2.0;
            else
                VOXEL.y--;
            break;
        case GLUT_KEY_UP:
            if ( ALT )
                TRK_offset.y += 0.5;
            else if ( CTRL )
                translation.y += 2.0;
            else
                VOXEL.y++;
            break;
        case GLUT_KEY_PAGE_DOWN:
            if ( ALT )
                TRK_offset.z -= 0.5;
            else
                VOXEL.z--;
            break;
        case GLUT_KEY_PAGE_UP:
            if ( ALT )
                TRK_offset.z += 0.5;
            else
                VOXEL.z++;
            break;

        default:
            doRedraw = false;
    }

    // check the bounds
    VOXEL.x = max( VOXEL.x, 0 );
    VOXEL.y = max( VOXEL.y, 0 );
    VOXEL.z = max( VOXEL.z, 0 );
    VOXEL.x = min( VOXEL.x, dim.x-1 );
    VOXEL.y = min( VOXEL.y, dim.y-1 );
    VOXEL.z = min( VOXEL.z, dim.z-1 );

    if ( doRedraw )
        glutPostRedisplay();
}



// MOUSE callback
// --------------
void GLUT__mouse( GLint button, GLint state, GLint x, GLint y )
{
    if (state == GLUT_DOWN)
    {
        if ( button == GLUT_LEFT_BUTTON && glutGetModifiers() != GLUT_ACTIVE_CTRL )
        {
            moving = 1;
            start.x = x;
            start.y = y;
        }
        // NOTE: does not work, issue with glutGetModifiers not getting CTRL
        // else if ( button == GLUT_LEFT_BUTTON && glutGetModifiers() == GLUT_ACTIVE_CTRL )
        // {
        //     moving = 2;
        //     start.x = x;
        //     start.y = y;
        // }
        else if ( (button == GLUT_MIDDLE_BUTTON) || (button == GLUT_LEFT_BUTTON && glutGetModifiers() == GLUT_ACTIVE_ALT) )
        {
            moving = 3;
            start.x = x;
            start.y = y;
        }
    }
    else if (state == GLUT_UP)
    {
        moving = 0;
    }
}


// MOTION callback
// ---------------
void GLUT__motion( GLint x, GLint y )
{
    if (moving==1)
    {
        OPENGL_utils::translate(id, 0,0,0, rot1);

        OPENGL_utils::rotateY(id,start.x-x,rot3);
        OPENGL_utils::matXMat(rot,rot1,rot2);
        OPENGL_utils::rotateX(id,start.y-y,rot1);
        OPENGL_utils::matXMat(rot2,rot1,rot);
        OPENGL_utils::matXMat(rot,rot3,rot2);

        OPENGL_utils::translate(id, 0,0,0, rot1);
        OPENGL_utils::matXMat(rot2,rot1,rot);

        start.x = x;
        start.y = y;
    }

    else if (moving==2)
    {
        zoom = zoom + (y-start.y)/2.0;
        start.y = y;
    }

    else if (moving==3)
    {
        translation.x = translation.x - (start.x-x)/3.0;
        translation.y = translation.y + (start.y-y)/3.0;
        start.x = x;
        start.y = y;
    }

    glutPostRedisplay();
}


// DISPLAY callback
// ----------------
void GLUT__display( void )
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glPushMatrix();
    glTranslatef(translation.x, translation.y, -zoom); // mouse translation + zoom
    glMultMatrixf(rot); // mouse rotation    
    glTranslatef( -pixdim.x*dim.x/2.0, -pixdim.y*dim.y/2.0, -pixdim.z*dim.z/2.0 ); // center the FOV
    glScalef( pixdim.x, pixdim.y, pixdim.z ); // account for voxel size

    /* ============= */
    /* Draw the AXES */
    /* ============= */
    if ( showAxes )
    {
        glLineWidth(2);
        glBegin(GL_LINES);
            glColor4f( 1,0,0,1); glVertex3f( 0,0,0 ); glVertex3f( 10,  0,  0 );
            glColor4f( 0,1,0,1); glVertex3f( 0,0,0 ); glVertex3f(  0, 10,  0 );
            glColor4f( 0,0,1,1); glVertex3f( 0,0,0 ); glVertex3f(  0,  0, 10 );
        glEnd();
    }

    /* =============== */
    /* Draw the TRACTS */
    /* =============== */
    if ( TRK_show )
    {
        glPushMatrix();
        glTranslatef(TRK_offset.x, TRK_offset.y, TRK_offset.z);

        glLineWidth(1.0f);

        float *ptr  = TRK_coords, *ptrc = TRK_colors;
        VECTOR<float> Vc( VOXEL.x+0.5, VOXEL.y+0.5, VOXEL.z+0.5 ); // voxel center
        float thr = 0.5*TRK_crop;
        for(int f=0; f < TRK_nTractsPlotted; f++)
        {
            glBegin(GL_LINE_STRIP);
            for(int i=0; i < TRK_nPoints[f]; i++)
            {
                // plot segment only if it's close to center of VOXEL
                if (
                      (
                        TRK_crop_mode && (
                        ( showPlane[0] && abs( (ptr[0]+TRK_offset.x) - Vc.x ) <= thr ) ||
                        ( showPlane[1] && abs( (ptr[1]+TRK_offset.y) - Vc.y ) <= thr ) ||
                        ( showPlane[2] && abs( (ptr[2]+TRK_offset.z) - Vc.z ) <= thr ) )
                      )
                      ||
                      (
                        !TRK_crop_mode && (
                        ( abs( (ptr[0]+TRK_offset.x) - Vc.x ) <= thr ) &&
                        ( abs( (ptr[1]+TRK_offset.y) - Vc.y ) <= thr ) &&
                        ( abs( (ptr[2]+TRK_offset.z) - Vc.z ) <= thr ) )
                      )
                    )
                {
                    glColor3f(  *ptrc++, *ptrc++, *ptrc++ );
                    glVertex3f( *ptr++,  *ptr++,  *ptr++  );
                }
                else
                {
                    glEnd();
                    glBegin(GL_LINE_STRIP);
                    ptr  += 3;
                    ptrc += 3;
                }
            }
            glEnd();
        }

        glPopMatrix();
    }

    /* ============== */
    /* Draw the PEAKS */
    /* ============== */
    if ( PEAKS_show || GLYPHS_show )
    {
        glDisable( GL_BLEND );
        glLineWidth( LINE_width );
        glPointSize( LINE_width );

        glPushMatrix();
        glTranslatef(.5,.5,.5);

        Vec3Df dir, col;
        int x,y,z,d,idx;
        float norms[PEAKS_n], normMax, b0, w;

        // plane YZ
        if ( showPlane[0]  )
        {
            x = (int)VOXEL.x;
            for(y=0; y<dim.y ;y++)
            for(z=0; z<dim.z ;z++)
            {
                if ( PEAKS_show )
                {
                    normMax = 0;
                    for(d=0; d<PEAKS_n; d++)
                    {
                        col.x = (*niiPEAKS->img)(x,y,z,3*d+0); // use "col" as tmp variable
                        col.y = (*niiPEAKS->img)(x,y,z,3*d+1);
                        col.z = (*niiPEAKS->img)(x,y,z,3*d+2);
                        if ( PEAKS_use_affine )
                        {
                            dir.x = col.x * ((float*)PEAKS_affine)[0] + col.y * ((float*)PEAKS_affine)[1] + col.z * ((float*)PEAKS_affine)[2];
                            dir.y = col.x * ((float*)PEAKS_affine)[3] + col.y * ((float*)PEAKS_affine)[4] + col.z * ((float*)PEAKS_affine)[5];
                            dir.z = col.x * ((float*)PEAKS_affine)[6] + col.y * ((float*)PEAKS_affine)[7] + col.z * ((float*)PEAKS_affine)[8];
                        }
                        else
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        norms[d] = dir.norm();
                        if ( norms[d] > normMax )
                            normMax = norms[d];
                    }

                    for(d=0; d<PEAKS_n; d++)
                    {
                        if ( norms[d] < PEAKS_thr*normMax )
                            continue;

                        col.x = (*niiPEAKS->img)(x,y,z,3*d+0); // use "col" as tmp variable
                        col.y = (*niiPEAKS->img)(x,y,z,3*d+1);
                        col.z = (*niiPEAKS->img)(x,y,z,3*d+2);
                        if ( PEAKS_use_affine )
                        {
                            dir.x = col.x * ((float*)PEAKS_affine)[0] + col.y * ((float*)PEAKS_affine)[1] + col.z * ((float*)PEAKS_affine)[2];
                            dir.y = col.x * ((float*)PEAKS_affine)[3] + col.y * ((float*)PEAKS_affine)[4] + col.z * ((float*)PEAKS_affine)[5];
                            dir.z = col.x * ((float*)PEAKS_affine)[6] + col.y * ((float*)PEAKS_affine)[7] + col.z * ((float*)PEAKS_affine)[8];
                        }
                        else
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        col.x = 0.5 * (PEAKS_flip[0]?-1:1) * dir.x / norms[d];
                        col.y = 0.5 * (PEAKS_flip[1]?-1:1) * dir.y / norms[d];
                        col.z = 0.5 * (PEAKS_flip[2]?-1:1) * dir.z / norms[d];

                        if ( PEAKS_doNormalize )
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        else
                        {
                            dir.x = col.x * norms[d] / normMax;
                            dir.y = col.y * norms[d] / normMax;
                            dir.z = col.z * norms[d] / normMax;
                        }

                        glColor3f( fabs(2.0*col.x), fabs(2.0*col.y), fabs(2.0*col.z) ); 
                        glBegin(GL_LINES);
                            glVertex3f( x-dir.x, y-dir.y, z-dir.z );
                            glVertex3f( x+dir.x, y+dir.y, z+dir.z );
                        glEnd();
                    }
                }
                if ( GLYPHS_show )
                {
                    b0 = (*niiDWI->img)(x,y,z,SCHEME_idxB0[0]);
                    if ( b0 > GLYPHS_b0_thr )
                    {
                        glBegin(GL_POINTS);
                        for(d=0; d < SCHEME_shells_idx[GLYPHS_shell].size() ;d++)
                        {
                            idx = SCHEME_shells_idx[GLYPHS_shell][d];
                            w = 0.5 * (float)(*niiDWI->img)(x,y,z,idx) / b0;
                            if ( GLYPHS_use_affine ) 
                            {
                                dir.x = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[0] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[1] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[2];
                                dir.y = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[3] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[4] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[5];
                                dir.z = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[6] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[7] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[8];
                                normMax = dir.norm();
                                dir.x *= w / normMax;
                                dir.y *= w / normMax;
                                dir.z *= w / normMax;
                            }
                            else
                            {
                                dir.x = w * SCHEME_dirs[idx].x;
                                dir.y = w * SCHEME_dirs[idx].y;
                                dir.z = w * SCHEME_dirs[idx].z;
                            }
                            normMax = dir.norm();
                            glColor3f( fabs(dir.x)/normMax, fabs(dir.y)/normMax, fabs(dir.z)/normMax );
                            glVertex3f( x+dir.x, y+dir.y, z+dir.z );
                            glVertex3f( x-dir.x, y-dir.y, z-dir.z );
                        }
                        glEnd();
                    }
                }
            }
        }

        // plane XZ
        if ( showPlane[1] )
        {
            y = (int)VOXEL.y;
            for(x=0; x<dim.x ;x++)
            for(z=0; z<dim.z ;z++)
            {
                if ( PEAKS_show )
                {
                    normMax = 0;
                    for(d=0; d<PEAKS_n; d++)
                    {
                        col.x = (*niiPEAKS->img)(x,y,z,3*d+0); // use "col" as tmp variable
                        col.y = (*niiPEAKS->img)(x,y,z,3*d+1);
                        col.z = (*niiPEAKS->img)(x,y,z,3*d+2);
                        if ( PEAKS_use_affine )
                        {
                            dir.x = col.x * ((float*)PEAKS_affine)[0] + col.y * ((float*)PEAKS_affine)[1] + col.z * ((float*)PEAKS_affine)[2];
                            dir.y = col.x * ((float*)PEAKS_affine)[3] + col.y * ((float*)PEAKS_affine)[4] + col.z * ((float*)PEAKS_affine)[5];
                            dir.z = col.x * ((float*)PEAKS_affine)[6] + col.y * ((float*)PEAKS_affine)[7] + col.z * ((float*)PEAKS_affine)[8];
                        }
                        else
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        norms[d] = dir.norm();
                        if ( norms[d] > normMax )
                            normMax = norms[d];
                    }

                    for(d=0; d<PEAKS_n; d++)
                    {
                        if ( norms[d] < normMax*PEAKS_thr )
                            continue;

                        col.x = (*niiPEAKS->img)(x,y,z,3*d+0); // use "col" as tmp variable
                        col.y = (*niiPEAKS->img)(x,y,z,3*d+1);
                        col.z = (*niiPEAKS->img)(x,y,z,3*d+2);
                        if ( PEAKS_use_affine )
                        {
                            dir.x = col.x * ((float*)PEAKS_affine)[0] + col.y * ((float*)PEAKS_affine)[1] + col.z * ((float*)PEAKS_affine)[2];
                            dir.y = col.x * ((float*)PEAKS_affine)[3] + col.y * ((float*)PEAKS_affine)[4] + col.z * ((float*)PEAKS_affine)[5];
                            dir.z = col.x * ((float*)PEAKS_affine)[6] + col.y * ((float*)PEAKS_affine)[7] + col.z * ((float*)PEAKS_affine)[8];
                        }
                        else
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        col.x = 0.5 * (PEAKS_flip[0]?-1:1) * dir.x / norms[d];
                        col.y = 0.5 * (PEAKS_flip[1]?-1:1) * dir.y / norms[d];
                        col.z = 0.5 * (PEAKS_flip[2]?-1:1) * dir.z / norms[d];

                        if ( PEAKS_doNormalize )
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        else
                        {
                            dir.x = col.x * norms[d] / normMax;
                            dir.y = col.y * norms[d] / normMax;
                            dir.z = col.z * norms[d] / normMax;
                        }

                        glColor3f( fabs(2.0*col.x), fabs(2.0*col.y), fabs(2.0*col.z) );
                        glBegin(GL_LINES);
                            glVertex3f( x-dir.x, y-dir.y, z-dir.z );
                            glVertex3f( x+dir.x, y+dir.y, z+dir.z );
                        glEnd();
                    }
                }

                if ( GLYPHS_show )
                {
                    b0 = (*niiDWI->img)(x,y,z,SCHEME_idxB0[0]);
                    if ( b0 > GLYPHS_b0_thr )
                    {
                        glBegin(GL_POINTS);
                        for(d=0; d < SCHEME_shells_idx[GLYPHS_shell].size() ;d++)
                        {
                            idx = SCHEME_shells_idx[GLYPHS_shell][d];
                            w = 0.5 * (float)(*niiDWI->img)(x,y,z,idx) / b0;
                            if ( GLYPHS_use_affine ) 
                            {
                                dir.x = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[0] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[1] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[2];
                                dir.y = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[3] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[4] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[5];
                                dir.z = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[6] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[7] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[8];
                                normMax = dir.norm();
                                dir.x *= w / normMax;
                                dir.y *= w / normMax;
                                dir.z *= w / normMax;
                            }
                            else
                            {
                                dir.x = w * SCHEME_dirs[idx].x;
                                dir.y = w * SCHEME_dirs[idx].y;
                                dir.z = w * SCHEME_dirs[idx].z;
                            }
                            normMax = dir.norm();
                            glColor3f( fabs(dir.x)/normMax, fabs(dir.y)/normMax, fabs(dir.z)/normMax );
                            glVertex3f( x+dir.x, y+dir.y, z+dir.z );
                            glVertex3f( x-dir.x, y-dir.y, z-dir.z );
                        }
                        glEnd();
                    }
                }
            }
        }

        // plane XY
        if ( showPlane[2] )
        {
            z = (int)VOXEL.z;
            for(y=0; y<dim.y ;y++)
            for(x=0; x<dim.x ;x++)
            {
                if ( PEAKS_show )
                {
                    normMax = 0;
                    for(d=0; d<PEAKS_n; d++)
                    {
                        col.x = (*niiPEAKS->img)(x,y,z,3*d+0); // use "col" as tmp variable
                        col.y = (*niiPEAKS->img)(x,y,z,3*d+1);
                        col.z = (*niiPEAKS->img)(x,y,z,3*d+2);
                        if ( PEAKS_use_affine )
                        {
                            dir.x = col.x * ((float*)PEAKS_affine)[0] + col.y * ((float*)PEAKS_affine)[1] + col.z * ((float*)PEAKS_affine)[2];
                            dir.y = col.x * ((float*)PEAKS_affine)[3] + col.y * ((float*)PEAKS_affine)[4] + col.z * ((float*)PEAKS_affine)[5];
                            dir.z = col.x * ((float*)PEAKS_affine)[6] + col.y * ((float*)PEAKS_affine)[7] + col.z * ((float*)PEAKS_affine)[8];
                        }
                        else
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        norms[d] = dir.norm();
                        if ( norms[d] > normMax )
                            normMax = norms[d];
                    }

                    for(d=0; d<PEAKS_n; d++)
                    {
                        if ( norms[d] < normMax*PEAKS_thr )
                            continue;

                        col.x = (*niiPEAKS->img)(x,y,z,3*d+0); // use "col" as tmp variable
                        col.y = (*niiPEAKS->img)(x,y,z,3*d+1);
                        col.z = (*niiPEAKS->img)(x,y,z,3*d+2);
                        if ( PEAKS_use_affine )
                        {
                            dir.x = col.x * ((float*)PEAKS_affine)[0] + col.y * ((float*)PEAKS_affine)[1] + col.z * ((float*)PEAKS_affine)[2];
                            dir.y = col.x * ((float*)PEAKS_affine)[3] + col.y * ((float*)PEAKS_affine)[4] + col.z * ((float*)PEAKS_affine)[5];
                            dir.z = col.x * ((float*)PEAKS_affine)[6] + col.y * ((float*)PEAKS_affine)[7] + col.z * ((float*)PEAKS_affine)[8];
                        }
                        else
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        col.x = 0.5 * (PEAKS_flip[0]?-1:1) * dir.x / norms[d];
                        col.y = 0.5 * (PEAKS_flip[1]?-1:1) * dir.y / norms[d];
                        col.z = 0.5 * (PEAKS_flip[2]?-1:1) * dir.z / norms[d];

                        if ( PEAKS_doNormalize )
                        {
                            dir.x = col.x;
                            dir.y = col.y;
                            dir.z = col.z;
                        }
                        else
                        {
                            dir.x = col.x * norms[d] / normMax;
                            dir.y = col.y * norms[d] / normMax;
                            dir.z = col.z * norms[d] / normMax;
                        }

                        glColor3f( fabs(2.0*col.x), fabs(2.0*col.y), fabs(2.0*col.z) );
                        glBegin(GL_LINES);
                            glVertex3f( x-dir.x, y-dir.y, z-dir.z );
                            glVertex3f( x+dir.x, y+dir.y, z+dir.z );
                        glEnd();
                    }
                }

                if( GLYPHS_show)
                {
                    b0 = (*niiDWI->img)(x,y,z,SCHEME_idxB0[0]);
                    if ( b0 > GLYPHS_b0_thr )
                    {
                        glBegin(GL_POINTS);
                        for(d=0; d < SCHEME_shells_idx[GLYPHS_shell].size() ;d++)
                        {
                            idx = SCHEME_shells_idx[GLYPHS_shell][d];
                            w = 0.5 * (float)(*niiDWI->img)(x,y,z,idx) / b0;
                            if ( GLYPHS_use_affine ) 
                            {
                                dir.x = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[0] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[1] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[2];
                                dir.y = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[3] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[4] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[5];
                                dir.z = SCHEME_dirs[idx].x * ((float*)GLYPHS_affine)[6] + SCHEME_dirs[idx].y * ((float*)GLYPHS_affine)[7] + SCHEME_dirs[idx].z * ((float*)GLYPHS_affine)[8];
                                normMax = dir.norm();
                                dir.x *= w / normMax;
                                dir.y *= w / normMax;
                                dir.z *= w / normMax;
                            }
                            else
                            {
                                dir.x = w * SCHEME_dirs[idx].x;
                                dir.y = w * SCHEME_dirs[idx].y;
                                dir.z = w * SCHEME_dirs[idx].z;
                            }

                            normMax = dir.norm();
                            glColor3f( fabs(dir.x)/normMax, fabs(dir.y)/normMax, fabs(dir.z)/normMax );
                            glVertex3f( x+dir.x, y+dir.y, z+dir.z );
                            glVertex3f( x-dir.x, y-dir.y, z-dir.z );
                        }
                        glEnd();
                    }
                }
            }
        }

        glPopMatrix();
    }

    /* =================== */
    /* Draw the SCALAR MAP */
    /* =================== */
    if ( showPlane[0] || showPlane[1] || showPlane[2] )
    {
        glDisable( GL_CULL_FACE );
        glEnable( GL_BLEND );
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

        // to avoid z-fighting
        glPolygonOffset( 1.0, 1.0 );
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glLineWidth( 3 );

        int x, y, z; // voxel coordinates NB: (0,0,0) -> corner of voxel
        float color;

        // plane YZ
        if ( showPlane[0]  )
        {
            glPushMatrix();
            glTranslatef(0.5,0,0);

            x = (int)VOXEL.x;
            for(y=0; y<dim.y ;y++)
            for(z=0; z<dim.z ;z++)
            {
                color = ( MAP(x,y,z) - MAP_min_view) / ( MAP_max_view - MAP_min_view );
                glColor4f(color,color,color,MAP_opacity);
                glBegin(GL_QUADS);
                    glVertex3f(x, y,   z);
                    glVertex3f(x, y,   z+1);
                    glVertex3f(x, y+1, z+1);
                    glVertex3f(x, y+1, z);
                glEnd();
            }
            // colored frame
            if ( showAxes )
            {
                glColor3f(1,0,0);
                glBegin(GL_LINE_STRIP);
                    glVertex3f(x,0,0);
                    glVertex3f(x,dim.y,0);
                    glVertex3f(x,dim.y,dim.z);
                    glVertex3f(x,0,dim.z);
                    glVertex3f(x,0,0);
                glEnd();
            }

            glPopMatrix();
        }

        // plane XZ
        if ( showPlane[1] )
        {
            glPushMatrix();
            glTranslatef(0,0.5,0);

            y = (int)VOXEL.y;
            for(x=0; x<dim.x ;x++)
            for(z=0; z<dim.z ;z++)
            {
                color = ( MAP(x,y,z) - MAP_min_view) / ( MAP_max_view - MAP_min_view );
                glColor4f(color,color,color,MAP_opacity);
                glBegin(GL_QUADS);
                    glVertex3f(x,   y, z);
                    glVertex3f(x,   y, z+1);
                    glVertex3f(x+1, y, z+1);
                    glVertex3f(x+1, y, z);
                glEnd();
            }
            // colored frame
            if ( showAxes )
            {
                glColor3f(0,1,0);
                glBegin(GL_LINE_STRIP);
                    glVertex3f(0,y,0);
                    glVertex3f(dim.x,y,0);
                    glVertex3f(dim.x,y,dim.z);
                    glVertex3f(0,y,dim.z);
                    glVertex3f(0,y,0);
                glEnd();
            }

            glPopMatrix();
        }

        // plane XY
        if ( showPlane[2] )
        {
            glPushMatrix();
            glTranslatef(0,0,0.5);

            z = (int)VOXEL.z;
            for(y=0; y<dim.y ;y++)
            for(x=0; x<dim.x ;x++)
            {
                color = ( MAP(x,y,z) - MAP_min_view) / ( MAP_max_view - MAP_min_view );
                glColor4f(color,color,color,MAP_opacity);
                glBegin(GL_QUADS);
                    glVertex3f(x,   y,   z);
                    glVertex3f(x+1, y,   z);
                    glVertex3f(x+1, y+1, z);
                    glVertex3f(x,   y+1, z);
                glEnd();
            }

            // colored frame
            if ( showAxes )
            {
                glColor3f(0,0,1);
                glBegin(GL_LINE_STRIP);
                    glVertex3f(0,0,z);
                    glVertex3f(dim.x,0,z);
                    glVertex3f(dim.x,dim.y,z);
                    glVertex3f(0,dim.y,z);
                    glVertex3f(0,0,z);
                glEnd();
            }

            glPopMatrix();
        }

        glEnable(GL_CULL_FACE);
        glDisable( GL_BLEND );
        glDisable(GL_POLYGON_OFFSET_FILL);
    }

    /* ====================== */
    /* Draw the CURRENT VOXEL */
    /* ====================== */
    if ( showAxes )
    {
        glPushMatrix();
        glTranslatef( VOXEL.x+0.5, VOXEL.y+0.5, VOXEL.z+0.5 );

        glEnable( GL_BLEND );
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        glLineWidth(1);
        glColor4f( 1,1,0,1 );
        glutWireCube( 1 );
        glColor4f( 1,1,0,0.25 );
        glutSolidCube( 1 );
        glDisable( GL_BLEND );

        glPopMatrix();
    }

    glPopMatrix();
    PrintConfig();
    glutSwapBuffers();
}


// INITIALIZATION
// --------------
void OpenGL_init( int argc, char** argv )
{
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA );
    ScreenX = 0.7*glutGet(GLUT_SCREEN_WIDTH);  if (ScreenX==0) ScreenX = 800;
    ScreenY = 0.7*glutGet(GLUT_SCREEN_HEIGHT); if (ScreenY==0) ScreenY = 600;
    glutInitWindowSize( ScreenX/2, ScreenY/2 );
    glutInitWindowPosition( 0.15*glutGet(GLUT_SCREEN_WIDTH), 0.15*glutGet(GLUT_SCREEN_HEIGHT) );
    glutCreateWindow( "COMMIT debugger" );
    glutReshapeWindow( ScreenX, ScreenY );

    // Projection and model matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // gluPerspective( 40.0f, (GLfloat)ScreenX / (GLfloat)ScreenY, 10.0f, 1000.0f );
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // gluLookAt(
    //     0.0, 0.0, 2.0*max(pixdim.x*dim.x,pixdim.y*dim.y) * (GLfloat)ScreenY/(GLfloat)ScreenX,
    //     0.0, 0.0, 0.0,
    //     0.0, 1.0, 0.0
    // );

    translation.x	= translation.y = 0;
    zoom			= 0;
    OPENGL_utils::identity( rot );
    OPENGL_utils::identity( id );

    // basic settings
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );

    glEnable( GL_DEPTH_TEST );
    glClearColor( 0.1, 0.1, 0.1, 0.0 );

    // lighting
    glShadeModel( GL_SMOOTH );
    glEnable( GL_NORMALIZE );

    GLfloat white[] = {.5f, .5f, .5f, 1.0f};
    glMaterialfv(GL_FRONT, GL_SPECULAR, white);
    GLfloat shininess[] = {32};
    glMaterialfv(GL_FRONT, GL_SHININESS, shininess);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    GLfloat global_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
    glEnable ( GL_COLOR_MATERIAL );	// use glColor3f() to colorize polygons

    // register CALLBACKS and open window
    glutKeyboardFunc( GLUT__keyboard );
    glutSpecialFunc(  GLUT__specialkey );
    glutDisplayFunc(  GLUT__display );
    glutReshapeFunc(  GLUT__reshape );
    glutMouseFunc(    GLUT__mouse );
    glutMotionFunc(   GLUT__motion );

    GLUT__createMenu();

    glutMainLoop();
}
