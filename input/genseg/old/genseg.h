// those have to be consistent with settings.h from search
#define C_OMEGA_R 7.2921151467064e-5
#define C_SIDDAY (2.*M_PI/C_OMEGA_R)
#define EPSILON 0.40909280422232891
#define C_SI 299792458

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#ifdef DEBUG
#define DEBUG_PRINT(x) printf x
#else
#define DEBUG_PRINT(x)
#endif

/* Default configuration file. */
/* Cann be overridden by the first command-line arg */
#define INI_FNAME "genseg.ini"

/* Buffer size for the outlier removal algorithm */
#define BUFSIZE 1<<15
#define MAX_LINE 8192

int fGrubbsOutliersMany (float *, float *, int, int, double, const char *);
