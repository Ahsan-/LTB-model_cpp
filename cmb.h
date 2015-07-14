#ifndef __PYX_HAVE__emcee_wrapper
#define __PYX_HAVE__emcee_wrapper


#ifndef __PYX_HAVE_API__emcee_wrapper

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(void) set_forth_EmceeInstance(int, int, int, int, MCMC *, MCMC_params *);

#endif /* !__PYX_HAVE_API__emcee_wrapper */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initemcee_wrapper(void);
#else
PyMODINIT_FUNC PyInit_emcee_wrapper(void);
#endif

#endif /* !__PYX_HAVE__emcee_wrapper */
