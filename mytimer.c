#include <stdlib.h>
#include <omp.h>

/**
 * @brief Gets the time
 * 
 * @return double Returns the time
 */
double mytimer(void)
{
   return omp_get_wtime();
}