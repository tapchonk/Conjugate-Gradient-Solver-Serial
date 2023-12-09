/*
* ddot.c - Compute the dot product of two vectors in parallel using pthreads and SIMD vector intrinsics, most notably AVX2 extensions and Intel's fma3 extension (fuse multiply add)
* instructions allows us to combine the operation of multiplying and adding to vectors together into a single instruction. This saves us a significant amount of cycles spent on 
* multiplication and addition. Furthermore, as we essentially combine two functions into one instruction, we can also make our final outcome more accurate due to the fact
* that with fuse multiply add, we now only have one rounding step as opposed to two.
*/


//pthreads > openmp (depends on the application though)
#include "ddot.h"
#include <immintrin.h>
#include <pthread.h>

//Mutex lock free mechanism

//Change the max thread here
#define MAX_THREAD 4
#define LOOP_FACTOR_DDOT 32

//Array to store the sum of each thread partition
float sumDDOTArray[MAX_THREAD];
int MAXDDOT;

//Define the struct for data to pass to the threads
struct argStruct {
  float *arrX;
  float *arrY;
  int part;
};


void* sumXYArray(void* arg1) {
  //Divide work between threads whilst also taking into account loop factor
  register struct argStruct *args = arg1;
  register int part = args->part;
  //Bitshifts are generally faster than division wherever possible
  register int  start = part*(MAXDDOT/MAX_THREAD), end = (((MAXDDOT/MAX_THREAD))>>5)<<5;
  register int i;
  register float *x = (float *)(args->arrX) + start;
  register float *y = (float *)(args->arrY) + start;
  register __m256 xVecSquareSum1 = _mm256_setzero_ps();
  register __m128 horSum;
  register float tempSum = 0.0;
  for (i = 0; i < end; i = i + LOOP_FACTOR_DDOT) {
    //Inline assembly to use the fmadd instruction
    asm (
      //Load 32 values from the x array into vector registers marked by ymm(0, 2, 4, 6)
      "leaq (%1), %%r14\n\t"
      "vmovups (%%r14), %%ymm0\n\t"
      "vmovups 32(%%r14), %%ymm2\n\t"
      "vmovups 64(%%r14), %%ymm4\n\t"
      "vmovups 96(%%r14), %%ymm6\n\t"
      //Load 32 values from the y array into vector registers marked by ymm(1, 3, 5, 7)
      "leaq (%2), %%r15\n\t"
      "vmovups (%%r15), %%ymm1\n\t"
      /*
       * Use fused multiply add instruction vfmadd231ps (the 231 is meant to denote (input 2 * input 3) + input 1, but it's backwards here)
       * The reason as to why it's backwards relates to the fact that all Intel x86 processors are little endian (should be at the top of the results of lscpu),
       * so the order of the bytes in the vector instructions is reversed and thus we are working with registers backwards
       */
      "vfmadd231ps %%ymm1, %%ymm0, %0\n\t"
      "vmovups 32(%%r15), %%ymm3\n\t"
      "vfmadd231ps %%ymm3, %%ymm2, %0\n\t"
      "vmovups 64(%%r15), %%ymm5\n\t"
      "vfmadd231ps %%ymm5, %%ymm4, %0\n\t"
      "vmovups 96(%%r15), %%ymm7\n\t"
      "vfmadd231ps %%ymm7, %%ymm6, %0\n\t"
      : "+x" (xVecSquareSum1)
      : "r" (x+i), "r" (y+i)
      : "ymm0", "ymm1", "ymm2", "ymm3", "ymm4", "ymm5", "ymm6", "ymm7", "r14", "r15" //There is no particular reason as to why I chose registers 14 and 15, I just like them :)
    ); //End of inline assembly, clobber vector registers
  }
  /*
  * This is an example of horizontal vector sum, we make use of this a lot as it reduces the number of additions from 8 to log2(8) = 3
  * It is usually significantly faster than hadd or just unloading the vector into memory and adding sequentially
  * In first vector sum we split the 256 bit vector into 128 bit vectors and add them together
  * One of the 128 bit vectors is the lower half of the vector, the other is the higher half
  * This is a diagram to help visualize this:
  * Original 256 bit vector:  [[7], [6], [5], [4], [3], [2], [1], [0]]
  * Higher 128 bit vector:    [[7], [6], [5], [4]]
  * Lower 128 bit vector:     [[3], [2], [1], [0]]
  * We then add the two vectors together to get the sum of the vector
  * horSum:                 [[7+3], [6+2], [5+1], [4+0]] 
  */
  horSum      = _mm_add_ps(_mm256_castps256_ps128(xVecSquareSum1), _mm256_extractf128_ps(xVecSquareSum1, 1));
  /*
  * In horSum we shift the higher elements down to the lower elements using (_mm_movehl_ps (move higher lower)) and sum them like so:
  * DC = "don't care"
  * Original 128 bit vector:  [[DC], [DC], [5+1]    , [4+0]    ]
  * Shifted 128 bit vector:   [[DC], [DC], [7+3]    , [6+2]    ]
  * horSum:                   [[DC], [DC], [5+1+7+3], [4+0+6+2]]
  */
  horSum      = _mm_add_ps(horSum, _mm_movehl_ps(horSum, horSum));
  /*
  * Original 128 bit vector:  [[DC], [DC], [DC], [4+0+6+2]]
  * Shifted 128 bit vect:     [[DC], [DC], [DC], [5+1+7+3]]
  * horSum:                   [[DC], [DC], [DC], [4+0+6+2+5+1+7+3]] be careful as the order of addition in not consistent with adding sequentially
  */
  horSum      = _mm_add_ps(horSum, _mm_shuffle_ps(horSum, horSum, 0x1));
  /*Sum is now in lowest end of the vector and we can extract it with _mm_cvtss_f32*/
  tempSum = tempSum + _mm_cvtss_f32(horSum);

  for (; i<(MAXDDOT/MAX_THREAD); i++) {
    tempSum = x[i]*y[i] + tempSum;
  }
  //insert local sum to the array of sums for each thread.
  sumDDOTArray[part] = tempSum;

  pthread_exit(NULL);
}

/**
 * @brief Compute the dot product of two vectors, Now with parallelisation and SIMD vector intrinsics!
 * 
 * @param n Number of vector elements
 * @param x Input vector
 * @param y Input vector
 * @param result Pointer to scalar result value
 * @return int 0 if no error
 */
__inline__ int ddot (const int n, const float * const x, const float * const y, float * const result) {
  //Declare struct for passing arguments to thread
  register int i;
  struct argStruct args[MAX_THREAD];
  for (i = 0; i < MAX_THREAD; i++) {
    args[i].arrX = (float *)x;
    args[i].arrY = (float *)y;
    args[i].part = i;
  }
  MAXDDOT = n;
  //Initialize threads
  pthread_t threads[MAX_THREAD];
  register float local_result = 0.0f;
  //Create threads
  for (i = 0; i < MAX_THREAD; i++) {
    pthread_create(&threads[i], NULL, sumXYArray, (void*)&args[i]);
  }
  //Join threads
  for (i = 0; i < MAX_THREAD; i++) {
    pthread_join(threads[i], NULL);
  }
  for (i = 0; i < MAX_THREAD; i++) {
    local_result = local_result + sumDDOTArray[i];
  }

  //Catches the rest of the work non-divisible by MAX_THREAD
  for (i = ((MAXDDOT/MAX_THREAD)*MAX_THREAD); i<MAXDDOT; i++) {
    local_result = local_result + x[i]*y[i];
  }

  *result = local_result;
  return 0;
}