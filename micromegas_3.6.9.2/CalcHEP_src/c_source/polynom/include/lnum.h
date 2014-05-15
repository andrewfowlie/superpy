#ifndef __LNUM__
#define __LNUM__


#if defined (NUM_DOUBLE)      /* double */

#define NUM_TYPE double
#include<math.h>
#define NUM_ONE  1.
#define NUM_ZERO 0.  
#define STR_NUM "lf"
#define NUM_STR ".0f"

#define DIV(a,b) ((a)>0. ? ((b)>0. ?  floor((a)/(b))   : -floor(-(a)/(b))) :\
			   ((b)>0. ?  -floor(-(a)/(b)) : floor((a)/(b))   ) )
#define REST(a,b)  ((a)-DIV(a,b)*(b))


#elif defined(NUM_LONG_LONG)  /* long long */
#define NUM_TYPE long long  
#define NUM_ONE  1
#define NUM_ZERO 0  

#if ! defined (NUM_STR)
#define   NUM_STR "lld"
#endif
#define STR_NUM  NUM_STR


#define DIV(a,b)   ((a)/(b)) 
#define REST(a,b)  ((a)%(b))

#else                         /* long */   

#define NUM_TYPE long
#define NUM_ONE  1
#define NUM_ZERO 0  
#define STR_NUM  "ld"
#define NUM_STR  "ld"

#define DIV(a,b)   ((a)/(b)) 
#define REST(a,b)  ((a)%(b))

#endif 

#endif
