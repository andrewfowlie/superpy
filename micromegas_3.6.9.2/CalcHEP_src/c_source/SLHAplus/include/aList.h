#ifndef __ALIST__ 
#define __ALIST__

#define aList3(type) type x1,type x2,type x3
#define aList4(type) aList3(type),type x4

#define aList6(type) aList4(type),type x5,type x6
#define aList9(type) aList6(type),type x7,type x8,type x9

#define aList10(type) aList9(type),type x10 
#define aList15(type) aList10(type),type x11,type x12,type x13,type x14,type x15

#define aList16(type) aList15(type),type x16 
#define aList25(type) aList16(type),type x17,type x18,type x19,type x20,type x21,type x22,type x23,type x24,type x25

#endif
