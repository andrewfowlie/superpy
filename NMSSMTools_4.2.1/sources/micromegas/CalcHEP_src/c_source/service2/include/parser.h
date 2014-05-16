#ifndef __PARSER_
#define __PARSER_

typedef void* ( * operation )(char * , int, void **);
typedef void* ( * rdelement)(char * s);
typedef void  ( *clear)(void *);

#define braketexpected       1
#define unexpectedcharacter  2
#define operationexpected    3
#define toolongidentifier    4
#define toomanyagruments     5
#define cannotread           6
#define cannotevaluate       7

#define unknownidentifier     8
#define unexpectedoperation   9
#define unknownfunction      10
#define wrongnumberofagr     11
#define typemismatch         12

#define naninoperation       13 

#define toolargenumber       14
#define indexuncompatibility 15
#define rangecheckerror      16


#define errortp  0
#define numbertp 1
#define polytp   2
#define rationtp 3
#define vectortp 4
#define indextp  5
#define tenstp   6
#define spintp   7
#define etenstp  8



extern int  rderrcode, rderrpos;


extern void * readExpression(char*source,rdelement rd,operation act,clear del);

extern char * errmesstxt(int ercodr);

#endif
