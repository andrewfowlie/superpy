
#include"SLHAplus.h"
#include"../../include/nType.h"


/* 
   This file contains diagonalizing routines presented in the format adopted 
   for CalcHEP. Input/output parameters in original Jacobi routines are  
   matrices and vectors. From other side  CalcHEP  can work only  with 
   external  functions whose parameters are single  numbers as return value 
   is also one  number. 
     
   This problem is solved by the following way. We write new interface for
   Jacobi routines  where  input matrix is passes as  a list of matrix
   elements each of them is a  parameter. Instead of output matrices new 
   programs   return identifier (id) of memory location of obtained results.
   We have auxiliary routines which return  eigenvalues  and matrix elements
   using id. 
*/

extern int FError;

static int idMax=-1;
static int * dim=NULL;
static REAL **  rV=NULL;
static REAL **  rU=NULL;
static COMPLEX ** cV=NULL;
static COMPLEX ** cU=NULL;
static REAL ** ev=NULL;
static int idCur=0;

/*idLim  maximum  id. id is increased automatically until  
  initDiagonal  call( normal case) or until the limit is reached.
  The limit is set for the case the user forget to call initDiagonal.
  in the begining of cycle.
 */ 
#define idLIM 200

/* numeration of independent elements of symmetry/hermit matrices */
#define sMT(i,j) ((i)*nDim-((i)*((i)+1))/2 +(j))


/* extendData allocates memory for new results of diagonalizing  */

static void extendData(int id, int nDim, int std)
{ int i;
  if(id>idMax)
  {
    rV=(REAL **)  realloc( rV,(id+1)*sizeof(REAL*));
    rU=(REAL **)  realloc( rU,(id+1)*sizeof(REAL*));
    cV=(COMPLEX **) realloc( cV,(id+1)*sizeof(COMPLEX*));
    cU=(COMPLEX **) realloc( cU,(id+1)*sizeof(COMPLEX*));  
    ev=(REAL **)  realloc(ev, (id+1)*sizeof(REAL*));
    dim=(int*)      realloc(dim,(id+1)*sizeof(int));
    for(i=idMax+1;i<=id;i++) 
    { rV[i]=NULL;rU[i]=NULL;cV[i]=NULL;cU[i]=NULL;ev[i]=NULL; dim[id]=0;}
    idMax=id;
  }
  
  if(std==1 || std==2)
  {  if(cV[id])    
     { 
       free(cV[id]);  cV[id]=NULL;
       if(cU[id]){ free(cU[id]); cU[id]=NULL;}
       free( ev[id]); ev[id]=NULL;
       dim[id]=0;
     }  
     if(nDim!=dim[id])
     {  
        if(rV[id]) { free(rV[id]);  rV[id]=NULL;}
        if(rU[id]) { free(rU[id]);  rU[id]=NULL;}
        if(ev[id]) { free( ev[id]); ev[id]=NULL;} 
        dim[id]=0;
     }
     
     if(dim[id]==0)
     { rV[id]=malloc(sizeof(REAL)*nDim*nDim);
       ev[id]=malloc(sizeof(REAL)*nDim);
       dim[id]=nDim;
     } 
     if(rU[id]==NULL && std==2) rU[id]=malloc(sizeof(REAL)*nDim*nDim);       
  }

  if(std==3 || std==4)
  {  if(rV[id])    
     { 
       free(rV[id]);  rV[id]=NULL;
       if(rU[id]){ free(rU[id]); rU[id]=NULL;}
       free( ev[id]);     ev[id]=NULL;
       dim[id]=0;
     }  
     if(nDim!=dim[id])
     {  
        if(cV[id]) { free(cV[id]);  cV[id]=NULL;}
        if(cU[id]){ free(cU[id]); cU[id]=NULL;}
        if(ev[id])     { free( ev[id]);     ev[id]=NULL;} 
        dim[id]=0;
     }
     
     if(dim[id]==0)
     { cV[id]=malloc(sizeof(COMPLEX)*nDim*nDim);
       ev[id]=malloc(sizeof(REAL)*nDim);
       dim[id]=nDim;
     } 
     if(cU[id]==NULL && std==4) cU[id]=malloc(sizeof(COMPLEX)*nDim*nDim);       
  }
}

/* Initialization of memory allocation for  storage results of 
   diagomalizing in compute memory
*/

int initDiagonal(void){ idCur=0; return 0;}

/* Adopted Jacobi routines. See explanations in arXiv:1008.0181 */

int rDiagonal(int nDim,...) 
{ va_list ap;
  REAL*MassM=malloc(sizeof(REAL)*nDim*nDim);
  int i,j,id;
   
  va_start(ap,nDim);
  for(i=0;i<nDim;i++)for(j=i;j<nDim;j++) MassM[sMT(i,j)]=va_arg(ap, REAL); 
//for(i=0;i<nDim;i++)for(j=i;j<nDim;j++) printf("%d %d %E\n",i,j, (doubble)MassM[sMT(i,j)]);
  va_end(ap);

  if(idCur>idLIM) idCur=0;
  id=idCur++;
  extendData(id,nDim,1);
  FError=FError|rJacobi(MassM, nDim, ev[id], rV[id]); 
  free(MassM);
//for(i=0;i<nDim;i++) printf(" %E ",(double)ev[id][i]);
//printf("\n");  
  
  if(!FError){ return id;}
  else return -1;
}

int cDiagonalH(int nDim,...) 
{ va_list ap;
  COMPLEX*MassM=malloc(sizeof(COMPLEX)*nDim*nDim);
  int i,j,id;
   
  va_start(ap,nDim);
  for(i=0;i<nDim;i++)for(j=i;j<nDim;j++) MassM[sMT(i,j)]=va_arg(ap, COMPLEX); 
  va_end(ap);

  if(idCur>idLIM) idCur=0;
  id=idCur++;
  extendData(id,nDim,3);
  FError=FError|cJacobiH(MassM, nDim, ev[id], cV[id]); 
  free(MassM);
  
  if(!FError){ return id;}
  else return -1;
}

int cDiagonalA(int nDim,...) 
{ va_list ap;
  COMPLEX*MassM=malloc(sizeof(COMPLEX)*nDim*nDim);
  int i,j,id;
   
  va_start(ap,nDim);
  for(i=0;i<nDim;i++)for(j=0;j<nDim;j++) MassM[i*nDim+j]=va_arg(ap, COMPLEX); 
  va_end(ap);

  if(idCur>idLIM) idCur=0;
  id=idCur++;
  extendData(id,nDim,4);
  FError=FError|cJacobiA(MassM, nDim, ev[id],  cU[id], cV[id]); 
  free(MassM);
  if(!FError){ return id;}
  else return -1;
}

int cDiagonalS(int nDim,...) 
{ va_list ap;
  COMPLEX*MassM=malloc(sizeof(COMPLEX)*nDim*nDim);
  int i,j,id;
   
  va_start(ap,nDim);
  for(i=0;i<nDim;i++)for(j=i;j<nDim;j++) MassM[sMT(i,j)]=va_arg(ap, COMPLEX); 
  va_end(ap);

  if(idCur>idLIM) idCur=0;
  id=idCur++;
  extendData(id,nDim,4);
  FError=FError|cJacobiS(MassM, nDim, ev[id], cV[id]); 
  free(MassM);
  
  if(!FError){ return id;}
  else return -1;
}



int rDiagonalA(int nDim,...) 
{ va_list ap;

  REAL* MassM=malloc(sizeof(REAL)*nDim*nDim);
  int i,j,id;
   
  va_start(ap,nDim);
  for(i=0;i<nDim;i++) for(j=0;j<nDim;j++)  MassM[i*nDim+j]=va_arg(ap,REAL);
  va_end(ap);
  
  if(idCur>idLIM) idCur=0;
  id=idCur++;
  extendData(id,nDim,2);

  FError=FError|rJacobiA(MassM, nDim, ev[id], rU[id],rV[id]); 
  free(MassM);  
  return id;
}

/* CalcHEP adopted output of eigenvalues and rotation matrices */

REAL MassArray(int id,  int i)
{ if(id>idMax||i<1 ||i>dim[id]) {FError=1; return 0;}
  return  ev[id][i-1];  
}

REAL MixMatrix(int id, int i,int j)
{ if(id>idMax|| i<1 || i>dim[id] ||j<1 ||j>dim[id]) {FError=1; return 0;}
  return  rV[id][j-1+dim[id]*(i-1)];  
}

REAL MixMatrixU(int id, int i,int j)
{ if(id>idMax || i<1 || i>dim[id] ||j<1 ||j>dim[id]||!rU[id]) {FError=1; return 0;}
  return  rU[id][j-1+dim[id]*(i-1)];  
}

COMPLEX cMixMatrix(int id, int i,int j)
{ if(id>idMax || i<1 || i>dim[id] ||j<1 ||j>dim[id]){FError=1; return 0;}
  return  cV[id][j-1+dim[id]*(i-1)];  
}

COMPLEX cMixMatrixU(int id, int i,int j)
{ if(id>idMax || i<1 || i>dim[id] ||j<1 ||j>dim[id]||!cU[id]) {FError=1; return 0;}
  return  cU[id][j-1+dim[id]*(i-1)];  
}

