#include"draw_ampl.h"
#include"crt.h"
#include"process.h"
#include<strings.h>

static int XX,YY;
static int cHeight,cWidth;


/*-------------------- moi!!!--------------------------*/
typedef void (*Pre_test)(vampl*,vertlink*,int,int*);
typedef int (*see_edge)(vampl*,edgeinvert*,int,int);
typedef int (*Pos_test)(vampl*,vertlink*,int,int,int);
typedef struct { int x,y;}  position;


/* S T A T I C S */

static int         quant;
static position    vertexPos[maxvert];
static position    edgePos[MAXINOUT];

static position    myedgePos[MAXINOUT];
static int         out_place;

static int         direction[maxvert];
static int         dir[maxvert][MAXVALENCE];
static int         dir_position[maxvert];
static int         order[maxvert][MAXVALENCE];

static int         outlegs[maxvert];
static vertlink  * line_array[maxvert];
static int         fromlegs[maxvert];
static int         linesize;
static vertlink    headline;

static int         LASTDIR;

#define FROM_UP (int)1
#define FROM_DOWN (int)0
#define FROM_LEFT (int)2

/* - - - - - - - -*/


int amplitudeFrameX(int nin, int nout)
{ if(nout==1) nout=2;
 cWidth=tg_textwidth("H");
 cHeight=tg_textheight("H");

XX=15*cWidth;
    quant=7/2*cWidth;
    quant=4*cWidth;
XX=(3+nout-1)*quant;
return XX;

}

int amplitudeFrameY(int nin, int nout)
{ if(nout==1) nout=2;
  cWidth=tg_textwidth("H");
  cHeight=tg_textheight("H");
  YY=8*cHeight;
    quant=7/2*cWidth;
    quant=4*cWidth;
  YY=quant*(nout-1)+2*cHeight;
  return YY;
}

static void super_outtextxy(int x, int y, char *name)
{ 
   if (name[0]=='~')
   { int cHeight=tg_textheight("H");
     tg_outtextxy(x,y,name+1);
     if(strlen(name+1)==2) tg_outtextxy(x,y-cHeight/2 ,"~~");
           else            tg_outtextxy(x,y-cHeight/2 ,"~");
   }else tg_outtextxy(x,y,name);
}


static void drawPropagator(int x1,int y1,int x2, int y2, int particle, int virtual)
{
  int spin2=prtclbase[particle-1].spin;
  int anti=prtclbase[particle-1].anti;
  int ln;
  int arrow=0;
 
  switch(spin2)
  { 
   case 0: ln = DottedLn;   break;
   case 1: ln = SolidLn;    break;
   case 2: ln = DashedLn;   break;
   case 3: ln = SolidLn;    break;
   case 4: ln = DashedLn;   break;                                                                                                                                                                        
  }
  tg_setlinestyle(ln, NormWidth); 
  
  
  if(particle < anti)  arrow=1;
  else if(particle > anti) { arrow=-1; particle=anti;}     
  

  switch(arrow)
  {
    case  0:  tg_line(x1,y1,x2,y2); break;
    case  1:  tg_arrowline(x1,y1,x2,y2);break;
    case -1:  tg_arrowline(x2,y2,x1,y1);break;
 } 

 if(virtual)
 { int xj=LeftText;
   int yj=BottomText;
   int xShift=1;
   int yShift=1;
    
   if(x1==x2) { yj=CenterText; if(arrow) xShift= cWidth/3+1;} else 
   if(y1==y2) { xj=CenterText; if(arrow) yShift= cHeight/3+1;}  else
              { if(arrow) yShift= cHeight/4+1;xShift= cWidth/4+1;}
   tg_settextjustify(xj,yj);
   super_outtextxy((x1+x2)/2+xShift,(y1+y2)/2-yShift,prtclbase[particle-1].name);  
 }

      
}

static void in_line_pre(vampl*v,vertlink*head,int fromlegno,int*dum){
   if(++linesize>1) 
   fromlegs[linesize-2]=fromlegno;
   line_array[linesize-1]=head;
}

static int see_taranov(vampl * v, vertlink *head, int fromlegno,
       Pre_test pre_test, see_edge leg_test, Pos_test pos_test){
       int i,j,Iout;
       int vertno,inleg;
       int Ivalence;   
       edgeinvert *leg;
       int legs_order[MAXVALENCE];
       
       vertno= head->vno;
       inleg= head->edno;
       Ivalence=v->valence[vertno];   
       j=0; 
       for(i=0;i<Ivalence;i++){
	  if(i==inleg) continue;
	  legs_order[j++]=i;
       }
       pre_test(v,head,fromlegno,legs_order);
       for (i=0; i<Ivalence ;i++){
	   if(i==inleg)continue;
           leg= &v->vertlist[vertno][i];
           if( (leg->prop & IN_PRTCL)   || 
               (leg->prop & OUT_PRTCL)  )
               {
                if(1==(Iout=leg_test(v,leg,vertno,i))) break;
               }
           else if (!(leg->prop & IN_PRTCL)   &&
                    !(leg->prop & OUT_PRTCL) )
               {
                if (1==(Iout=see_taranov(v,&leg->link,i,
					 pre_test,
					 leg_test,
					 pos_test))) break;
               }
       }
       pos_test(v,head,Iout,i,leg->prop);
       return Iout;
}

static void out_count_pre(vampl*v,vertlink*head,int fromlegno,int*dum){
   int vertno;
   vertno=head->vno;
   outlegs[vertno]=0;
/*   return 0; */
}

static int out_count_leg(vampl*v,edgeinvert*ed,int vertno,int legno){
      outlegs[vertno]++;
      return 0;
}

static edgeinvert *Leg(vampl*v,vertlink*head){
      return &v->vertlist[head->vno][head->edno];
}

static int out_count_pos(vampl*v,vertlink*head,int Iout,int fleg,int prop){
      int previous_vertex;
      int vertno;
      edgeinvert *ed;
      int legs[MAXVALENCE];
      int i,j,k;
      int valen;
      int inleg;
      ed= Leg(v,head);
      vertno=head->vno;
      previous_vertex=ed->link.vno;
      outlegs[previous_vertex]+=outlegs[vertno];
      valen=v->valence[vertno];
      inleg=head->edno;
      j=0;
      for(i=0;i<valen;i++){ 
	 if(i==inleg) continue;
	 if(v->vertlist[vertno][i].prop&OUT_PRTCL) legs[j++]=1;
         else legs[j++]=outlegs[v->vertlist[vertno][i].link.vno];
      }
      for(i=0;i<valen-1;i++) order[vertno][i]=i;
      for(i=0;i<valen-2;i++)
	for(j=0;j<valen-i-2;j++)
	 if(legs[order[vertno][j]]<legs[order[vertno][j+1]]){
	    k=order[vertno][j];
	    order[vertno][j]=order[vertno][j+1];
	    order[vertno][j+1]=k;
	    }
      return 0;
}



static void count_outlegs(vampl*v){
   int i,j,vertno;
   for(i=0;i<linesize;i++){
    vertno=line_array[i]->vno;
    for(j=0;j<v->valence[vertno];j++)
       { 
	  if(j==fromlegs[i] || j==line_array[i]->edno) continue;
	  if((v->vertlist[vertno][j].prop&OUT_PRTCL)==0)
	     see_taranov(v,&v->vertlist[vertno][j].link,0,
			 out_count_pre,out_count_leg,out_count_pos);
	     else
	     outlegs[vertno]++;
       }
  }
}


static int in_line_leg(vampl*v,edgeinvert*ed,int vertno,int legno){
   if(ed->prop&IN_PRTCL)  return 1;
   if(ed->prop&OUT_PRTCL) return 2;
   return 0;
}

static int in_line_pos(vampl*v,vertlink*head,int Iout,int fleg,int prop){
   if (prop&IN_PRTCL) fromlegs[linesize-1]=fleg;
   if(Iout!=1){linesize--;}
   return 0;
}


static int in_line(vampl*v,vertlink*head_line){
   int i;
   edgeinvert *ed;
   for(i=0;i<v->outno;i++){
     ed= Leg(v,&v->outer[i]);
     if(ed->prop&IN_PRTCL){ head_line=&v->outer[i]; break; }
   }
   if(nin == 1) {
      linesize=1;
      fromlegs[linesize-1]=0;
      line_array[linesize-1]=head_line;
      return linesize;
      }
   linesize = 0;
   see_taranov(v,head_line,0,in_line_pre,in_line_leg,in_line_pos);
   return linesize;
}

static void coord_edges(vampl*v,position*edgePos){
    int i; 
    for(i=0;i<nin;i++)
    { edgePos[i].x=1+2*cWidth;
      edgePos[i].y=quant/2+i*quant;
    }
    for(i=nin;i<v->outno;i++)
    { edgePos[i].x=XX-1-2*cWidth;
      edgePos[i].y=cHeight/2+(i-nin)*quant;
    }
}

static void decay_pre(vampl*v,vertlink*head,int dumm,int*legs_order){
      int i,tmp; 
      int vertno;
      int previous_vertex;
      edgeinvert*ed;

      ed= Leg(v,head);
      vertno=head->vno;
      previous_vertex=ed->link.vno;

      direction[vertno]=dir[previous_vertex][dir_position[previous_vertex]];

      if (direction[vertno] == FROM_UP) {

	 for(i=0;i<(v->valence[vertno]-1)/2;i++){
	    tmp=order[vertno][i];
	    order[vertno][i] = order[vertno][v->valence[vertno]-1-i];
	    order[vertno][v->valence[vertno]-1-i]=tmp;
	 }
      }
      for(i=1;i<v->valence[vertno]-1;i++) 
	 legs_order[i]=i;/*order[vertno][i];*/

      dir[vertno][0]=FROM_DOWN;
      if((v->vertlist[vertno][legs_order[0]].prop&(IN_PRTCL|OUT_PRTCL))==0 ) dir[vertno][0]=FROM_LEFT;
      for(i=1;i<v->valence[vertno]-2;i++) dir[vertno][i]=FROM_LEFT;
      dir[vertno][v->valence[vertno]-2]=FROM_UP;

      dir_position[vertno] = 0;
      vertexPos[vertno].x=XX-1-2*cWidth-quant;

      dir_position[previous_vertex]++;
/*      return 0; */
}

static void put_edge(vampl*v,int vertno, int legno, int direct){
   int i;
   LASTDIR=direct;
   for(i=0;i<MAXINOUT;i++){
      if(v->outer[i].vno == vertno && v->outer[i].edno==legno)
       {
          out_place++;
          /*------------ stawim out-particle!!! -----------*/
	  edgePos[i]=myedgePos[nin+out_place-1];
	  if(direct==FROM_DOWN) vertexPos[vertno].y=edgePos[i].y+quant/2;
	  if(direct==FROM_LEFT) vertexPos[vertno].y=edgePos[i].y;
	  break;
       }
   }
}


static int decay_leg(vampl*v,edgeinvert*ed,int vertno,int legno){

   put_edge(v, vertno, legno, dir[vertno][dir_position[vertno]]);
   dir_position[vertno]++;
   vertexPos[vertno].x=MIN(vertexPos[vertno].x,XX-1-2*cWidth-quant);
   return 0;
}

static int decay_pos(vampl*v,vertlink*head,int dum1,int dum2,int dum3)
{
   int previous_vertex;
   int vertno;

   edgeinvert*ed;

   ed= Leg(v,head);
   vertno=head->vno;

   previous_vertex=ed->link.vno;

   if(direction[vertno]==FROM_LEFT) 
   {
      vertexPos[previous_vertex].x=
      MIN(vertexPos[vertno].x-quant,vertexPos[previous_vertex].x);
      vertexPos[previous_vertex].y=vertexPos[vertno].y;
   }
   else
   {
      vertexPos[previous_vertex].x=
          MIN(vertexPos[vertno].x,vertexPos[previous_vertex].x);
      if(direction[vertno]==FROM_DOWN) 
	     vertexPos[previous_vertex].y=myedgePos[nin+out_place-1].y+quant/2;
   }
   return 0;
}


static int Decay(vampl*v,vertlink*head,int dire, position*vertexPos)
{
   int previous_vertex;
   edgeinvert*ed = Leg(v,head);
  
   previous_vertex=ed->link.vno;
   dir_position[previous_vertex]=0;
   dir[previous_vertex][0]=dire;
   see_taranov(v,head,0,decay_pre,decay_leg,decay_pos);

   return 0;
}


static void getPositions(vampl*v,position*vertexPos,position*edgePos)
{
  int line_x;
  int i; 
  int ii,iva,DIR;
  int j,vertno;
  int xshift;
  
  for (i=0;i<maxvert;i++) line_array[i]=(vertlink*)malloc(sizeof(vertlink));
  in_line(v,&headline);
  for(i=0;i<linesize;i++){
    vertno=line_array[i]->vno;
    vertexPos[vertno].x=XX;
  }
  count_outlegs(v);
  coord_edges(v,myedgePos);
  out_place =0;
  if (nin==1) {
	vertno=line_array[0]->vno;
	vertexPos[vertno].x=myedgePos[nin].x-quant;
	Decay(v,line_array[0],FROM_DOWN,vertexPos);
  }
  else
  for(i=0;i<linesize;i++){
    vertno=line_array[i]->vno;
    vertexPos[vertno].x=myedgePos[nin].x-quant;
    ii=0;
    iva=v->valence[vertno]-2;
    for(j=0;j<v->valence[vertno];j++) { 
       if(j==fromlegs[i] || j==line_array[i]->edno) continue;
       DIR=FROM_LEFT;
       if(iva==1 && i==linesize-1 && linesize > 1)DIR=FROM_UP;
       if(iva==2)  DIR=FROM_DOWN;
       if(++ii==2) DIR=FROM_UP;
       DIR=FROM_LEFT;                        
       if((v->vertlist[vertno][j].prop&OUT_PRTCL)==0){
	if(DIR==FROM_UP && ii==1){
	   if(LASTDIR==FROM_UP)
	    vertexPos[vertno].y=myedgePos[nin+out_place-1].y+quant/2;
	   if(LASTDIR==FROM_LEFT)
	    vertexPos[vertno].y=myedgePos[nin+out_place-1].y+2*quant/2;
	   if(LASTDIR==FROM_DOWN)
	    vertexPos[vertno].y=myedgePos[nin+out_place-1].y+3*quant/2;
	} else if(DIR==FROM_UP && ii==2){
	;
	}
	  Decay(v,&v->vertlist[vertno][j].link,DIR,vertexPos);
       }
       else{
	  put_edge(v,vertno,j,DIR);
       }
       }
    }
  line_x=XX;
  for(i=0;i<linesize;i++){
    vertno=line_array[i]->vno;
    line_x=MIN(line_x,vertexPos[vertno].x);
    }
  for(i=0;i<linesize;i++){
    vertno=line_array[i]->vno;
    vertexPos[vertno].x=line_x;
    }

  if(linesize > 1){
    vertno=line_array[0]->vno;
    edgePos[0].y=vertexPos[vertno].y;
                       edgePos[0].x=vertexPos[vertno].x-quant;
    vertno=line_array[linesize-1]->vno;
    edgePos[nin-1].y=vertexPos[vertno].y;
                       edgePos[nin-1].x=vertexPos[vertno].x-quant;
  }
  else {
    vertno=line_array[0]->vno;
    if(nin==1){
    edgePos[0].y=vertexPos[vertno].y;
    edgePos[0].x=vertexPos[vertno].x-quant;
    }
    else if(nin==2) {
    edgePos[0].y=vertexPos[vertno].y-quant/2;
    edgePos[0].x=vertexPos[vertno].x-quant;
      edgePos[1].y=vertexPos[vertno].y+quant/2;
      edgePos[1].x=vertexPos[vertno].x-quant;
      }
  }
  xshift=(XX-edgePos[nin+nout-1].x-edgePos[0].x)/2;
  for(i=0;i<(nin+nout);i++)
        edgePos[i].x=edgePos[i].x+xshift;
  for(i=0;i<v->size;i++) vertexPos[i].x=vertexPos[i].x+xshift;

  if(nout==1)
  { int yshift=quant;
    for(i=0;i<v->size;i++) vertexPos[i].y=vertexPos[i].y+quant/2;
    for(i=0;i<(nin+nout);i++)
           edgePos[i].y=edgePos[i].y+quant/2;
  }         

}


void drawAmpitudeDiagram(vampl* diagr,int mode, int x, int y)
{
  int i,j;

  cHeight=tg_textheight("H");
  cWidth=tg_textwidth("H");
       
  getPositions(diagr, vertexPos,edgePos);
  if(mode&1)
  for(i=0;i<diagr->size;i++) 
  { vertexPos[i].x *=-1;
    for(j=0;j<diagr->valence[i];j++)
    { edgeinvert *v=&diagr->vertlist[i][j];
      if(v->link.vno==nullvert) edgePos[v->link.edno].x *=-1;  
    }
  }    
  for(i=0;i<diagr->size;i++) 
  { position in=vertexPos[i];
    for(j=0;j<diagr->valence[i];j++)
    { edgeinvert *v=&diagr->vertlist[i][j];
      position out;
      if(v->link.vno==nullvert) 
      { out=edgePos[v->link.edno];
        drawPropagator(x+in.x,y+in.y,x+out.x,y+out.y,v->partcl,0);
        if(v->prop&IN_PRTCL) 
        { if(mode&1) tg_settextjustify(LeftText,CenterText); else tg_settextjustify(RightText,CenterText);
          tg_outtextxy(x+out.x,y+out.y,prtclbase[prtclbase[v->partcl-1].anti -1].name); 
        }else 
        { 
          if(mode&1)tg_settextjustify(RightText,CenterText);else  tg_settextjustify(LeftText,CenterText); 
          super_outtextxy(x+out.x,y+out.y,prtclbase[v->partcl-1].name);
/*
          sprintf(num,"%d",v->link.edno+1);
          tg_outtextxy(x+out.x+2*cWidth,y+out.y+cHeight/3,num);
*/          
        }
      }   
      else if(i<v->link.vno)
      {  out=vertexPos[v->link.vno];
         drawPropagator(x+in.x,y+in.y,x+out.x,y+out.y,v->partcl,1);  
      }
    }
  }
}
