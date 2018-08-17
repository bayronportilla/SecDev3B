#include <stdio.h>

int counterLines(char *infile){

 int NDAT,c;
 static FILE *pf;

 if((pf=fopen(infile,"r"))==NULL)
   printf("no puedo abrir archivo %s\n",infile);

 NDAT=0;

 while((c=fgetc(pf))!= EOF)
   {
     if(c == '\n')
       {
 ++NDAT;
       }
   }

 fclose(pf);

 return NDAT;
}
