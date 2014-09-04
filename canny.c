/****************************
** JUAMMY LORA              **
**  ASSIGNMENT 1 - CANNY    **
**   9-3-14   CAP4453       **
******************************/
#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#include <stdlib.h>
#define  PICSIZE 256
#define  MAXMASK 100

         int    pic[PICSIZE][PICSIZE];
         double outpic1[PICSIZE][PICSIZE];
         double outpic2[PICSIZE][PICSIZE];
         int    edgeflag[PICSIZE][PICSIZE];
         double mask[MAXMASK][MAXMASK];

         double conv[PICSIZE][PICSIZE];

main(argc,argv)
int argc;
char **argv;
{
        int     i,j,p,q,s,t,mr,centx,centy;
        double  maskx, masky, sum,sig,maxival,minival,maxval,ZEROTOL;
        FILE    *fo1, *fo2, *fo3, *fp1, *fopen();
        char    *foobar;
        float   percnt;

        //this first argument reads in the INPUT file
        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");

		//this argument reads in the name of the output file
        argc--; argv++;
        foobar = *argv;
        fo1=fopen(foobar,"wb");
        //Added File header to magnitude output file
        fprintf(fo1, "P5\n 256 256\n 255\n");

        //this argument reads in the sigma value
        argc--; argv++;
        foobar = *argv;
        sig = atof(foobar);

        //this argument reads in the %%PERCENT%% value
        argc--; argv++;
        foobar = *argv;
        percnt = atof(foobar);

        //lets create the PEAKS output and insert the Header
        fo2 = fopen("peaks.pgm", "wb");
        fprintf(fo2, "P5\n 256 256\n 255\n");

        //lets now create the DOUBLE-THRESHOLD output file & insert Header
        fo3 = fopen("double-thrshld.pgm", "wb");
        fprintf(fo3, "P5\n 256 256 \n 255\n");


        mr = (int)(sig * 3);
        centx = (MAXMASK / 2);
        centy = (MAXMASK / 2);

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                }
        }
//HERE WE TAKE THE FIRST DERIVATIVE OF THE GAUSSIAN FOR THE X MASK
        for (p=-mr;p<=mr;p++)
        {  for (q=-mr;q<=mr;q++)
           {
              maskx = (  (2-(((p*p)+(q*q))/(sig*sig))) * (p*(exp((-1/2)*((p/sig)*(p/sig)))))  );
              (mask[p+centy][q+centx]) = maskx;
           }
        }
//HERE WE TAKE THE FIRST DERIVATIVE OF THE GAUSSIAN FOR THE Y MASK
        for (p=-mr;p<=mr;p++)
        {  for (q=-mr;q<=mr;q++)
           {
              masky = (  (2-(((p*p)+(q*q))/(sig*sig))) * (p*(exp((-1/2)*((p/sig)*(p/sig)))))  );
              (mask[p+centy][q+centx]) = masky;
           }
        }

//THIS IS THE CONVOLUTION FOR THE MASK X
        for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sum = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sum += pic[i+p][j+q] * mask[p+centy][q+centx];
                }
             }
             outpic1[i][j] = sum;
             conv[i][j] = sum;
          }
        }

//THIS IS THE CONVOLUTION FOR THE MASK Y
        for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sum = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sum += pic[i+p][j+q] * mask[p+centy][q+centx];
                }
             }
             outpic2[i][j] = sum;//maybe change the outpic
             conv[i][j] = sum;
          }
        }
        maxval  = 0;
        maxival = 0;
        minival = 255;
//HERE IS THE MAGNITUDE FORMULA
//this takes the magnitude

        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             conv[i][j]=sqrt((double)((outpic1[i][j]*outpic1[i][j]) +
                                      (outpic2[i][j]*outpic2[i][j])));
             if (conv[i][j] > maxival)
                maxival = conv[i][j];

           }
        }

        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             if (outpic1[i][j] > maxival)
                maxival = outpic1[i][j];
             if (outpic1[i][j] < minival)
                minival = outpic1[i][j];
           }
        }

        if (fabs(maxival) > fabs(minival))
           maxval = fabs(maxival);
        else
           maxval = fabs(minival);

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
          {
             outpic1[i][j] = ((((outpic1[i][j]) / maxval) + 1) * 127);
             fprintf(fo1,"%c",(char)((int)(outpic1[i][j])));
          }
        }

//        for (i=mr;i<=255-mr;i++)
//        {  for (j=mr;j<=255-mr;j++)
//           {
//                         outpic2[i][j] = 0;
//             if (conv[i][j] > ZEROTOL)
//             {
//               for (p=-1;p<=1;p++)
//               {
//                 for (q=-1;q<=1;q++)
//                 {
//                   if (conv[i+p][j+q] < -ZEROTOL)
//                   {
//                     outpic2[i][j] = 255;
//                   }
//                 }
//               }
//             }
//             else if ((fabs)(conv[i][j]) < ZEROTOL)
//             {
//                     if (((conv[i+1][j] > ZEROTOL) &&
//                          (conv[i-1][j] < -ZEROTOL))   ||
//                         ((conv[i+1][j] < -ZEROTOL) &&
//                          (conv[i-1][j] > ZEROTOL)))
//                     {
//                       outpic2[i][j] = 255;
//                     }
//                     else if (((conv[i][j+1] > ZEROTOL) &&
//                               (conv[i][j-1] < -ZEROTOL))   ||
//                              ((conv[i][j+1] < -ZEROTOL) &&
//                               (conv[i][j-1] > ZEROTOL)))
//                     {
//                       outpic2[i][j] = 255;
//                     }
//                     else if (((conv[i+1][j+1] > ZEROTOL) &&
//                               (conv[i-1][j-1] < -ZEROTOL))   ||
//                              ((conv[i+1][j+1] < -ZEROTOL) &&
//                               (conv[i-1][j-1] > ZEROTOL)))
//                     {
//                       outpic2[i][j] = 255;
//                     }
//                     else if (((conv[i+1][j-1] > ZEROTOL) &&
//                               (conv[i-1][j+1] < -ZEROTOL))   ||
//                              ((conv[i+1][j-1] < -ZEROTOL) &&
//                               (conv[i-1][j+1] > ZEROTOL)))
//                     {
//                       outpic2[i][j] = 255;
//                     }
//             }
//           }
//        }

        //this is for display purposes
        for (i=0;i<256;i++){
            for (j=0;j<256;j++){
                  if (outpic2[i][j] == 255) outpic2[i][j]=0;
             else outpic2[i][j]=255;
             fprintf(fo2,"%c",(char)((int)(outpic2[i][j])));
          }
        }
}
