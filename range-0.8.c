/* range-0.8.c                                       10/2007 jlong@jimlong.org
 *
 * RANGE - RAndom Network GEnerator
 * print to stdout a valid transcription network in the NEMO (NEtwork
 * MOtif) language, which if piped into nemo2sbml will output SBML for
 * input to a biochemical simulator, such as COPASI.
 *
 * compile: gcc -o range-0.8 range-0.8.c -lm
 *
 * usage ./range-0.8 <number of nodes in network> | ./nemo2sbml
 *
 * Copyright (C) 2007, University of Alaska Fairbanks
 * Biotechnology Computing Research Group
 * Author: James Long
 *-------------------------------------------------------------------------------
 * RANGE BSD License
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    * Redistributions of source code must retain the above copyright notice, 
 *      this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright notice,
 *      this list of conditions and the following disclaimer in the documentation
 *      and/or other materials provided with the distribution.
 *    * Neither the name of the University of Alaska Fairbanks nor the names of 
 *      its contributors may be used to endorse or promote products derived from 
 *      this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE 
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ------------------------------------------------------------------------------
 * RANGE GPL License
 * 
 * This project consists of free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This resource is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_GENES     16000 /* max genes in network */
#define MIN_GENES        32
#define VERSION        "0.8"

int *bin, bufSz, didTML=0, done, maxDegree, maxGenesInDOR, *nodeDegree,
    numGenes, nextGene=0, numAtRegulatedGeneDegree, poolSize, 
    regulatedGeneDegree, *tfList, tfListSize;
char *buf, tmBuf[65536];

void AddMotifs(void);
int BuildDOR(void);
void CheckBufSz(int);
double P(int);
void UpdatePowerLawDistributionForDOR(void);
void UpdatePowerLawDistributionForTML(void);

int main(int argc, char *argv[])
{
  int i, j, k, currentMaxDegree, didDOR, didTML, neededBufSz=0, option;
  char c1, c2, buf2[32];
  
  if(argc == 1)
  {
    printf("usage: range [options] <number of genes in network, >= %d, <= %d>\n", MIN_GENES, MAX_GENES);
    printf("              -h --help\n");
    printf("              -v print version\n");
    return 1;
  }
  
  numGenes = atoi(argv[1]);
  
  /* options parsing */
  while((option = getopt(argc, argv, "hv")) > 0)
  {
    switch(option)
    {
      case 'h':
        printf("print to stdout a valid random transcription\n");
        printf("network in the NEMO (NEtwork MOtif) language\n");
        printf("usage: range [options] <number of genes in network, >= %d, <= %d>\n", MIN_GENES, MAX_GENES);
        printf("              -h --help\n");
        printf("              -v print version\n");
        return 0;
     
      case 'v':
        printf("ver %s\n", VERSION);
        return 0;
       
      default: break;
    }
  }
  
  if((argc > 1) && (numGenes <= MAX_GENES) && (numGenes >= MIN_GENES))
  {
    /* Compute the maximum degree of a node. For our power-law distribution,
     * the probability that a node has k edges is P(k) = 0.62/(k^2). The max
     * degree then is that k for which we expect to find just one node with
     * degree k. 
     */
     
    maxDegree = (int)ceil(sqrt(0.62*(double)numGenes));

    /* keep track of node degrees */
    nodeDegree = (int *) malloc(numGenes*sizeof(int)); 
    if(!nodeDegree)
    {
      fprintf(stderr, "range: malloc error, returning...\n");
      return 1;
    }
    
    for(i=0; i<numGenes; i++)
      nodeDegree[i] = 0;

    /* power-law bins */
    bin = (int *) malloc((maxDegree+1)*sizeof(int)); 
    if(!bin)
    {
      fprintf(stderr, "range: malloc error, returning...\n");
      return 1;
    }
    
    for(i=0; i<maxDegree+1; i++)
      bin[i] = 0;
    
    /* list of transcription factors */
    tfList = (int *) malloc(sizeof(int)); /* will reallocate as we go */
    if(!tfList)
    {
      fprintf(stderr, "range: malloc error, returning...\n");
      return 1;
    }
    tfListSize = 1;
    tfList[tfListSize-1] = numGenes - tfListSize;
      
    /* High-level view of the algorithm for construction the network:
     * Iteratively create a network backbone composed of DORs with
     * network motifs attached to every other regulated gene of the DOR.
     */
    
    maxGenesInDOR = poolSize = maxDegree;
    regulatedGeneDegree = 1;
    numAtRegulatedGeneDegree = (int)ceil((double)maxGenesInDOR*P(regulatedGeneDegree));

    printf("[\n");
    strcpy(tmBuf, ",\nTMLIST(\n");
    didDOR = 0;
    while(1)
    {
      if(numAtRegulatedGeneDegree>0)
      {
        if(BuildDOR())
        {
          if(didDOR)
            printf(",\n");
          
          if(buf[strlen(buf)-4]==',') buf[strlen(buf)-4] = ' ';
          
          printf("%s", buf);

          UpdatePowerLawDistributionForDOR();
             
          /* now march through buf and hang a motif at every other gene,
           * being careful not to violate current power-law distribution
           */
          AddMotifs();
          
          didDOR = 1;
        }
        else
          didDOR = 0;
          
        free(buf);
        if(done) break;
      }
      else
        didDOR = 0;
        
      if(nextGene>=tfList[tfListSize-1]) break;
      if(numAtRegulatedGeneDegree<=0)    break;
      if(maxGenesInDOR<2)                break;
    }

    /* print the TMLIST */
    printf("%s", tmBuf);
    if(didTML) printf("\n)\n");

    /* now let the transcription factors regulate each other,
     * and G0 regulate up to maxDegree-1 of them 
     */
    bin[nodeDegree[0]]--; 
    printf(",\nGLIST(\n");

    currentMaxDegree = maxDegree;
    for(i=0; i<tfListSize-1; i++)
    {
      if(drand48() > 0.5) c1 = '+';
      else                c1 = '-';
      if(drand48() > 0.5) c2 = '+';
      else                c2 = '-';
      
      if((nodeDegree[0] < maxDegree-1) && i%2)
      {
        printf("  G%d(P%d%c,P0%c", tfList[i], tfList[i+1], c1, c2);
        nodeDegree[0]++;
      }
      else
        printf("  G%d(P%d%c", tfList[i], tfList[i+1], c1);
      
      if(bin[nodeDegree[tfList[i+1]]] > 0)
         bin[nodeDegree[tfList[i+1]]]--;
         
          nodeDegree[tfList[i+1]]++;
      bin[nodeDegree[tfList[i+1]]]++;
          nodeDegree[tfList[i]]  = 1;
      bin[nodeDegree[tfList[i]]] = 1;

      j = 1;
      while((nodeDegree[tfList[i]] < currentMaxDegree) && bin[nodeDegree[tfList[i]]]
           < (int)ceil((double)numGenes*P(nodeDegree[tfList[i]])))
      {
        if(drand48() > 0.5) c1 = '+';
        else                c1 = '-';
        
        printf(",P%d%c", tfList[tfListSize-i-1-j], c1);
        if(bin[nodeDegree[tfList[tfListSize-i-1-j]]] > 0)
           bin[nodeDegree[tfList[tfListSize-i-1-j]]]--;
           
            nodeDegree[tfList[tfListSize-i-1-j]]++;
        bin[nodeDegree[tfList[tfListSize-i-1-j]]]++;
        bin[nodeDegree[tfList[i]]]--;
            nodeDegree[tfList[i]]++;
        bin[nodeDegree[tfList[i]]]++;

        j++;
      }
      printf("),\n"); 
      currentMaxDegree--;
    }
    bin[nodeDegree[0]]++;
    
    /* finish up any genes that got missed */
    bufSz = 128;
    buf = (char *) malloc(bufSz);
    if(!buf)
    {
      fprintf(stderr, "range: malloc error, exiting...\n");
      exit(1);
    }
    
    didTML = 0;
    strcpy(buf, ",\nTMLIST(\n");
    
    while(nextGene < tfList[i])
    {
      if(drand48() > 0.5) c1 = '+';
      else                c1 = '-';
      
      if(nextGene == tfList[i]-1)
      {
        printf("  G%d(P%d%c),\n", nextGene, nextGene, c1);
        if(bin[nodeDegree[nextGene]] > 0) bin[nodeDegree[nextGene]]--;
        
            nodeDegree[nextGene] += 2;
        bin[nodeDegree[nextGene]]++;
        
      }
      else
      {
        printf("  G%d(P%d%c),\n", nextGene, nextGene+1, c1);
        nodeDegree[nextGene] = 1;
        bin[1]++;
        
        if(bin[nodeDegree[nextGene+1]] > 0) bin[nodeDegree[nextGene+1]]--;
        
            nodeDegree[nextGene+1]++;
        bin[nodeDegree[nextGene+1]]++;
            nodeDegree[nextGene] += 1;
        bin[nodeDegree[nextGene]]++;
      }
     
      /* add a sim? */
      j = 0; 
      /* find an unfilled bin */
      while((bin[maxDegree-j] >= 
            (int)floor((double)numGenes*P(maxDegree-j))) && 
            (j < maxDegree))
        j++;
      
      if(((nextGene+1) < (tfList[i]-maxDegree-j-2)) && (maxDegree-j-2 > 0))
      {
        neededBufSz += 32;
        CheckBufSz(neededBufSz);
        sprintf(buf2, "  P%d(%cG%d", nextGene, c1, nextGene+1);
        strcat (buf, buf2);
        nextGene += 2;
        didTML = 1;
          
        for(k=0; k<maxDegree-j-2; k++)
        {
          neededBufSz += 8;
          CheckBufSz(neededBufSz);
          sprintf(buf2, ",G%d", nextGene++);
          strcat (buf, buf2);
          bin[1]++;
        }
        
        bin[maxDegree-j] += 1;
          
        neededBufSz += 16;
        CheckBufSz(neededBufSz);

        while((bin[maxDegree-j] >= 
            (int)floor((double)numGenes*P(maxDegree-j))) && 
            (j < maxDegree))
          j++;
        if(((nextGene+1) < (tfList[i]-maxDegree-j-2)) && (maxDegree-j-2 > 0))
          strcat(buf, "),\n");
        else
          strcat(buf, ")\n)\n");
      }
      else
      {
        bin[2]++;
        nextGene++;
      }
    }

    printf("  G%d(P%d%c)\n)\n", nextGene, tfList[0], c1);
    if(didTML) 
    {
      UpdatePowerLawDistributionForTML();
      printf("%s", buf);
    }
    
    free(buf);
    printf("]\n");
    
    return 0;
  }
  else
  {
    printf("usage: range [options] <number of genes in network, >= %d, <= %d>\n", MIN_GENES, MAX_GENES);
    printf("              -h --help\n");
    printf("              -v print version\n");
    return 1;
  }
}

/* add a motif only if it does not violate the power-law distribution */
void AddMotifs(void)
{
  int i, j, k, gene=0, geneNum, must_break;
  char buffer[128], c1, c2, c3, *p, *p2, tempBuf[32];
  double d;
  
  p2 = buf;
  while(p=strstr(p2, "G"))
  {
    gene++;

    if(gene%2)
    {
      /* which gene is this? */
      i = 0;
      buffer[0] = 0x0;
      while(isdigit((p+1)[i]))
        buffer[i] = (p+1)[i++];
      buffer[i] = 0x0;
      
      geneNum = atoi(buffer);

      /* determine the minimum amount that we may increment the node's
       * degree without exceeding the power-law distribution
       */
      for(i=2;i<maxDegree-nodeDegree[geneNum]-(tfListSize/4);i++)
      {
        must_break = 0;
        if(bin[nodeDegree[geneNum]+i] > 
          (int)ceil((double)numGenes*P(nodeDegree[geneNum]+i)))
          continue; /* full bin */
        else
        {
          switch(i)
          {
            case 2:  /* add FFL */
                       
                     if(bin[2] + 2 > (int)ceil((double)numGenes*P(2)))
                       break;
            
                     /* make sure one of the genes is not a TF */
                     for(j=0; j<tfListSize; j++)
                       for(k=0; k<i; k++)
                         if(nextGene+k == tfList[j]) must_break = 1;
                     if(must_break) break;
                     
                     d = drand48();
                     if(     d < 0.125) {c1 = '+'; c2 = '+'; c3 = '+';}
                     else if(d < 0.250) {c1 = '+'; c2 = '+'; c3 = '-';}
                     else if(d < 0.375) {c1 = '+'; c2 = '-'; c3 = '+';}
                     else if(d < 0.500) {c1 = '+'; c2 = '-'; c3 = '-';}
                     else if(d < 0.625) {c1 = '-'; c2 = '+'; c3 = '+';}
                     else if(d < 0.750) {c1 = '-'; c2 = '+'; c3 = '-';}
                     else if(d < 0.875) {c1 = '-'; c2 = '-'; c3 = '+';}
                     else               {c1 = '-'; c2 = '-'; c3 = '-';}
                     
                     if(didTML) strcat(tmBuf, ",\n");
                     sprintf(buffer, "  P%d(%cG%d%cG%d%c)", geneNum, c1, nextGene, c2, nextGene+1, c3);
                     strcat(tmBuf, buffer);
                     
                     /* maintain power law distribution bins */
                     if(bin[nodeDegree[geneNum]] > 0) bin[nodeDegree[geneNum]]--;
                     
                         nodeDegree[geneNum] += 2;
                     bin[nodeDegree[geneNum]]++;
                         nodeDegree[nextGene] += 2;
                     bin[nodeDegree[nextGene]]++;
                         nodeDegree[nextGene+1] += 2;
                     bin[nodeDegree[nextGene+1]]++;
                     nextGene += 2;
                     poolSize += 2;
                     
                     didTML = 1;
                     goto outtaHere;
                     
            default: /* case > 2 */
            
                     /* make sure one of the genes is not a TF */
                     for(j=0; j<tfListSize; j++)
                       for(k=0; k<i; k++)
                         if(nextGene+k == tfList[j]) must_break = 1;
                     if(must_break) break;
            
                     if(1)
                     {
                       /* first try a multi-FFL */

                       if(bin[i] + 1   > (int)ceil((double)numGenes*P(i)))
                         goto sim;

                       if(bin[2] + i-1 > (int)ceil((double)numGenes*P(2)))
                         goto sim;

                       d = drand48();
                       if(     d < 0.125) {c1 = '+'; c2 = '+'; c3 = '+';}
                       else if(d < 0.250) {c1 = '+'; c2 = '+'; c3 = '-';}
                       else if(d < 0.375) {c1 = '+'; c2 = '-'; c3 = '+';}
                       else if(d < 0.500) {c1 = '+'; c2 = '-'; c3 = '-';}
                       else if(d < 0.625) {c1 = '-'; c2 = '+'; c3 = '+';}
                       else if(d < 0.750) {c1 = '-'; c2 = '+'; c3 = '-';}
                       else if(d < 0.875) {c1 = '-'; c2 = '-'; c3 = '+';}
                       else               {c1 = '-'; c2 = '-'; c3 = '-';}
                       
                       if(didTML) strcat(tmBuf, ",\n");
                       sprintf(buffer, "  P%d(%cG%d%c(", geneNum, c1, nextGene, c2);
                       
                       /* maintain power law distribution bins */
                       if(bin[nodeDegree[geneNum]] > 0) bin[nodeDegree[geneNum]]--;
                       
                           nodeDegree[geneNum] += i;
                       bin[nodeDegree[geneNum]]++;
                           nodeDegree[nextGene] = i;
                       bin[nodeDegree[nextGene]]++;
                       nextGene++;
                       
                       for(j=0; j<i-1; j++)
                       {
                         if(j) strcat(buffer, ",");
                         sprintf(tempBuf, "G%d", nextGene);
                         strcat(buffer, tempBuf);
                         
                             nodeDegree[nextGene] += 2;
                         bin[nodeDegree[nextGene]]++;
                         nextGene++;
                       }
                       strcat(tmBuf, buffer);
                       sprintf(tempBuf, ")%c)", c3);
                       strcat(tmBuf, tempBuf);
        
                       poolSize += i;
                     }
                     else
                     {
                       /* add sim */
sim:
                       if(bin[1] + i > (int)ceil((double)numGenes*P(1)))
                         break;
            
                       if(drand48() > 0.5) c1 = '+';
                       else                c1 = '-';
                       
                       if(didTML) strcat(tmBuf, ",\n");
                       sprintf(buffer, "  P%d(%cG%d,", geneNum, c1, nextGene);
                       
                       /* maintain power law distribution bins */
                       if(bin[nodeDegree[geneNum]] > 0) bin[nodeDegree[geneNum]]--;
                       
                           nodeDegree[geneNum] += i;
                       bin[nodeDegree[geneNum]]++;
                           nodeDegree[nextGene] += 1;
                       bin[nodeDegree[nextGene]]++;
                       nextGene++;
                       
                       for(j=0; j<i-1; j++)
                       {
                         if(j) strcat(buffer, ",");
                         sprintf(tempBuf, "G%d", nextGene);
                         strcat(buffer, tempBuf);
                         
                             nodeDegree[nextGene] += 1;
                         bin[nodeDegree[nextGene]]++;
                         nextGene++;
                       }
                       strcat(tmBuf, buffer);
                       strcat(tmBuf, ")");
                       
                       poolSize += i;
                     }
                     
                     didTML = 1;
                     goto outtaHere;
          }
        }
      }
    }
outtaHere:
    p2 = p+1;
  }
}

int BuildDOR(void)
{
  int i, j, didOne=0, moreDOR, neededBufSz=0, numTFsUsed=0, tmp, totGenes=0;
  char c1, buf2[32];
  void *p;
  
  bufSz = 128;
  buf = (char *) malloc(bufSz);
  if(!buf)
  {
    fprintf(stderr, "range: malloc error, exiting...\n");
    exit(1);
  }
  
  done = 0;
  strcpy(buf, "DOR(\n");

  /* add regulated genes to DOR */
  while(numAtRegulatedGeneDegree>0)
  {
    for(i=0; i<numAtRegulatedGeneDegree; i++)
    {
      neededBufSz += 32;
      CheckBufSz(neededBufSz);
      if(((nextGene+1) < tfList[tfListSize-regulatedGeneDegree]) &&
        (bin[nodeDegree[tfList[tfListSize-regulatedGeneDegree]]+1] <
        (int)ceil((double)numGenes*
           P(nodeDegree[tfList[tfListSize-regulatedGeneDegree]]+1))))
      {
        if(didOne && buf[strlen(buf)-2]!=',')
        {
          buf[strlen(buf)-1] = ',';
          strcat(buf, "\n");
        }
        sprintf(buf2, "  G%d(", nextGene++);
        strcat (buf, buf2);
        didOne = 1;
        totGenes++;
      }
      else
      {
        strcat(buf, ")\n");
        done = 1;
        
        if(totGenes < 2) nextGene--;
        if(totGenes>1)
          return 1;
        else
          return 0;
      }
          
      /* fill in the TFs for the gene */
      for(j=0; j<regulatedGeneDegree; j++)
      {
        neededBufSz += 32;
        CheckBufSz(neededBufSz);
        
        tmp = tfListSize - regulatedGeneDegree + j;
        if(nextGene >= tfList[tmp]) /* won't happen for j=0 */
        {
          strcat(buf, ")\n)\n");
          done = 1;
          
          if(totGenes < 2) nextGene--;
          if(totGenes>1)
            return 1;
          else
            return 0;
        }
            
        if(j && (bin[nodeDegree[tfList[tmp]]+1] > (int)ceil((double)numGenes*
           P(nodeDegree[tfList[tmp]]+1))))
          break;
        
        if(j) strcat(buf, ",");
          
        if(drand48() > 0.5) c1 = '+';
        else                c1 = '-';

        sprintf(buf2, "P%d%c", tfList[tmp], c1);
        strcat (buf, buf2);
      }
      
      if(j <= numTFsUsed) moreDOR = 0;
      else                moreDOR = 1;
      numTFsUsed = j;

      if(regulatedGeneDegree == 1) /* many motifs are of degree two */
      {
        if(((i+1) < numAtRegulatedGeneDegree ||
          ((int)ceil((double)poolSize * P(regulatedGeneDegree+1))/2 > 0)) &&
          (nextGene != numGenes - tfListSize-1) && moreDOR)
          strcat(buf, "),\n");
        else
        {
          strcat(buf, ")\n");
          break;
        }
      }
      else
      {
        if(((i+1) < numAtRegulatedGeneDegree || 
          ((int)ceil((double)poolSize * P(regulatedGeneDegree+1)) > 0)) &&
          (nextGene != numGenes - tfListSize-1) && moreDOR)
          strcat(buf, "),\n");
        else
        {
          strcat(buf, ")\n");
          break;
        }
      }
    }

    regulatedGeneDegree++;
    numAtRegulatedGeneDegree = (int)floor((double)poolSize * P(regulatedGeneDegree));
    if(regulatedGeneDegree == 2) /* many motifs are of degree two */
      numAtRegulatedGeneDegree /= 2;
    
    if(nextGene == numGenes - tfListSize-1)
    {
      strcat(buf, ")\n");
      done = 1;
      
      if(totGenes < 2) nextGene--;
      if(totGenes>1)
        return 1;
      else
        return 0;
    }
        
    tfListSize++;
    p = (int *) realloc(tfList, tfListSize*sizeof(int));
    if(!p)
    {
      fprintf(stderr, "range: realloc error, exiting...\n");
      exit(1);
    }
    tfList = p;
    tfList[tfListSize-1] = numGenes - tfListSize;
  }

  strcat(buf, ")\n");
  maxGenesInDOR--;
  regulatedGeneDegree = 1;
  
  /* maxGenesInDOR instead of poolSize for regulatedGeneDegree=1 keeps number of
   * genes with degree one low, as motifs added later will increase the number
   */
  numAtRegulatedGeneDegree = (int)floor((double)maxGenesInDOR*P(regulatedGeneDegree));

  if(totGenes < 2)
  {
    nextGene--;
    done = 1;
  }
  
  if(totGenes>1)
    return 1;
  else
    return 0;
}

void CheckBufSz(int size)
{
  void *p;
  
  if(size >= bufSz)
  {
    p = (char *) realloc(buf, size+128);
    if(!p)
    {
      fprintf(stderr, "range: CheckBufSz realloc error, exiting...\n");
      exit(1);
    }
    buf = p;
    bufSz = size+128;
  }
}

double P(int degree)
{
  return 0.7/(double)(degree*degree);
}

void UpdatePowerLawDistributionForDOR(void)
{
  int i, j, degree, geneNum;
  char buffer[128], *p, *p2;
  
  p2 = buf;
  while(p=strstr(p2, "G"))
  {
    /* what is the degree of this node? */
    i = degree = 0;
    while((p+1)[i] != 'G' && i<strlen(p+1))
    {
      if((p+1)[i] == 'P')
      {
        degree++;
        j = 0;
        buffer[0] = 0x0;
        while(isdigit((p+1)[i+j+1]))
          buffer[j] = (p+1)[i+j+++1];
        buffer[j] = 0x0;
        geneNum = atoi(buffer);

        if(bin[nodeDegree[geneNum]] > 0) bin[nodeDegree[geneNum]]--;
        
            nodeDegree[geneNum]++;
        bin[nodeDegree[geneNum]]++;
      }
      i++;
    }
    
    /* which gene is this? */
    i = 0;
    buffer[0] = 0x0;
    while(isdigit((p+1)[i]))
      buffer[i] = (p+1)[i++];
    buffer[i] = 0x0;
    
    geneNum = atoi(buffer);
    
    /* update power-law distribution */
    nodeDegree[geneNum] = degree;
    bin[degree]++;
    
    p2 = p+1;
  }
}

void UpdatePowerLawDistributionForTML(void)
{
  int i, j, degree, firstG, geneNum, motif, numG, numSigns, proteinNum;
  char buffer[128], *p, *p2;

  p2 = buf;
  while(p=strstr(p2, "P"))
  {
    /* what protein is this? */
    i = 0;
    buffer[0] = 0x0;
    while(isdigit((p+1)[i]))
      buffer[i] = (p+1)[i++];
    buffer[i] = 0x0;
        
    proteinNum = atoi(buffer);
    
    /* what degree gets added to this node? 
     * first, is this motif an FFL, a multiFFL, or a sim? 
     */
    i = numG = numSigns = 0;
    while((p+1)[i] != 'P' && i<strlen(p+1))
    {
      if(((p+1)[i] == '+') || ((p+1)[i] == '-'))
        numSigns++;
      if((p+1)[i] == 'G')
        numG++;
      i++;
    }
    
    if(numSigns==3 && numG==2)     motif = 0;
    else if(numSigns==3 && numG>2) motif = 1;
    else if(numSigns==1)           motif = 2;
    else fprintf(stderr, "range: error, unknown motif, failing...\n"); 
 
    switch(motif)
    {
      case 0: /* FFL */
    
              i = 0;
              while((p+1)[i] != 'P' && i<strlen(p+1))
              {
                if((p+1)[i] == 'G')
                {
                  j = 0;
                  buffer[0] = 0x0;
                  while(isdigit((p+1)[i+j+1]))
                    buffer[j] = (p+1)[i+j+++1];
                  buffer[j] = 0x0;
        
                  geneNum = atoi(buffer);
                  
                      nodeDegree[geneNum] += 2;
                  bin[nodeDegree[geneNum]]++;
                }
                i++;
              }
              if(bin[nodeDegree[proteinNum]] > 0) bin[nodeDegree[proteinNum]]--;
        
                  nodeDegree[proteinNum] += 2;
              bin[nodeDegree[proteinNum]]++;

              break;
      
      case 1: /* mFFL */
      
              i = degree = firstG = 0;
              while((p+1)[i] != 'P' && i<strlen(p+1))
              {
                if((p+1)[i] == 'G')
                {
                  degree++;
                  j = 0;
                  buffer[0] = 0x0;
                  while(isdigit((p+1)[i+j+1]))
                    buffer[j] = (p+1)[i+j+++1];
                  buffer[j] = 0x0;
        
                  geneNum = atoi(buffer);
                  if(firstG)
                  {
                        nodeDegree[geneNum] += 2;
                    bin[nodeDegree[geneNum]]++;
                  }
                  else
                    firstG = geneNum; /* should never be 0 */
                }
                i++;
              }
              if(bin[nodeDegree[proteinNum]] > 0) bin[nodeDegree[proteinNum]]--;
        
                  nodeDegree[proteinNum] += degree;
              bin[nodeDegree[proteinNum]]++;
        
                  nodeDegree[firstG] = degree;
              bin[nodeDegree[firstG]]++;

              break;
      
      case 2: /* sim */
    
              i = degree = 0;
              while((p+1)[i] != 'P' && i<strlen(p+1))
              {
                if((p+1)[i] == 'G')
                {
                  degree++;
                  j = 0;
                  buffer[0] = 0x0;
                  while(isdigit((p+1)[i+j+1]))
                    buffer[j] = (p+1)[i+j+++1];
                  buffer[j] = 0x0;
        
                  geneNum = atoi(buffer);
        
                  nodeDegree[geneNum] = 1;
                  bin[1]++;
                }
                i++;
              }
              if(bin[nodeDegree[proteinNum]] > 0) bin[nodeDegree[proteinNum]]--;
   
                  nodeDegree[proteinNum] += degree;
              bin[nodeDegree[proteinNum]]++;

              break;
      
      default: fprintf(stderr, "range: error2, unknown motif, failing...\n"); 
    }
    p2 = p+1;
  }
}
