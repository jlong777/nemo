/* range.c                                            10/2007 jlong@jimlong.org
 *
 * RANGE - RAndom Network GEnerator
 * print to stdout a valid transcription network in the NEMO (NEtwork
 * MOtif) language, which if piped into nemo2sbml will output SBML for
 * input to a biochemical simulator, such as COPASI.
 *
 * compile: gcc -o range range.c -lm
 *
 * usage ./range [options] <number of nodes in network>
 *               -d print the node degree distribution as well
 *               -h --help
 *               -n print the node degree for each node
 *               -s <seedval> set the seed for drand48
 *               -v print version
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

#define MAX_GENES     16000  /* max genes in network */
#define MIN_GENES       100
#define VERSION        "1.7" /* code unchanged from 1.4 */

void AddMotifs(void);
void BuildDOR(void);
double P(int);

int *bin, degreeTF, *idealNodeDegree, firstDORGene, firstTF, 
    *masterGeneRegulators, maxDegree, maxNonFatTailDegree, nextGene,
    *nodeDegree, numDOR, numGenes, numSkipped=0, whichDOR;
char c1, c2, c3, Network[1000000];

int main(int argc, char *argv[])
{
  int i, j, k, numGenesAccomodated, option, printDistribution=0, 
      printNodeDegrees=0, tmp;
  long seedval;

  if(argc == 1)
  {
    printf("usage: range [options] <number of genes in network, >= %d, <= %d>\n", MIN_GENES, MAX_GENES);
    printf("              -d print the node degree distribution as well\n");
    printf("              -h --help\n");
    printf("              -n print the node degree for each node\n");
    printf("              -s <seedval> set the seed for drand48\n");
    printf("              -v print version\n");
    return 1;
  }

  seedval = 123456789;
  srand48(seedval);

  /* options parsing */
  while((option = getopt(argc, argv, "s:dhnv")) > 0)
  {
    switch(option)
    {
      case 'd':
        printDistribution = 1;
        break;

      case 'h':
        printf("print to stdout a valid random transcription\n");
        printf("network in the NEMO (NEtwork MOtif) language\n");
        printf("usage: range [options] <number of genes in network, >= %d, <= %d>\n", MIN_GENES, MAX_GENES);
        printf("              -d print the node degree distribution as well\n");
        printf("              -h --help\n");
        printf("              -n print the node degree for each node\n");
        printf("              -s <seedval> set the seed for drand48\n");
        printf("              -v print version\n");
        return 0;

      case 'n':
        printNodeDegrees = 1;
        break;

      case 's':
        for(i=0; i<strlen(optarg); i++)
          if(!isdigit(optarg[i]))
          {
            fprintf(stderr, "range: the \"seedval\" argument (%s) must be a number, returning...\n", optarg);
            return 1;
          }
        seedval = atol(optarg);
        srand48(seedval);
        break;

      case 'v':
        printf("ver %s\n", VERSION);
        return 0;

      default: 
        break;
    }
  }

  numGenes = atoi(argv[argc-1]);

  if((numGenes <= MAX_GENES) && (numGenes >= MIN_GENES))
  {
    /* Compute the maximum degree of a node. For our power-law distribution,
     * the probability that a node has k edges is P(k) = 0.62/(k^2). The max
     * degree then is that k for which we expect to find just one node with
     * degree k. The node with this degree is our master regulator, regulating
     * the first transcription factor (TF) in each DOR, and being regulated
     * by regulated DOR genes, at least one per DOR.
     */

    maxDegree = (int)floor(sqrt(0.62*(double)numGenes) + 0.5); /* round to nearest int */

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

    /* now compute the ideal node degrees based on a power-law distribution,
     * will use this to compute number of DORs, & hence length of "fat-tail"
     */
    idealNodeDegree = (int *) malloc(numGenes*sizeof(int));
    if(!idealNodeDegree)
    {
      fprintf(stderr, "range: malloc error, returning...\n");
      return 1;
    }

    for(i=0; i<numGenes; i++)
      idealNodeDegree[i] = 1;

    k = 0;
    for(i=maxDegree; i>0; i--)
      for(j=0; j<(int)floor(numGenes*P(i)+0.5); j++)
      {
        if(k == numGenes) break;
        idealNodeDegree[k++] = i;
      }

    /* High-level view of the algorithm for construction the network:
     * Iteratively create a network backbone composed of DORs with network
     * motifs attached to the regulated genes of the DOR. The first TF for each
     * DOR, along with the master regulator, constitute the "fat-tail" of the
     * power-law distribution. Each first TF in a DOR will regulate one gene
     * of the DOR, the other TFs in the DOR, and itself.
     */

    /* compute the max number of DORs needed to accomodate all the genes */
    numDOR = -1;
    for(i=4; i<maxDegree; i++) /* degree of DOR TFs */
    {
      numGenesAccomodated = 1;   /* the master regulator */
      for(j=1; j<numGenes; j++)
      {
        if(i >= idealNodeDegree[j]) break;
        tmp = idealNodeDegree[j] - 4; /* number of regulated TFs */
        numGenesAccomodated += i*tmp*(i-2) + tmp + 1;

        if((numGenesAccomodated > numGenes+maxDegree-j) &&
          (j<maxDegree)) /* need extra for master regulator feedback */
        {
          if(i > degreeTF) 
          {
            if(j > numDOR)
            {
              numDOR = j; 
              degreeTF = i;
              maxNonFatTailDegree = idealNodeDegree[j] - 1;
            }
          }
          break;
        }
      }
    }

    /* do less because often the motifs use up the genes */
    if(numDOR > 3) numDOR /=2;

    if(numDOR < 0)
    {
      numDOR = 3;
      degreeTF = 4;
      maxNonFatTailDegree = idealNodeDegree[numDOR] - 1;
    }

    /* there are maxDegree-numDOR that regulate the master */
    masterGeneRegulators = (int *) malloc((maxDegree-numDOR)*sizeof(int));
    if(!masterGeneRegulators)
    {
      fprintf(stderr, "range: malloc error, returning...\n");
      return 1;
    }

    /* first issue a GLIST of the DOR TFs, build one for master later */
    printf("[\nGLIST(\n");

    nextGene = 1;

    /* construct the network */
    whichDOR   = 0;
    Network[0] = 0x0;
    for(i=0; i<numDOR; i++)
    {
      if(nextGene >= numGenes-idealNodeDegree[i+1]-4)
      {
        /* this reduces the master gene degree, fortunately its rare */
        numSkipped = numDOR - i;
        break;
      }
      whichDOR++;
      if(i) printf(",\n");

      if(drand48() > 0.5) c1 = '+';
      else                c1 = '-';
      if(drand48() > 0.5) c2 = '+';
      else                c2 = '-';

      firstTF = nextGene;
      printf("  G%d(P0%c,P%d%c)", firstTF, c1, firstTF, c2); /* 1st TF in DOR */
      nodeDegree[firstTF] += 3;
      nextGene++;

      for(j=0; j<idealNodeDegree[i+1]-4; j++)
      {
        if(drand48() > 0.5) c1 = '+';
        else                c1 = '-';

        printf(",\n");
        printf("  G%d(P%d%c)", nextGene, firstTF, c1); /* other TFs in DOR */
        nodeDegree[firstTF]++;
        nodeDegree[nextGene]++;
        nextGene++;
      }
      BuildDOR();
      AddMotifs();
    }
    printf("\n)\n");

    printf("%s", Network);

    /* finally issue the master regulator */
    if(drand48() > 0.5) c1 = '+';
    else                c1 = '-';

    printf(",\nGLIST(\n  G0(P%d%c", masterGeneRegulators[0], c1);
    for(i=1; i<maxDegree-numDOR-numSkipped; i++)
    {
      if(drand48() > 0.5) c1 = '+';
      else                c1 = '-';

      printf(",P%d%c", masterGeneRegulators[i], c1);
    }
    nodeDegree[0] = maxDegree - 2*numSkipped; /* rare that it gets reduced */
    printf(")");

    /* any leftovers? */
    if(drand48() > 0.5) c1 = '+';
    else                c1 = '-';

    if(nextGene < numGenes) printf(",\n  G%d(P%d%c)", nextGene, nextGene, c1);
    nodeDegree[nextGene] += 2;
    nextGene++;

    while(nextGene < numGenes)
    {
      if(drand48() > 0.5) c1 = '+';
      else                c1 = '-';

      printf(",\n  G%d(P%d%c)", nextGene, nextGene, c1);
      nodeDegree[nextGene] += 2;
      nextGene++;
    }

    printf("\n)\n]\n");

    if(printDistribution)
    {
      for(i=1; i<maxDegree+1; i++)
        for(j=0; j<numGenes; j++)
          if(nodeDegree[j] == i)
            bin[i]++;

      printf("\nNode Degree Distribution:\n\n");
      for(i=1; i<maxDegree+1; i++)
        printf("%3d = %4d\n", i, bin[i]);
    }

    if(printNodeDegrees)
    {
      printf("\nNode Degrees:\n");
      for(i=0; i<numGenes; i++)
      {
        if(!(i%5)) printf("\n");
        printf("%7d = %3d", i, nodeDegree[i]);
      }
    }
    printf("\n");
  }
  else
  {
    printf("usage: range [options] <number of genes in network, >= %d, <= %d>\n", MIN_GENES, MAX_GENES);
    printf("              -d print the node degree distribution as well\n");
    printf("              -h --help\n");
    printf("              -n print the node degree for each node\n");
    printf("              -s <seedval> set the seed for drand48\n");
    printf("              -v print version\n");
    return 0;
  }


  return 0;
}

void AddMotifs(void)
{
  /* increase degree by as much as maxNonFatTailDegree - the nodeDegree */

  int i, j, k, gene, pick;
  char buf[4096];
  double fp;

  if(nextGene+1 >= numGenes)
    return;

  gene = firstDORGene;

  strcat(Network, ",\nTMLIST(\n");

  /* add random motifs, keeping the node degree out of the "fat-tail" */
  for(i=0; i<idealNodeDegree[whichDOR]-4; i++)
  {
    for(j=0; j<degreeTF-1; j++)
    {
      if(!i && !j) /* at a gene who regulates the master, don't add motif */
      {
        break;
      }
      else
      {
        /* if a master regulator, skip */
        for(k=0; k<maxDegree-numDOR; k++)
          if(masterGeneRegulators[k] == gene)
            break;

        /* randomly add an FFL, a multiFFL, or a sim */
        fp = drand48();
        if((fp < 0.4) || (nextGene == numGenes-2))      /* FFL */
        {
          if(nextGene > numGenes-2) break;

          if(drand48() > 0.5) c1 = '+';
          else                c1 = '-';
          if(drand48() > 0.5) c2 = '+';
          else                c2 = '-';
          if(drand48() > 0.5) c3 = '+';
          else                c3 = '-';

          sprintf(buf, "  P%d(%cG%d%cG%d%c", gene, c1, nextGene,   c2,
                                                       nextGene+1, c3);
          strcat(Network, buf);
          nodeDegree[nextGene] += 2;
          nodeDegree[nextGene+1] +=2;
          nodeDegree[gene] += 2;
          nextGene += 2;
          gene++;
        }
        else if(fp < 0.8) /* multiFFL */
        {
	         /* pick from an exp dist the number of genes */
          pick = 3 + (int)floor(-1.0*log(1-drand48()));
          if(pick > maxNonFatTailDegree-nodeDegree[gene]) 
             pick = maxNonFatTailDegree-nodeDegree[gene]; /* clamp */

          if(pick < 3) pick = 3;

          if(nextGene > numGenes-pick) break;

          if(drand48() > 0.5) c1 = '+';
          else                c1 = '-';
          if(drand48() > 0.5) c2 = '+';
          else                c2 = '-';
          if(drand48() > 0.5) c3 = '+';
          else                c3 = '-';

          sprintf(buf, "  P%d(%cG%d%c(", gene, c1, nextGene, c2);
          strcat(Network, buf);
          nodeDegree[nextGene]++;
          nodeDegree[gene]++;
          nextGene++;

	         for(k=0; k<pick-1; k++)
       	  {
            if(nextGene >= numGenes) break;

            if(k) strcat(Network, ",");

       	    sprintf(buf, "G%d", nextGene);
            strcat(Network, buf);
            nodeDegree[nextGene]++;
            nodeDegree[nextGene-k-1]++;
            nodeDegree[gene]++;
            nextGene++;
          }
          sprintf(buf, ")%c", c3);
          strcat(Network, buf);
          gene++;
        }
        else               /* sim */
        {
          /* pick from an exp dist the number of genes */
          pick = 2 + (int)floor(-10.0*log(1-drand48()));
          if(pick > maxNonFatTailDegree-nodeDegree[gene]) 
             pick = maxNonFatTailDegree-nodeDegree[gene]; /* clamp */

          if(pick < 2) pick = 2;

          if(nextGene > numGenes-pick) break;

          if(drand48() > 0.5) c1 = '+';
          else                c1 = '-';

          sprintf(buf, "  P%d(%c", gene, c1);
          strcat(Network, buf);

	         for(k=0; k<pick; k++)
       	  {
            if(nextGene >= numGenes) break;

            if(k) strcat(Network, ",");

       	    sprintf(buf, "G%d", nextGene);
            strcat(Network, buf);
            nodeDegree[nextGene]++;
            nodeDegree[gene]++;
            nextGene++;
	         }
          gene++;
        }
        
        if((nextGene < numGenes-1) && 
          !((i == idealNodeDegree[whichDOR]-5) && (j == degreeTF-2) &&
          (maxNonFatTailDegree-nodeDegree[gene] < 2)))
          strcat(Network, "),\n");
        else
          strcat(Network, ")\n");
      }
    }
  }

  if(nextGene > numGenes-2) goto outahere;

  /* if a master regulator, skip */
  for(k=0; k<maxDegree-numDOR; k++)
    if(masterGeneRegulators[k] == gene)
      goto outahere;

  /* at the last gene, can increase degree by as
   * much as maxNonFatTailDegree - the nodeDegree
   */
	
  /* add a sim */
  k = maxNonFatTailDegree-nodeDegree[gene];
  if(k < 2) goto outahere;
  sprintf(buf, "  P%d(%c", gene, c1);
  strcat(Network, buf);
  for(i=0; i<k; i++)
  {
    if(nextGene >= numGenes) break;

    if(i) strcat(Network, ",");

    sprintf(buf, "G%d", nextGene);
    strcat(Network, buf);
    nodeDegree[nextGene]++;
    nodeDegree[gene]++;
    nextGene++;
  }

  strcat(Network, ")\n");
outahere:
  strcat(Network, ")\n");
}

void BuildDOR(void)
{
  int i, j, gene, haveLastParen, lastGene;
  static int masterGeneRegulatorIndex=0;
  char buf[32], buf01[1024];

  if(nextGene >= numGenes-1) return;

  strcat(Network, ",\nDOR(\n");

 /* how many genes does the first DOR have that regulates the master? 
  * there are maxDegree-numDOR that regulate the master; the other numDOR-1
  * DORs each contribute one, so first DOR will contribute
  * (maxDegree-numDOR) - (numDOR-1) = maxDegree - 2*numDOR + 1
  */

  if(drand48() > 0.5) c1 = '+';
  else                c1 = '-';

  if(drand48() > 0.5) c2 = '+';
  else                c2 = '-';

  /* the 1st gene is regulated by the 1st and 2nd TFs */
  gene = 1;
  masterGeneRegulators[masterGeneRegulatorIndex++] = nextGene;
  firstDORGene = nextGene;
  sprintf(buf, "  G%d(P%d%c", nextGene, firstTF, c1); /* 1st TF */
  strcpy(buf01, buf);
  
  nodeDegree[firstTF]++;
  nodeDegree[nextGene]++;
  lastGene = nextGene++;
  gene++;
  haveLastParen = 0;

  for(i=0; i<idealNodeDegree[whichDOR]-4; i++) /* the remaining TFs */
  {
    for(j=0; j<degreeTF-1; j++)                /* edges for TFs */
    {
      if(drand48() > 0.5) c1 = '+';
      else                c1 = '-';

      if(!j)
      {
        if(nextGene >= numGenes) break;
        sprintf(buf, ", P%d%c", firstTF+i+1, c1);
        if(!i) strcat(buf01, buf);

        /* increase the degree of the firstDORGene as much as possible */
        if((i && nodeDegree[firstDORGene] < maxNonFatTailDegree) &&
                (nodeDegree[firstTF+i+1]  < maxNonFatTailDegree))
        {
          strcat(buf01, buf);
          nodeDegree[firstDORGene]++;
          nodeDegree[firstTF+i+1]++;
        }
        if(i) strcat(buf, ")");

        nodeDegree[firstTF+i+1]++;
        nodeDegree[lastGene]++;
        haveLastParen = 1;
      }
      else if(j<degreeTF-2)
      {
        if(nextGene >= numGenes) break;

        if((whichDOR == 1) && (gene <= maxDegree - 2*numDOR + 1))
          masterGeneRegulators[masterGeneRegulatorIndex++] = nextGene;

        sprintf(buf, "  G%d(P%d%c)", nextGene, firstTF+i+1, c1);
        nodeDegree[firstTF+i+1]++;
        nodeDegree[nextGene]++;
        gene++;
        nextGene++;
        haveLastParen = 1;
      }
      else
      {
        if(nextGene >= numGenes) break;

        if((whichDOR == 1) && (gene <= maxDegree - 2*numDOR + 1))
          masterGeneRegulators[masterGeneRegulatorIndex++] = nextGene;

        sprintf(buf, "  G%d(P%d%c",  nextGene, firstTF+i+1, c1);
        nodeDegree[firstTF+i+1]++;
        nodeDegree[nextGene]++;
        lastGene = nextGene++;
        gene++;
        haveLastParen = 0;
      }
      
      if(i || j)
        strcat(Network, buf);

      if((i && !j) || ((j<degreeTF-2) && j))
        strcat(Network, ",\n");
    }
  }
  
  if(!haveLastParen) strcat(Network, "),\n");
  
  strcat(Network, buf01);
  strcat(Network, ")\n)\n");
}

double P(int degree)
{
  return 0.62/(double)(degree*degree);
}

