%{
/* NEMO - NEtwork MOtif language
 *
 * commands:
 * $ yacc -d nemo.y (or bison -y -d nemo.y)
 * $ lex nemo.lex   (or flex nemo.lex)
 * $ gcc -o nemo2sbml lex.yy.c y.tab.c -ll -ly -lm -lsbml 
 * (or 
 * $ gcc -o nemo2sbml lex.yy.c y.tab.c -lfl -lm -lsbml for flex/bison)
 * $ ./nemo2sbml
 *
 * sample input:
 * [DOR(G1(P1+,P2-), G2(P1-,P3+)), TMLIST(P1(+G3-G4+), P1(+G6,G7,G8), P2(-G9+(G10,G11,G12)+)), GLIST(G5(P4-,P5+), G13(P13+,P7+,P8-))]
 * ^d to stop
 *
 * or save input to a file and do
 * $ cat file | ./nemo2sbml
 *
 * if using range to generate the network, do something like
 * $ ./range 1000 | ./nemo2sbml
 *
 * Lexer for a grammar to describe transcription network motifs as enumerated in
 * "An Introduction to Systems Biology - Design Principles of Biological 
 * Circuits" by Uri Alon, as well as their interconnection hierarchy.
 *
 * Copyright (C) 2007, University of Alaska Fairbanks
 * Biotechnology Computing Research Group
 * Author: James Long
 *-------------------------------------------------------------------------------
 * NEMO BSD License
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
 * NEMO GPL License
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

#include "y.tab.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int add_GENE(char *);
int add_PROTEIN(char *);
int proteinsOK(void);

int lineNum=1;
%}

%%

DOR     {
          strcpy(yylval.string, yytext);
          return(DOR);
        }
G[0-9]* {
          if(add_GENE(yytext))
          {
            strcpy(yylval.string, yytext);
            return(GENE);
          }
        }
GLIST   {
          strcpy(yylval.string, yytext);
          return(GLIST);
        }
TMLIST  {
          strcpy(yylval.string, yytext);
          return(TMLIST);
        }
P[0-9]* {
          /* naming convention: 'P' is followed by 
           * the number of the gene that makes it.
           */
          if(add_PROTEIN(yytext))
          {
            strcpy(yylval.string, yytext);
            return(PROTEIN);
          }
        }
[\]]    {
          if(proteinsOK())
            return yytext[0];
          else
          {
            yytext[0] = 'E';  /* dummy variable */
            return yytext[0]; /* throws a syntax error */
          }
        }
[0-9]*  {
          strcpy(yylval.string, yytext);
          return(DIGITS);
        }
abs     {
          strcpy(yylval.string, yytext);
          return(ABS);
        }
arccos  {
          strcpy(yylval.string, yytext);
          return(ARCCOS);
        }
arcsin  {
          strcpy(yylval.string, yytext);
          return(ARCSIN);
        }
arctan  {
          strcpy(yylval.string, yytext);
          return(ARCTAN);
        }
ceiling {
          strcpy(yylval.string, yytext);
          return(CEILING);
        }
cos     {
          strcpy(yylval.string, yytext);
          return(COS);
        }
exp     {
          strcpy(yylval.string, yytext);
          return(EXP);
        }
floor   {
          strcpy(yylval.string, yytext);
          return(FLOOR);
        }
ln      {
          strcpy(yylval.string, yytext);
          return(LN);
        }
log     {
          strcpy(yylval.string, yytext);
          return(LOG);
        }
power   {
          strcpy(yylval.string, yytext);
          return(POWER);
        }
root    {
          strcpy(yylval.string, yytext);
          return(ROOT);
        }
sin     {
          strcpy(yylval.string, yytext);
          return(SIN);
        }
tan     {
          strcpy(yylval.string, yytext);
          return(TAN);
        }
[ \t]   { ;/* white space doesn't count */}
\n      { lineNum++; }
.       { /* anything but +, -, :, (, ), [, and F throws a syntax error */
          return yytext[0];
        }
%%

/* node for linked list of coding regions and proteins, used to check:
 * 1) A particular coding region (GENE) may only appear once, when all of the 
 * motifs its Ps regulate, and Ps that are input to it, are enumerated.
 * 2) A particular protein (P) may appear more than once, however, 
 * no protein may appear who does not have a gene that makes it.
 */
struct node
{
  char gene[16];
  char prot[16];
  struct node *next;
};

struct node *list = NULL;

/* add gene and make sure it only appears once */
int add_GENE(char *gene)
{
  char gene_suffix[16];
  char prot_suffix[16];
  struct node *pt;
  
  strcpy(gene_suffix, strstr(gene, "G")+1);

  if(!list) /* create the first node and add entry */
  {
    list = (struct node *) malloc(sizeof(struct node));
    if(!list)
    {
      fprintf(stderr, "add_GENE: malloc error, failing...\n");
      return 0;
    }
    strcpy(list->gene, gene);
    list->next = NULL;
    list->prot[0] = 0x0;
    return 1;
  }
  
  /* march down the list */
  pt = list;
  while(strcmp(gene, pt->gene))
  {
    if(pt->prot[0]) /* has the gene's protein been entered already? */
    {
      strcpy(prot_suffix, strstr(pt->prot, "P")+1);
      if(!strcmp(gene_suffix, prot_suffix))
      {
        strcpy(pt->gene, gene);
        return 1;
      }
    }
    if(pt->next)
      pt = pt->next; 
    else /* at the last node and no match, so add new GENE at head */
    {
      pt = (struct node *) malloc(sizeof(struct node));
      if(!pt)
      {
        fprintf(stderr, "add_GENE: malloc error, failing...\n");
        return 0;
      }
      strcpy(pt->gene, gene);
      pt->next = list;
      pt->prot[0] = 0x0;
      list = pt;
      return 1;
    }
  }

  /* must be a match if we got this far */
  fprintf(stderr, "Error: %s must only appear once.\n", gene);
  return 0;
}


int add_PROTEIN(char *prot)
{
  char gene_suffix[16];
  char prot_suffix[16];
  struct node *pt;
  
  strcpy(prot_suffix, strstr(prot, "P")+1);

  if(!list) /* create the first node and add entry */
  {
    list = (struct node *) malloc(sizeof(struct node));
    if(!list)
    {
      fprintf(stderr, "add_PROTEIN: malloc error, failing...\n");
      return 0;
    }
    strcpy(list->prot, prot);
    list->gene[0] = 0x0;
    list->next = NULL;
    return 1;
  }

  /* march down the list */
  pt = list;
  while(strcmp(prot, pt->prot))
  {
    if(pt->gene[0]) /* has the protein's gene been entered already? */
    {
      strcpy(gene_suffix, strstr(pt->gene, "G")+1);
      if(!strcmp(prot_suffix, gene_suffix))
      {
        strcpy(pt->prot, prot);
        return 1;
      }
    }
    if(pt->next)
      pt = pt->next; 
    else /* at the last node and no match, so add new P */
    {
      pt = (struct node *) malloc(sizeof(struct node));
      if(!pt)
      {
        fprintf(stderr, "add_PROTEIN: malloc error, failing...\n");
        return 0;
      }
      strcpy(pt->prot, prot);
      pt->gene[0] = 0x0;
      pt->next = list;
      list = pt;
      return 1;
    }
  }
  return 1;
}

/* free the gene and protein list */
void free_list(void)
{
  struct node *p0, *p1;

  p0 = list;
  while(p0)
  {
    p1 = p0->next;
    free(p0);
    p0 = p1;
  }

  list = NULL;
}

/* make sure each protein mentioned has a gene that makes it */
int proteinsOK(void)
{
  struct node *pt;
  
  /* for each node's protein, look for a matching gene */
  pt = list;
  
next_node:
  while(pt)
  {
    if(pt->prot[0])
    {
      if(pt->gene[0]) /* OK */
      {
        pt = pt->next;
        goto next_node;
      }
      else
      {
        fprintf(stderr, "Error: %s has no parent GENE, failing...\n", pt->prot);
        free_list();
        return 0; 
      }
    }
    pt = pt->next;
  }
  free_list();
  return 1;
}

void yyerror(char *s)
{
  static int timesCalled=0;
  
  timesCalled++;
  if(timesCalled > 10) return;
  
  fprintf(stderr, "line %4d: %s at '%s'\n", lineNum, s, yytext);
}


