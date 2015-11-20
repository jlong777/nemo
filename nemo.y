%{
/* NEMO - NEtwork MOtif language
 * A grammar to recognize transcription network motifs as enumerated in
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

#define VERSION "1.7"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sbml/SBMLTypes.h>
#include <sys/types.h>
#include <regex.h>
#include <unistd.h>


#define BUFSZ      256
#define GENES    16384
#define NLT_SZ    8192
#define NON_LINEAR     /* comment this out if you don't want non-linear terms in the Hill function */
#define RES          2 /* the max number of regex subexpressions allowed */
#define SBML_LEVEL   2
#define SBML_VERSION 1

int check_dor(char *);
void mark_neighbors(int);
char * randomGeneralizedHill(char *, char *);
char * insertNonLinearTerms(char *, char *);
char * explicitKineticLaw(char *, char *, char *);
void xgmmlXML(char *, char *);

int cytoBufSz, firstP=1, kineticLawInfo=0, num_files=0, num_genes, num_sgn=0,
    parameterIndex=0, parseInfo=0, rand_func=0, tot_genes=0, user_func=0, xgmml=0;
char *cytoBuf, docbuf[64], genes[GENES][BUFSZ], followingGene[BUFSZ], *kineticLawString, *kLSp,
     marked[GENES], modelname[BUFSZ], output[BUFSZ], *p, *pt, *pgp, pg[BUFSZ], protein[BUFSZ],
     returnString[NLT_SZ], sgn0, sgn1, sgn2, temp[BUFSZ], *tmp, tmpCytoBuf[BUFSZ], 
     transcriptionFactors[8*BUFSZ], xgmmlTmp[BUFSZ];

extern FILE *yyin, *yyout;
FILE *cyto_graph;
const char *sid = "Cell"; /* compartment name */
regex_t preg;

/* SBML */
Compartment_t              *compart;
KineticLaw_t               *kl;
Model_t                    *model;
ModifierSpeciesReference_t *msr;
Parameter_t                *dparam;
Parameter_t                *param;
Reaction_t                 *react;
SBMLDocument_t             *doc;
Species_t                  *species;
SpeciesReference_t         *reactant;
Unit_t                     *unit;
UnitDefinition_t           *unitdef;
%}

%union
{
  char string[32];
  char *string_pt;
}

%token <string> ABS ARCCOS ARCSIN ARCTAN CEILING COS DIGITS DOR EXP FLOOR
                GENE GLIST LN LOG POWER PROTEIN ROOT SIN TAN TEN TMLIST
%type <string_pt> constant dor expr gene ff_loop gene_list multi_out pg sim sim_list
                  start protein p_error p_list sgn term tmotif tmotif_list tr_group
%left '+' '-'
%left '*' '/'
%nonassoc UMINUS

%%
start       : start '[' tr_group ']'                                                          {
                                                                                                /* output new SBML file */
                                                                                                if(output[0])
                                                                                                  sprintf(docbuf, "%s_%d", output, num_files++);
                                                                                                else
                                                                                                  sprintf(docbuf, "regulatoryNetwork_%dgenes_%d", tot_genes, num_files++);
                                                                                                
                                                                                                Model_setId(model, docbuf);
                                                                                                
                                                                                                if(rand_func && user_func)
                                                                                                  sprintf(modelname, "Synthetic Network: user specified input functions, and randomized parameters in generalized Hill Functions (nemo2sbml ver %s)", VERSION);
                                                                                                else if(rand_func)
                                                                                                  sprintf(modelname, "Synthetic Network: randomized parameters in generalized Hill Functions (nemo2sbml ver %s)", VERSION);
                                                                                                else if(user_func)
                                                                                                  sprintf(modelname, "Synthetic Network: user specified input functions (nemo2sbml ver %s)", VERSION);
                                                                                                else /* probably shouldn't get here */
                                                                                                  sprintf(modelname, "Synthetic Network: no input functions (nemo2sbml ver %s)", VERSION);
                                                                                                  
                                                                                                Model_setName(model, modelname);
                                                                                                
                                                                                                strcat(docbuf, ".xml");
                                                                                                if(writeSBML(doc, docbuf))
                                                                                                  printf("SBML document written: %s\n", docbuf);
                                                                                                else
                                                                                                  fprintf(stderr, "nemo2sbml: Error, failed to write SBML document %s\n", docbuf);

                                                                                                if(xgmml)
                                                                                                {
                                                                                                  /* output xgmml file for cytoscape */
                                                                                                  sprintf(docbuf, "cytoscapeGraph_%dgenes_%d.xgmml", tot_genes, num_files-1);
                                                                                                  cyto_graph = fopen(docbuf, "w");
                                                                                                  if(!cyto_graph)
                                                                                                  {
                                                                                                    fprintf(stderr, "nemo2sbml: Error, failed to open %s for writing, continuing\n", docbuf);
                                                                                                  }
                                                                                                  else
                                                                                                  {
                                                                                                    fprintf(cyto_graph, "<?xml version=\"1.0\"?>\n");
                                                                                                    fprintf(cyto_graph, "<graph label=\"%s\" id=\"0\" xmlns=\"http://www.cs.rpi.edu/XGMML\">\n", docbuf);
                                                                                                    fprintf(cyto_graph, "%s", cytoBuf);
                                                                                                    fprintf(cyto_graph, "</graph>\n");
                                                                                                    fclose (cyto_graph);
                                                                                                    printf("XGMML document written: %s\n", docbuf);
                                                                                                  }
                                                                                                }
                                                                                                
                                                                                                firstP = 1;
                                                                                                tot_genes = 0;
                                                                                                parameterIndex = 0;
                                                                                                rand_func = user_func = 0;
                                                                                                
                                                                                                //SBMLDocument_free(doc); /* why does this cause a segfault? */
                                                                                                
                                                                                                doc   = SBMLDocument_createWith(SBML_LEVEL, SBML_VERSION);
                                                                                                model = SBMLDocument_createModel(doc);
 
                                                                                                compart = Model_createCompartment(model);
                                                                                                Compartment_setId(compart, sid);
                                                                                                Compartment_setVolume(compart, 1.0);

                                                                                                /* cell volume is 10^-15 liters = 1 cubic micrometer = 1 femtoliter*/
                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "volume");
                                                                                                UnitDefinition_setName(unitdef, "femtoliter");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), 1, -15));
  
                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "microM_cell");
                                                                                                UnitDefinition_setName(unitdef, "microMole/cell");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("mole"), 1, -6));
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), -1, -15));

                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "hill_coeff");
                                                                                                UnitDefinition_setName(unitdef, "Hill Coefficient");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("dimensionless"), 1, 0)); /* 1 to 4 */

                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "s_fl");
                                                                                                UnitDefinition_setName(unitdef, "sec/femtoliter");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("second"), 1, 0));
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), -1, -15));

                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "s_mole");
                                                                                                UnitDefinition_setName(unitdef, "sec/microMole");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("second"), 1, 0));
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("mole"), -1, -6));

                                                                                                species = Model_createSpecies(model);
                                                                                                Species_setId(species, "devNull");
                                                                                                Species_setName(species, "devNull");
                                                                                                Species_setCompartment(species, sid);
                                                                                                Species_setInitialConcentration(species, 0.0);
                                                                                                Species_setBoundaryCondition(species, 1);
                                                                                                Species_setConstant(species, 1);
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed network:     [%s]\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "["); strcat(tmp, $3); strcat(tmp, "]");  
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | '[' tr_group ']'                                                                {
                                                                                                /* output new SBML file */
                                                                                                if(output[0])
                                                                                                  sprintf(docbuf, "%s_%d", output, num_files++);
                                                                                                else
                                                                                                  sprintf(docbuf, "regulatoryNetwork_%dgenes_%d", tot_genes, num_files++);
                                                                                                  
                                                                                                Model_setId(model, docbuf);
                                                                                                
                                                                                                if(rand_func && user_func)
                                                                                                  sprintf(modelname, "Synthetic Network: user specified input functions, and randomized parameters in generalized Hill Functions (nemo2sbml ver %s)", VERSION);
                                                                                                else if(rand_func)
                                                                                                  sprintf(modelname, "Synthetic Network: randomized parameters in generalized Hill Functions (nemo2sbml ver %s)", VERSION);
                                                                                                else if(user_func)
                                                                                                  sprintf(modelname, "Synthetic Network: user specified input functions (nemo2sbml ver %s)", VERSION);
                                                                                                else /* probably shouldn't get here */
                                                                                                  sprintf(modelname, "Synthetic Network: no input functions (nemo2sbml ver %s)", VERSION);
                                                                                                
                                                                                                Model_setName(model, modelname);
                                                                                                
                                                                                                strcat(docbuf, ".xml");
                                                                                                if(writeSBML(doc, docbuf))
                                                                                                  printf("SBML document written: %s\n", docbuf);
                                                                                                else
                                                                                                  fprintf(stderr, "nemo2sbml: Error, failed to write SBML document %s\n", docbuf);
 
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  /* output xgmml file for cytoscape */
                                                                                                  sprintf(docbuf, "cytoscapeGraph_%dgenes_%d.xgmml", tot_genes, num_files-1);
                                                                                                  cyto_graph = fopen(docbuf, "w");
                                                                                                  if(!cyto_graph)
                                                                                                  {
                                                                                                    fprintf(stderr, "nemo2sbml: Error, failed to open %s for writing, continuing\n", docbuf);
                                                                                                  }
                                                                                                  else
                                                                                                  {
                                                                                                    fprintf(cyto_graph, "<?xml version=\"1.0\"?>\n");
                                                                                                    fprintf(cyto_graph, "<graph label=\"%s\" id=\"0\" xmlns=\"http://www.cs.rpi.edu/XGMML\">\n", docbuf);
                                                                                                    fprintf(cyto_graph, "%s", cytoBuf);
                                                                                                    fprintf(cyto_graph, "</graph>\n");
                                                                                                    fclose (cyto_graph);
                                                                                                    printf("XGMML document written: %s\n", docbuf);
                                                                                                  }
                                                                                                }
                                                                                                
                                                                                                firstP = 1;
                                                                                                tot_genes = 0;
                                                                                                parameterIndex = 0;
                                                                                                rand_func = user_func = 0;
                                                                                                
                                                                                                //SBMLDocument_free(doc); /* why does this cause a segfault? */
                                                                                                
                                                                                                doc   = SBMLDocument_createWith(2, 1);
                                                                                                model = SBMLDocument_createModel(doc);
 
                                                                                                compart = Model_createCompartment(model);
                                                                                                Compartment_setId(compart, sid);
                                                                                                Compartment_setVolume(compart, 1.0);

                                                                                                /* cell volume is 10^-15 liters = 1 cubic micrometer = 1 femtoliter*/
                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "volume");
                                                                                                UnitDefinition_setName(unitdef, "femtoliter");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), 1, -15));
  
                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "microM_cell");
                                                                                                UnitDefinition_setName(unitdef, "microMole/cell");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("mole"), 1, -6));
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), -1, -15));

                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "hill_coeff");
                                                                                                UnitDefinition_setName(unitdef, "Hill Coefficient");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("dimensionless"), 1, 0)); /* 1 to 4 */

                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "s_fl");
                                                                                                UnitDefinition_setName(unitdef, "sec/femtoliter");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("second"), 1, 0));
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), -1, -15));

                                                                                                unitdef = Model_createUnitDefinition(model);
                                                                                                UnitDefinition_setId(unitdef, "s_mole");
                                                                                                UnitDefinition_setName(unitdef, "sec/microMole");
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("second"), 1, 0));
                                                                                                UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("mole"), -1, -6));

                                                                                                species = Model_createSpecies(model);
                                                                                                Species_setId(species, "devNull");
                                                                                                Species_setName(species, "devNull");
                                                                                                Species_setCompartment(species, sid);
                                                                                                Species_setInitialConcentration(species, 0.0);
                                                                                                Species_setBoundaryCondition(species, 1);
                                                                                                Species_setConstant(species, 1);
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed network:     [%s]\n", $2);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($2) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, "["); strcat(tmp, $2); strcat(tmp, "]");  
                                                                                                $$ = tmp;
                                                                                                free($2);
                                                                                              }
            ; 
tr_group    : tr_group ',' dor                                                                {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tr_group:    %s, %s\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, ","); strcat(tmp, $3); 
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | tr_group ',' GLIST '(' gene_list ')'                                            {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tr_group:    %s, GLIST(%s)\n", $1, $5);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($5) + 4);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, ","); strcat(tmp, $3);
                                                                                                strcat(tmp, "("); strcat(tmp, $5);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($5);
                                                                                              }
            | tr_group ',' TMLIST '(' tmotif_list ')'                                         {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tr_group:    %s, TMLIST(%s)\n", $1, $5);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($5) + 4);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, ","); strcat(tmp, $3);
                                                                                                strcat(tmp, "("); strcat(tmp, $5);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($5);
                                                                                              }
            | dor                                                                             {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tr_group:    %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            | GLIST '(' gene_list ')'                                                         {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tr_group:    GLIST(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | TMLIST '(' tmotif_list ')'                                                      {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tr_group:    TMLIST(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            ;
dor         : DOR '(' gene_list ')'                                                           {
                                                                                                /* the graph of transcription factors and genes in a dor
                                                                                                 * (dense overlapping regulon) is connected, and consists
                                                                                                 * of at least two genes; check for that here.
                                                                                                 */
                                                                                                if(check_dor($3))
                                                                                                {
                                                                                                  if(parseInfo)
                                                                                                    printf("parsed dor:         DOR(%s)\n", $3);
                                                                                                    
                                                                                                  tmp = (char *) malloc(strlen($3) + 6);
                                                                                                  if(!tmp)
                                                                                                  {
                                                                                                    yyerror("malloc error, exiting...");
                                                                                                    YYABORT;
                                                                                                  }
                                                                                                
                                                                                                  strcpy(tmp, "DOR("); strcat(tmp, $3); strcat(tmp, ")");
                                                                                                  $$ = tmp;
                                                                                                  free($3);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("DOR graph error.");
                                                                                                  YYABORT;
                                                                                                }
                                                                                              }
            ;
gene_list   : gene_list ',' gene '(' p_list')'                                                {
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, transcriptionFactors);
                                                                                                  xgmmlXML($3, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                /* instantiate Kinetic Law */
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($3, transcriptionFactors);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                                                                                firstP  = 1;
                                                                                                num_sgn = 0;
                                                                                                transcriptionFactors[0] = 0x0;

                                                                                                if(parseInfo)
                                                                                                  printf("parsed gene_list:   %s, %s(%s)\n", $1, $3, $5);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($5) + 4);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, ","); strcat(tmp, $3); 
                                                                                                strcat(tmp, "("); strcat(tmp, $5);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($5);
                                                                                              }
            | gene_list ',' gene '(' p_list ':' 'F' '(' expr ')' ')'                          {
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, transcriptionFactors);
                                                                                                  xgmmlXML($3, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                /* instantiate Kinetic Law */
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($3, transcriptionFactors, $9);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;
                                                                                                firstP  = 1;
                                                                                                num_sgn = 0;
                                                                                                transcriptionFactors[0] = 0x0;

                                                                                                if(parseInfo)
                                                                                                  printf("parsed gene_list:   %s, %s(%s:F(%s))\n", $1, $3, $5, $9);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($5) + strlen($9) + 8);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, ","); strcat(tmp, $3); 
                                                                                                strcat(tmp, "("); strcat(tmp, $5);  strcat(tmp, ":F(");
                                                                                                strcat(tmp, $9);  strcat(tmp, "))");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($5);
                                                                                                free($9);
                                                                                              }
            | gene '(' p_list')'                                                              {
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, transcriptionFactors);
                                                                                                  xgmmlXML($1, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                /* instantiate Kinetic Law */
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($1, transcriptionFactors);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                                                                                firstP  = 1;
                                                                                                num_sgn = 0;
                                                                                                transcriptionFactors[0] = 0x0;

                                                                                                if(parseInfo)
                                                                                                  printf("parsed gene_list:   %s(%s)\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | gene '(' p_list ':' 'F' '(' expr ')' ')'                                        {
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, transcriptionFactors);
                                                                                                  xgmmlXML($1, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                /* instantiate Kinetic Law */
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($1, transcriptionFactors, $7);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;
                                                                                                firstP  = 1;
                                                                                                num_sgn = 0;
                                                                                                transcriptionFactors[0] = 0x0;

                                                                                                if(parseInfo)
                                                                                                  printf("parsed gene_list:   %s(%s:F(%s))\n", $1, $3, $7);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($7) + 7);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);    strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ":F("); strcat(tmp, $7);  strcat(tmp, "))");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($7);
                                                                                              }
            ;
tmotif_list : tmotif_list ',' tmotif                                                          {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tmotif_list: %s, %s\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, ","); strcat(tmp, $3);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | tmotif                                                                          {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tmotif_list: %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            ;
tmotif      : ff_loop                                                                         {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tmotif:      %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            | multi_out                                                                       {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tmotif:      %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            | sim                                                                             {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed tmotif:      %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            ;
p_list      : p_list ',' protein '+'                                                          {
                                                                                                if(firstP)
                                                                                                {
                                                                                                  strcpy(transcriptionFactors, "+");
                                                                                                  firstP = 0;
                                                                                                }
                                                                                                else
                                                                                                  strcat(transcriptionFactors, "+");

                                                                                                strcat(transcriptionFactors, $3);
                                                                                                strcat(transcriptionFactors, ",");

                                                                                                if(parseInfo)
                                                                                                  printf("parsed p_list:      %s, %s+\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, ","); strcat(tmp, $3); 
                                                                                                strcat(tmp, "+");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | p_list ',' protein '-'                                                          {
                                                                                                if(firstP)
                                                                                                {
                                                                                                  strcpy(transcriptionFactors, "-");
                                                                                                  firstP = 0;
                                                                                                }
                                                                                                else
                                                                                                  strcat(transcriptionFactors, "-");

                                                                                                strcat(transcriptionFactors, $3);
                                                                                                strcat(transcriptionFactors, ",");

                                                                                                if(parseInfo)
                                                                                                  printf("parsed p_list:      %s, %s-\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, ","); strcat(tmp, $3); 
                                                                                                strcat(tmp, "-");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | p_list ',' p_error                                                              {
                                                                                                yyerror("Error: PROTEIN must be followed by '+' or '-'\n"); YYABORT;
                                                                                              }
            | protein '+'                                                                     {
                                                                                                if(firstP)
                                                                                                {
                                                                                                  strcpy(transcriptionFactors, "+");
                                                                                                  firstP = 0;
                                                                                                }
                                                                                                else
                                                                                                  strcat(transcriptionFactors, "+");
 
                                                                                                strcat(transcriptionFactors, $1);
                                                                                                strcat(transcriptionFactors, ",");
  
                                                                                                if(parseInfo)
                                                                                                  printf("parsed p_list:      %s+\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "+");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            | protein '-'                                                                     {
                                                                                                if(firstP)
                                                                                                {
                                                                                                  strcpy(transcriptionFactors, "-");
                                                                                                  firstP = 0;
                                                                                                }
                                                                                                else
                                                                                                  strcat(transcriptionFactors, "-");
  
                                                                                                strcat(transcriptionFactors, $1);
                                                                                                strcat(transcriptionFactors, ",");

                                                                                                if(parseInfo)
                                                                                                  printf("parsed p_list:      %s-\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "-");
                                                                                                $$ = tmp;
                                                                                              }
            | p_error                                                                         {
                                                                                                yyerror("Error: %s must be followed by '+' or '-'\n", $1); YYABORT;
                                                                                              }
            ;
ff_loop     : protein '(' sgn gene sgn gene sgn ')'                                           { /* instantiate 2 Kinetic Laws */
                                                                                                sprintf(temp, "%c%s;", sgn0, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($4, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($4, temp);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                              
                                                                                                sprintf(temp, "%cP%s;%c%s;", sgn1, $4+1, sgn2, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($6, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($6, temp);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed ff_loop:     %s(%s%s%s%s%s)\n", $1, $3, $4, $5, $6, $7);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($4) + strlen($5) + 
                                                                                                                      strlen($6) + strlen($7) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3); 
                                                                                                strcat(tmp, $4); strcat(tmp, $5);  strcat(tmp, $6);
                                                                                                strcat(tmp, $7); strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($4);
                                                                                                free($5);
                                                                                                free($6);
                                                                                                free($7);
                                                                                              }
            | protein '(' 'F' '(' expr ')' ':' sgn gene sgn gene sgn ')'                      { /* instantiate 2 Kinetic Laws */
                                                                                                sprintf(temp, "%c%s;", sgn0, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($9, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($9, temp, $5);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;

                                                                                                sprintf(temp, "%cP%s;%c%s;", sgn1, $9+1, sgn2, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($11, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($11, temp);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed ff_loop:     %s(F(%s):%s%s%s%s%s)\n", $1, $5, $8, $9, $10, $11, $12);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($5) + strlen($8) + strlen($9) + 
                                                                                                                      strlen($10)+ strlen($11)+ strlen($12)+ 7);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);   strcat(tmp, "(F("); strcat(tmp, $5);
                                                                                                strcat(tmp, "):"); strcat(tmp, $8);    strcat(tmp, $9);
                                                                                                strcat(tmp, $10);  strcat(tmp, $11);   strcat(tmp, $12);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($5);
                                                                                                free($8);
                                                                                                free($9);
                                                                                                free($10);
                                                                                                free($11);
                                                                                                free($12);
                                                                                              }
            | protein '(' sgn gene sgn gene sgn ':' 'F' '(' expr ')' ')'                      {
                                                                                                /* instantiate 2 Kinetic Laws */
                                                                                                sprintf(temp, "%c%s;", sgn0, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($4, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($4, temp);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                              
                                                                                                sprintf(temp, "%cP%s;%c%s;", sgn1, $4+1, sgn2, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($6, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($6, temp, $11);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed ff_loop:     %s(%s%s%s%s%s:F(%s))\n", $1, $3, $4, $5, $6, $7, $11);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($4) + strlen($5) + 
                                                                                                                      strlen($6) + strlen($7) + strlen($11)+ 7);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, "(");   strcat(tmp, $3);
                                                                                                strcat(tmp, $4);  strcat(tmp, $5);    strcat(tmp, $6);
                                                                                                strcat(tmp, $7);  strcat(tmp, ":F("); strcat(tmp, $11);
                                                                                                strcat(tmp, "))");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($4);
                                                                                                free($5);
                                                                                                free($6);
                                                                                                free($7);
                                                                                                free($11);
                                                                                              }
            | protein '(' 'F' '(' expr ')' ':' sgn gene sgn gene sgn ':' 'F' '(' expr ')' ')' { /* instantiate 2 Kinetic Laws */
                                                                                                sprintf(temp, "%c%s;", sgn0, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($9, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($9, temp, $5);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                              
                                                                                                sprintf(temp, "%cP%s;%c%s;", sgn1, $9+1, sgn2, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($11, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($11, temp, $16);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed ff_loop:     %s(F(%s):%s%s%s%s%s:F(%s))\n", $1, $5, $8, $9, $10, $11, $12, $16);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($5) + strlen($8) + strlen($9) + 
                                                                                                                      strlen($10)+ strlen($11)+ strlen($12)+ strlen($16)+ 11);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);    strcat(tmp, "(F("); strcat(tmp, $5);
                                                                                                strcat(tmp, "):");  strcat(tmp, $8);    strcat(tmp, $9);
                                                                                                strcat(tmp, $10);   strcat(tmp, $11);   strcat(tmp, $12);
                                                                                                strcat(tmp, ":F("); strcat(tmp, $16);   strcat(tmp, "))");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($5);
                                                                                                free($8);
                                                                                                free($9);
                                                                                                free($10);
                                                                                                free($11);
                                                                                                free($12);
                                                                                                free($16);
                                                                                              }
            ;
multi_out   : protein '(' sgn gene sgn pg ')' sgn ')'                                         { /* instantiate >= 2 Kinetic Laws */
                                                                                                sprintf(temp, "%s%s;", $3, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($4, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($4, temp);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                                                                                
                                                                                                /* 1 Kinetic Law for each gene in pg */
                                                                                                strcpy(pg, $6+1);
                                                                                                pgp = pg;
                                                                                                pt = strsep(&pgp, ",");
                                                                                                do
                                                                                                {
                                                                                                  sprintf(temp, "%sP%s;%s%s;", $5, $4+1, $8, $1);
                                                                                                  
                                                                                                  if(xgmml)
                                                                                                  {
                                                                                                    strcpy(xgmmlTmp, temp);
                                                                                                    xgmmlXML(pt, xgmmlTmp);
                                                                                                  }
                                                                                                  
                                                                                                  kl = KineticLaw_create();
                                                                                                  react = Model_createReaction(model);
                                                                                                  if(p=strstr(pt, ":F("))
                                                                                                  {
                                                                                                    p[strlen(p)-1] = 0x0; /* remove last ")" */
                                                                                                    kLSp = explicitKineticLaw(pt, temp, p+3);
                                                                                                    if(kLSp)
                                                                                                    {
                                                                                                      KineticLaw_setFormula(kl, kLSp);
                                                                                                      free(kLSp);
                                                                                                    }
                                                                                                    else
                                                                                                    {
                                                                                                      yyerror("NULL kineticLawString, exiting...");
                                                                                                      YYABORT;
                                                                                                    }
                                                                                                  
                                                                                                    user_func = 1;
                                                                                                  }
                                                                                                  else
                                                                                                  {
                                                                                                    kLSp = randomGeneralizedHill(pt, temp);
                                                                                                    if(kLSp)
                                                                                                    {
                                                                                                      KineticLaw_setFormula(kl, kineticLawString);
                                                                                                      free(kineticLawString);
                                                                                                    }
                                                                                                    else
                                                                                                    {
                                                                                                      yyerror("NULL kineticLawString, exiting...");
                                                                                                      YYABORT;
                                                                                                    }
                                                                                                      
                                                                                                    rand_func = 1;
                                                                                                  }
                                                                                                
                                                                                                  Reaction_setKineticLaw(react, kl);
                                                                                                
                                                                                                  pt = strsep(&pgp, ",");
                                                                                                }
                                                                                                while(pt);
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed multi_out:   %s(%s%s%s%s)%s)\n", $1, $3, $4, $5, $6, $8);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($4) + strlen($5) + 
                                                                                                                      strlen($6) + strlen($8) + 4);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, $4);  strcat(tmp, $5);  strcat(tmp, $6);
                                                                                                strcat(tmp, ")"); strcat(tmp, $8);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($4);
                                                                                                free($5);
                                                                                                free($6);
                                                                                                free($8);
                                                                                              }
            | protein '(' 'F' '(' expr ')' ':' sgn gene sgn pg ')' sgn ')'                    { /* instantiate >= 2 Kinetic Laws */
                                                                                                sprintf(temp, "%s%s;", $8, $1);
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($9, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($9, temp, $5);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;
                                                                                                
                                                                                                /* 1 Kinetic Law for each gene in pg */
                                                                                                strcpy(pg, $11+1);
                                                                                                pgp = pg;
                                                                                                pt = strsep(&pgp, ",");
                                                                                                do
                                                                                                {
                                                                                                  sprintf(temp, "%sP%s;%s%s;", $10, $9+1, $13, $1);
                                                                                                  
                                                                                                  if(xgmml)
                                                                                                  {
                                                                                                    strcpy(xgmmlTmp, temp);
                                                                                                    xgmmlXML(pt, xgmmlTmp);
                                                                                                  }
                                                                                                  
                                                                                                  kl = KineticLaw_create();
                                                                                                  react = Model_createReaction(model);
                                                                                                  if(p=strstr(pt, ":F("))
                                                                                                  {
                                                                                                    p[strlen(p)-1] = 0x0; /* remove last ")" */
                                                                                                    kLSp = explicitKineticLaw(pt, temp, p+3);
                                                                                                    if(kLSp)
                                                                                                    {
                                                                                                      KineticLaw_setFormula(kl, kLSp);
                                                                                                      free(kLSp);
                                                                                                    }
                                                                                                    else
                                                                                                    {
                                                                                                      yyerror("NULL kineticLawString, exiting...");
                                                                                                      YYABORT;
                                                                                                    }
                                                                                                
                                                                                                    user_func = 1;
                                                                                                  }
                                                                                                  else
                                                                                                  {
                                                                                                    kLSp = randomGeneralizedHill(pt, temp);
                                                                                                    if(kLSp)
                                                                                                    {
                                                                                                      KineticLaw_setFormula(kl, kineticLawString);
                                                                                                      free(kineticLawString);
                                                                                                    }
                                                                                                    else
                                                                                                    {
                                                                                                      yyerror("NULL kineticLawString, exiting...");
                                                                                                      YYABORT;
                                                                                                    }
                                                                                                
                                                                                                    rand_func = 1;
                                                                                                  }
                                                                                                  
                                                                                                  Reaction_setKineticLaw(react, kl);
                                                                                                
                                                                                                  pt = strsep(&pgp, ",");
                                                                                                }
                                                                                                while(pt);
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed multi_out:   %s(F(%s):%s%s%s%s)%s)\n", $1, $5, $8, $9, $10, $11, $13);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($5) + strlen($8) + strlen($9) + 
                                                                                                                      strlen($10)+ strlen($11)+ strlen($13)+ 8);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);   strcat(tmp, "(F("); strcat(tmp, $5);
                                                                                                strcat(tmp, "):"); strcat(tmp, $8);    strcat(tmp, $9);
                                                                                                strcat(tmp, $10);  strcat(tmp, $11);   strcat(tmp, ")");
                                                                                                strcat(tmp, $13);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($5);
                                                                                                free($8);
                                                                                                free($9);
                                                                                                free($10);
                                                                                                free($11);
                                                                                                free($13);
                                                                                              }
            ;
pg          : '(' gene                                                                        {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed pg:          (%s\n", $2);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($2) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, "("); strcat(tmp, $2);
                                                                                                $$ = tmp;
                                                                                                free($2);
                                                                                              }
            | '(' gene ':' 'F' '(' expr ')'                                                   {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed pg:          (%s:F(%s)\n", $2, $6);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($2) + strlen($6) + 6);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, "("); strcat(tmp, $2); strcat(tmp, ":F(");
                                                                                                strcat(tmp, $6);  strcat(tmp, ")"); 
                                                                                                $$ = tmp;
                                                                                                free($2);
                                                                                                free($6);
                                                                                              }
            | pg ',' gene                                                                     {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed pg:          %s, %s\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, ","); strcat(tmp, $3);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | pg ',' gene ':' 'F' '(' expr ')'                                                {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed pg:          %s, %s:F(%s)\n", $1, $3, $7);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($7) + 6);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);    strcat(tmp, ","); strcat(tmp, $3);
                                                                                                strcat(tmp, ":F("); strcat(tmp, $7);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($7);
                                                                                              }
            ;
sim         : protein '(' sim_list gene ')'                                                   { /* instantiate Kinetic Law */
                                                                                                sprintf(temp, "%c", sgn0); 
                                                                                                strcat(temp, $1);
                                                                                                strcat(temp, ";");
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($4, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($4, temp);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                                                                                num_sgn = 0;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed sim:         %s(%s%s)\n", $1, $3, $4);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($4) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, $4);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($4);
                                                                                              }
            | protein '(' sim_list gene ':' 'F' '(' expr ')' ')'                              { /* instantiate Kinetic Law */
                                                                                                sprintf(temp, "%c", sgn0); 
                                                                                                strcat(temp, $1);
                                                                                                strcat(temp, ";");
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($4, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($4, temp, $8);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;
                                                                                                num_sgn = 0;

                                                                                                if(parseInfo)
                                                                                                  printf("parsed sim:         %s(%s%s:F(%s))\n", $1, $3, $4, $8);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + strlen($4) + strlen($8) + 7);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, "(");   strcat(tmp, $3);
                                                                                                strcat(tmp, $4);  strcat(tmp, ":F("); strcat(tmp, $8);
                                                                                                strcat(tmp, "))");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                                free($4);
                                                                                                free($8);
                                                                                              }
            ;
sim_list    : sim_list gene ','                                                               { /* instantiate Kinetic Law */
                                                                                                sprintf(temp, "%c", sgn0); 
                                                                                                strcat(temp, protein);
                                                                                                strcat(temp, ";");
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($2, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($2, temp);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                                                                                num_sgn = 0;

                                                                                                if(parseInfo)
                                                                                                  printf("parsed sim_list:    %s%s,\n", $1, $2);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($2) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, $2); strcat(tmp, ",");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($2);
                                                                                              }
            | sim_list gene ':' 'F' '(' expr ')' ','                                          { /* instantiate Kinetic Law */
                                                                                                sprintf(temp, "%c", sgn0); 
                                                                                                strcat(temp, protein);
                                                                                                strcat(temp, ";");
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($2, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($2, temp, $6);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;
                                                                                                num_sgn = 0;

                                                                                                if(parseInfo)
                                                                                                  printf("parsed sim_list:    %s%s:F(%s),\n", $1, $2, $6);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($2) + strlen($6) + 6);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, $2); strcat(tmp, ":F(");
                                                                                                strcat(tmp, $6); strcat(tmp, "),");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($2);
                                                                                                free($6);
                                                                                              }
            | sgn gene ','                                                                    { /* instantiate Kinetic Law */
                                                                                                strcpy(temp, $1);
                                                                                                strcat(temp, protein);
                                                                                                strcat(temp, ";");
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($2, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = randomGeneralizedHill($2, temp);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kineticLawString);
                                                                                                  free(kineticLawString);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                rand_func = 1;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed sim_list:    %s%s,\n", $1, $2);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($2) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, $2); strcat(tmp, ",");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($2);
                                                                                              }
            | sgn gene ':' 'F' '(' expr ')' ','                                               { /* instantiate Kinetic Law */
                                                                                                strcpy(temp, $1);
                                                                                                strcat(temp, protein);
                                                                                                strcat(temp, ";");
                                                                                                
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  strcpy(xgmmlTmp, temp);
                                                                                                  xgmmlXML($2, xgmmlTmp);
                                                                                                }
                                                                                                
                                                                                                kl = KineticLaw_create();
                                                                                                react = Model_createReaction(model);
                                                                                                kLSp = explicitKineticLaw($2, temp, $6);
                                                                                                if(kLSp)
                                                                                                {
                                                                                                  KineticLaw_setFormula(kl, kLSp);
                                                                                                  free(kLSp);
                                                                                                }
                                                                                                else
                                                                                                {
                                                                                                  yyerror("NULL kineticLawString, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                Reaction_setKineticLaw(react, kl);
                                                                                                user_func = 1;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed sim_list:    %s%s:F(%s),\n", $1, $2, $6);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($2) + strlen($6) + 5);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, $2); strcat(tmp, ":F(");
                                                                                                strcat(tmp, $6); strcat(tmp, "),");
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($2);
                                                                                                free($6);
                                                                                              }
            ;
sgn         : '+'                                                                             {
                                                                                                switch(num_sgn)
                                                                                                {
                                                                                                  case 0: sgn0 = '+'; break;
                                                                                                  case 1: sgn1 = '+'; break;
                                                                                                  case 2: sgn2 = '+'; break;
                                                                                                }
                                                                                                num_sgn++;
                                                                                                num_sgn %= 3;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed sgn:         +\n");
                                                                                                  
                                                                                                tmp = (char *) malloc(2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, "+");
                                                                                                $$ = tmp;
                                                                                              }
            | '-'                                                                             {
                                                                                                switch(num_sgn)
                                                                                                {
                                                                                                  case 0: sgn0 = '-'; break;
                                                                                                  case 1: sgn1 = '-'; break;
                                                                                                  case 2: sgn2 = '-'; break;
                                                                                                }
                                                                                                num_sgn++;
                                                                                                num_sgn %= 3;
                                                                                                
                                                                                                if(parseInfo)
                                                                                                  printf("parsed sgn:         -\n");

                                                                                                tmp = (char *) malloc(2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, "-");
                                                                                                $$ = tmp;
                                                                                              }
            ;
gene        : GENE                                                                            {
                                                                                                if(xgmml)
                                                                                                {
                                                                                                  if((strlen(cytoBuf) + 64) > cytoBufSz)
                                                                                                  {
                                                                                                    tmp = realloc(cytoBuf, cytoBufSz + 256);
                                                                                                    if(!tmp)
                                                                                                    {
                                                                                                      yyerror("realloc error for cytoBuf, exiting...");
                                                                                                      YYABORT;
                                                                                                    }
                                                                                                    else
                                                                                                    {
                                                                                                      cytoBuf = tmp;
                                                                                                      cytoBufSz += 256;
                                                                                                    }
                                                                                                  }
                                                                                                  sprintf(tmpCytoBuf, "  <node id=\"%s\" label=\"%s\"/>\n", $1, $1);
                                                                                                  strcat(cytoBuf, tmpCytoBuf);
                                                                                                }

                                                                                                tot_genes++;
                                                                                                if(parseInfo)
                                                                                                  printf("parsed GENE:        %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                              }
            ;
protein     : PROTEIN                                                                         {
                                                                                                if(!Model_getSpeciesById(model, $1))
                                                                                                {
                                                                                                  species = Model_createSpecies(model);
                                                                                                  Species_setId(species, $1);
                                                                                                  Species_setName(species, $1);
                                                                                                  Species_setCompartment(species, sid);
                                                                                                  Species_setInitialConcentration(species, 1.0);
                                                                                                }
                                                                                                
                                                                                                strcpy(protein, $1);
                                                                                               
                                                                                                if(parseInfo)
                                                                                                  printf("parsed PROTEIN:     %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                              }
            ;
p_error     : PROTEIN ','                                                                     {
                                                                                                tmp = (char *) malloc(strlen($1) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, ",");
                                                                                                $$ = tmp;
                                                                                              }
            | PROTEIN ')'                                                                     {
                                                                                                tmp = (char *) malloc(strlen($1) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                              }
            ;
expr        : expr '+' expr                                                                   {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        %s+%s\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "+"); strcat(tmp, $3);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | expr '-' expr                                                                   {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        %s-%s\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "-"); strcat(tmp, $3);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | expr '*' expr                                                                   {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        %s*%s\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "*"); strcat(tmp, $3);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                                free($3);
                                                                                              }
            | expr '/' expr                                                                   {
                                                                                                if($3 == 0) {yyerror("Error: division by zero\n"); YYABORT;}
                                                                                                else
                                                                                                {
                                                                                                  if(parseInfo)
                                                                                                    printf("parsed expr:        %s/%s\n", $1, $3);
                                                                                                  
                                                                                                  tmp = (char *) malloc(strlen($1) + strlen($3) + 2);
                                                                                                  if(!tmp)
                                                                                                  {
                                                                                                    yyerror("malloc error, exiting...");
                                                                                                    YYABORT;
                                                                                                  }
                                                                                                
                                                                                                  strcpy(tmp, $1); strcat(tmp, "/"); strcat(tmp, $3);
                                                                                                  $$ = tmp;
                                                                                                  free($1);
                                                                                                  free($3);
                                                                                                }
                                                                                              }
            | '-' expr %prec UMINUS                                                           {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        (-%s)\n", $2);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($2) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, "-"); strcat(tmp, $2);
                                                                                                $$ = tmp;
                                                                                                free($2);
                                                                                              }
            | '(' expr ')'                                                                    {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        (%s)\n", $2);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($2) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, "("); strcat(tmp, $2); strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($2);
                                                                                              }
            | term                                                                            {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            ;
term        : protein                                                                         {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed term:        %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            | constant                                                                        {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed term:        %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                                free($1);
                                                                                              }
            | ABS '(' expr ')'                                                                {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        abs(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | ARCCOS '(' expr ')'                                                             {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        arccos(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | ARCSIN '(' expr ')'                                                             {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        arcsin(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | ARCTAN '(' expr ')'                                                             {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        arctan(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | CEILING '(' expr ')'                                                            {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        ceiling(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | COS '(' expr ')'                                                                {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        cos(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | EXP '(' expr ')'                                                                {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        exp(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | FLOOR '(' expr ')'                                                              {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        floor(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | LN '(' expr ')'                                                                 {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        ln(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | LOG '(' DIGITS ',' expr ')'                                                     {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        log(%s,%s)\n", $3, $5);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($3) + strlen($5) + 4);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ","); strcat(tmp, $5);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($5);
                                                                                              }
            | POWER '(' expr ',' expr ')'                                                     {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        power(%s,%s)\n", $3, $5);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($3) + strlen($5) + 4);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ","); strcat(tmp, $5);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                                free($5);
                                                                                              }
            | ROOT '(' DIGITS ',' expr ')'                                                    {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        root(%s,%s)\n", $3, $5);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($3) + strlen($5) + 4);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);  strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ","); strcat(tmp, $5);  strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($5);
                                                                                              }
            | SIN '(' expr ')'                                                                {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        sin(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            | TAN '(' expr ')'                                                                {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed expr:        tan(%s)\n", $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 3);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "("); strcat(tmp, $3);
                                                                                                strcat(tmp, ")");
                                                                                                $$ = tmp;
                                                                                                free($3);
                                                                                              }
            ;
constant    : DIGITS '.' DIGITS                                                               {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed constant:    %s,%s\n", $1, $3);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + strlen($3) + 2);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1); strcat(tmp, "."); strcat(tmp, $3);
                                                                                                $$ = tmp;
                                                                                              }
            | DIGITS                                                                          {
                                                                                                if(parseInfo)
                                                                                                  printf("parsed constant:    %s\n", $1);
                                                                                                  
                                                                                                tmp = (char *) malloc(strlen($1) + 1);
                                                                                                if(!tmp)
                                                                                                {
                                                                                                  yyerror("malloc error, exiting...");
                                                                                                  YYABORT;
                                                                                                }
                                                                                                
                                                                                                strcpy(tmp, $1);
                                                                                                $$ = tmp;
                                                                                              }
            ;
%%

int main(int argc, char **argv)
{
  int i, j=0, option;
  long seedval;

  seedval = 123456789;
  srand48(seedval);
  
  /* options parsing */
  while((option = getopt(argc, argv, "s:hkpvx")) > 0)
  {
    switch(option)
    {
      case 'h':
        printf("compile into Systems Biology Markup Language a\n");
        printf("network in the NEMO (NEtwork MOtif) language\n");
        printf("usage: nemo2sbml [options] <input file> <output file>\n");
        printf("                 -h --help\n");
        printf("                 -k print kinetic law info\n");
        printf("                 -p print parse info\n");
        printf("                 -s <seedval>, set the seed for drand48, default = 123456789\n");
        printf("                 -v print version\n");
        printf("                 -x output an XGMML file for cytoscape\n");
        return 0;
        
      case 'k':
        kineticLawInfo = 1;
        break;
        
      case 'p':
        parseInfo = 1;
        break;

      case 's':
        for(i=0; i<strlen(optarg); i++)
        {
          if(!isdigit(optarg[i]))
          {
            fprintf(stderr, "nemo2sbml: -s: \"%s\" must be an integer argument >= 0, returning...\n", optarg);
            return 1;
          }
        }
        seedval = atol(optarg);
        srand48(seedval);
        break;
      
      case 'v':
        printf("ver %s\n", VERSION);
        return 0;
        
      case 'x':
        xgmml = 1;
        
        cytoBuf = (char *) malloc(BUFSZ);
        if(!cytoBuf)
        {
          fprintf(stderr, "nemo2sbml: malloc error for cytoBuf, returning...\n");
          return 1;
        }
        
        cytoBufSz = BUFSZ;
        break;
       
      default: break;
    }
  }
  
  output[0] = 0x0;
  
  if(argv[optind] != NULL)
  {
    yyin = fopen(argv[optind], "r");
    if(yyin == NULL)
    {
      fprintf(stderr, "nemo2sbml: unable to open input file %s...\n", argv[optind]);
      return 1;
    }
    
    if(argv[optind+1] != NULL)
    {
      if(!strcmp(argv[optind], argv[optind+1]))
      {
        fprintf(stderr, "nemo2sbml: output file \"%s\" must be different than input file\n", argv[optind+1]);
        return 1;
      }
      
      yyout = fopen(argv[optind+1], "w");
      if(yyout == NULL)
      {
        fprintf(stderr, "nemo2sbml: unable to open output file %s...\n", argv[optind+1]);
        return 1;
      }
      strcpy(output, argv[optind+1]);
    }
  }


  doc = SBMLDocument_createWith(SBML_LEVEL, SBML_VERSION);
  model = SBMLDocument_createModel(doc);

  compart = Model_createCompartment(model);
  Compartment_setId(compart, sid);
  Compartment_setVolume(compart, 1.0); /* 1 femtoliter, see below */

  /* cell volume is 10^-15 liters = 1 cubic micrometer = 1 femtoliter */
  unitdef = Model_createUnitDefinition(model);
  UnitDefinition_setId(unitdef, "volume");
  UnitDefinition_setName(unitdef, "femtoliter");
  UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), 1, -15));
  
  unitdef = Model_createUnitDefinition(model);
  UnitDefinition_setId(unitdef, "microM_cell");
  UnitDefinition_setName(unitdef, "microMole/cell");
  UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("mole"), 1, -6));
  UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), -1, -15));

  unitdef = Model_createUnitDefinition(model);
  UnitDefinition_setId(unitdef, "hill_coeff");
  UnitDefinition_setName(unitdef, "Hill Coefficient");
  UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("dimensionless"), 1, 0)); /* 1 to 4 */

  unitdef = Model_createUnitDefinition(model);
  UnitDefinition_setId(unitdef, "s_fl");
  UnitDefinition_setName(unitdef, "sec/femtoliter");
  UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("second"), 1, 0));
  UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("litre"), -1, -15));

  unitdef = Model_createUnitDefinition(model);
  UnitDefinition_setId(unitdef, "s_mole");
  UnitDefinition_setName(unitdef, "sec/microMole");
  UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("second"), 1, 0));
  UnitDefinition_addUnit(unitdef, Unit_createWith(UnitKind_forName("mole"), -1, -6));

  species = Model_createSpecies(model);
  Species_setId(species, "devNull");
  Species_setName(species, "devNull");
  Species_setCompartment(species, sid);
  Species_setInitialConcentration(species, 0.0);
  Species_setBoundaryCondition(species, 1);
  Species_setConstant(species, 1);

  do
  {
    yyparse();
  }
  while(!feof(yyin));
  
  return 0;
}

int check_dor(char *dor_text)
{
  /* For the first GENE, check that at every other GENE can be reached, i.e.
   * there is a path between any two GENEs through the Ps, thus the graph is 
   * connected. Algorithm is a recursive Depth First Search
   */

  int i, j;
  char *buf, *p;
  size_t nmatch=RES;
  regmatch_t  pmatch[RES];

  /* init */
  for(i=0; i<GENES; i++)
  {
    genes[i][0] = 0x0;
    marked[i] = 0x0;
  }

  regcomp(&preg, "P[0-9]*", REG_EXTENDED);

  buf = (char *) malloc(strlen(dor_text)+8);
  if(buf == NULL)
  {
    fprintf(stderr, "check_dor: malloc error, exiting...\n");
    return 0;
  }

  /* grab all the p_lists */
  j = 0;
  for(i=0; i<strlen(dor_text); i++) /* remove spaces and tabs */
    if(dor_text[i] != ' ' && dor_text[i] != '\t')
      buf[j++] = dor_text[i];
  buf[j] = 0x0;

  i = 0;
  if(p = strtok(buf, ")"))
  {
    if(!regexec(&preg, p, nmatch, pmatch, 0)) /* make sure a P lives here */
    {
      strcpy(genes[i++], p);
      //printf("%d: %s\n", i-1, genes[i-1]);
    }
  }

  do
  {
    if(p = strtok(NULL, ")"))
    {
      if(!regexec(&preg, p, nmatch, pmatch, 0))
      {
        strcpy(genes[i++], p);
        //printf("%d: %s\n", i-1, genes[i-1]);
      }
    }
  }
  while(p);

  num_genes = i;
  
  if(num_genes < 2)
  {
    fprintf(stderr, "check_dor: error, num_genes < 2, exiting...\n");
    free(buf);
    regfree(&preg);
    return 0;
  }

  mark_neighbors(0); /* recursive depth first search */

  /* every GENE (node) should have been visited */
  for(i=0; i<num_genes; i++)
  {
    if(!marked[i]) /* failed */
    {
      fprintf(stderr, "check_dor: error, DOR graph is not connected\n");
      free(buf);
      regfree(&preg);
      return 0;
    }
  }

  free(buf);
  regfree(&preg);
  return 1;
}

/* mark neigboring nodes recursively */
void mark_neighbors(int gene_idx)
{
  int j;
  long ioff, joff;
  char p_bufi[BUFSZ], p_bufj[BUFSZ];
  size_t in=RES, jn=RES;
  regmatch_t  ip[RES], jp[RES];

  if(marked[gene_idx]) return; /* this node has been visited */
  marked[gene_idx] = 0x1;

  /* for every P of gene, mark all the other genes that P regulates */
  ioff = 0;
  while(!regexec(&preg, genes[gene_idx]+ioff, in, ip, 0))
  {
    strncpy(p_bufi, genes[gene_idx]+ioff+ip[0].rm_so, ip[0].rm_eo-ip[0].rm_so);
    p_bufi[ip[0].rm_eo-ip[0].rm_so] = 0x0;

    /* look for a match to p_bufi in other p_lists */
    for(j=0; j<num_genes; j++)
    {
      if(j==gene_idx) continue;

      joff = 0;
      while(!regexec(&preg, genes[j]+joff, jn, jp, 0))
      {
        strncpy(p_bufj, genes[j]+joff+jp[0].rm_so, jp[0].rm_eo-jp[0].rm_so);
        p_bufj[jp[0].rm_eo-jp[0].rm_so] = 0x0;

        if(!strcmp(p_bufi, p_bufj)) /* a match */
          mark_neighbors(j); /* recursive depth first search */

        joff += jp[0].rm_eo;
      }
    }
    ioff += ip[0].rm_eo;
  }
}

/* 
 Return a Kinetic Law string. This particular implementation is a generalized
 Hill Function with randomized parameters*, but can be replaced with the 
 biochemical model of your choice, as long as the code in the grammar is also
 altered accordingly to pass the proper parameters.
 tfs (transcription factors string) format: ([+-]P[,;])*
 
 *see Likhoshvai V., Ratushny A., "Generalized Hill Function Method for 
 Modeling Molecular Processes", Journal of Bioinformatics and Computaional
 Biology, vol 5, issue 2B, pg 521-531, April 2007
 The implementation here is slightly modified from that in the article.
 Comment out the "#define NON_LINEAR" statement above to use linear terms
 only in the Hill function
*/
char * randomGeneralizedHill(char *geneRegulated, char *tfs)
{
  int denom_sz, numer_sz;
  char buf[32], *denom, *numer, *p=0x0, savP[32], *sav_tfs;
  KineticLaw_t  *dl;
  Reaction_t *degrad;

  denom = (char *) malloc(NLT_SZ);
  numer = (char *) malloc(NLT_SZ);
  if(denom==NULL || numer==NULL)
  {
    fprintf(stderr, "randomGeneralizedHill: malloc error, returning NULL...\n");
    return NULL;
  }
  denom_sz = numer_sz = NLT_SZ;

  /* degradation */
  dl = KineticLaw_create();
  degrad = Model_createReaction(model);
  strcpy(denom, "P"); /* borrow the denom array as a buffer */
  strcat(denom, strstr(geneRegulated, "G")+1);
  strcat(denom, "_degrad");
  Reaction_setId(degrad, denom);
  strcpy(denom, "P"); 
  strcat(denom, strstr(geneRegulated, "G")+1);
  strcat(denom, " degradation");
  Reaction_setName(degrad, denom);
  
  sprintf(buf, "dc_%d", parameterIndex);
  dparam = Parameter_createWith(buf, 0.01+drand48()/10, "dimensionless"); /* drand48 is uniform rand [0.0 - 1.0) */
  KineticLaw_addParameter(dl, dparam);
  Reaction_setReversible(degrad, 0); /* balance of degradation is below */
  
  /* synthesis */
  strcpy(denom, "P"); 
  strcat(denom, strstr(geneRegulated, "G")+1);
  strcat(denom, "_synthesis");
  Reaction_setId(react, denom);

  strcpy(denom, "P"); 
  strcat(denom, strstr(geneRegulated, "G")+1);
  strcat(denom, " synthesis");
  Reaction_setName(react, denom);
  Reaction_setReversible(react, 0);

  strcpy(denom, "P"); 
  strcat(denom, strstr(geneRegulated, "G")+1);
  reactant = SpeciesReference_createWith(denom, 1.0, 1);
  Reaction_addProduct(react, reactant);
  Reaction_addReactant(react, SpeciesReference_createWith("devNull", 1.0, 1));
  
  /* balance of degradation */
  Reaction_addProduct(degrad, SpeciesReference_createWith("devNull", 1.0, 1));
  Reaction_addReactant(degrad, reactant);
  strcpy(numer, buf); /* borrow the numer array as a buffer */
  strcat(numer, "*");
  strcat(numer, denom);
  KineticLaw_setFormula(dl, numer);
  Reaction_setKineticLaw(degrad, dl);
  
  
  if(!Model_getSpeciesById(model, denom)) /* has this species been created yet? */
  {
    species = Model_createSpecies(model);
    Species_setId(species, denom);
    Species_setName(species, denom);
    Species_setCompartment(species, sid);
    Species_setInitialConcentration(species, 1.0);
  }
  
  /* save tfs string */
  sav_tfs = (char *) malloc(strlen(tfs)+8);
  if(sav_tfs == NULL)
  {
    fprintf(stderr, "randomGeneralizedHill: malloc error, returning NULL Kinetic Law for %s\n", geneRegulated);
    return NULL;
  }
  strcpy(sav_tfs, tfs);

  /* build the numerator and denominator strings */
  p = strtok(tfs, " ,;)");
  if(p)
  {
    strcpy (numer, "(");
    strcpy (denom, "/(1+power(");

    sprintf(buf, "B_%d", parameterIndex);
    strcat (numer, buf);

    param = Parameter_createWith(buf, 0.0001+drand48(), "microM_cell"); /* 0.0001 ~ 1.0 */
    KineticLaw_addParameter(kl, param);

    if(strstr(p, "+")) /* activator */
    {
      strcat (numer, "*power(");
      sprintf(buf, "%s", strstr(p, "P"));
      strcpy(savP, "+");
      strcat(savP, buf);
      msr = ModifierSpeciesReference_createWith(buf);
      Reaction_addModifier(react, msr);
      strcat (numer, buf);
      strcat (denom, buf);
      strcat (numer, "/");
      strcat (denom, "/");
      sprintf(buf, "K_%d", parameterIndex);
      strcat (numer, buf);
      strcat (denom, buf);
      strcat (numer, ", ");
      strcat (denom, ", ");

      param = Parameter_createWith(buf, 0.5+drand48(), "microM_cell"); /* 0.5 ~ 1.5 */
      KineticLaw_addParameter(kl, param);

      sprintf(buf, "n_%d", parameterIndex);
      strcat (numer, buf);
      strcat (denom, buf);

      param = Parameter_createWith(buf, 1.0+floor(4*drand48()), "hill_coeff"); /* 1 - 4 */
      KineticLaw_addParameter(kl, param);

      strcat (numer, ")");
      strcat (denom, ")");
      
      if((denom_sz - strlen(denom)) < NLT_SZ)
      {
        denom = realloc(denom, denom_sz + NLT_SZ);
        if(denom == NULL)
        {
          fprintf(stderr, "randomGeneralizedHill: realloc error, returning NULL Kinetic Law for %s\n", geneRegulated);
          return NULL;
        }
        denom_sz += NLT_SZ;
      }
      if((numer_sz - strlen(numer)) < NLT_SZ)
      {
        numer = realloc(numer, numer_sz + NLT_SZ);
        if(numer == NULL)
        {
          fprintf(stderr, "randomGeneralizedHill: realloc error, returning NULL Kinetic Law for %s\n", geneRegulated);
          return NULL;
        }
        numer_sz += NLT_SZ;
      }
#ifdef NON_LINEAR
      strcat (numer, insertNonLinearTerms(savP, sav_tfs));
      strcat (denom, insertNonLinearTerms(savP, sav_tfs));
#endif
    }
    else               /* repressor */
    {
      sprintf(buf, "%s", strstr(p, "P"));
      strcpy(savP, "-");
      strcat(savP, buf);
      msr = ModifierSpeciesReference_createWith(buf);
      Reaction_addModifier(react, msr);
      strcat (denom, buf);
      strcat (denom, "/");
      sprintf(buf, "K_%d", parameterIndex);
      strcat (denom, buf);
      strcat (denom, ", ");

      param = Parameter_createWith(buf, 0.5+drand48(), "microM_cell");
      KineticLaw_addParameter(kl, param);

      sprintf(buf, "n_%d", parameterIndex);
      strcat (denom, buf);

      param = Parameter_createWith(buf, 1.0+floor(4*drand48()), "hill_coeff");
      KineticLaw_addParameter(kl, param);
      strcat (denom, ")");
      
      if((denom_sz - strlen(denom)) < NLT_SZ)
      {
        denom = realloc(denom, denom_sz + NLT_SZ);
        if(denom == NULL)
        {
          fprintf(stderr, "randomGeneralizedHill: realloc error, returning NULL Kinetic Law for %s\n", geneRegulated);
          return NULL;
        }
        denom_sz += NLT_SZ;
      }
#ifdef NON_LINEAR
      strcat (denom, insertNonLinearTerms(savP, sav_tfs));
#endif
    }
    
    parameterIndex++;
  }
  else
  {
    fprintf(stderr, "randomGeneralizedHill: NULL tfs, returning NULL Kinetic Law for %s\n", geneRegulated);
    return NULL;
  }

  do
  {
    p = strtok(NULL, " ,;)");
    if(p)
    {
      strcat (numer, "+");
      strcat (denom, "+power(");

      sprintf(buf, "B_%d", parameterIndex);
      strcat (numer, buf);

      param = Parameter_createWith(buf, 0.0001+drand48(), "microM_cell");
      KineticLaw_addParameter(kl, param);

      if(strstr(p, "+")) /* activator */
      {
        strcat (numer, "*power(");
        sprintf(buf, "%s", strstr(p, "P"));
        strcpy(savP, "+");
        strcat(savP, buf);
        msr = ModifierSpeciesReference_createWith(buf);
        Reaction_addModifier(react, msr);
        strcat (numer, buf);
        strcat (denom, buf);
        strcat (numer, "/");
        strcat (denom, "/");
        sprintf(buf, "K_%d", parameterIndex);
        strcat (numer, buf);
        strcat (denom, buf);
        strcat (numer, ", ");
        strcat (denom, ", ");

        param = Parameter_createWith(buf, 0.5+drand48(), "microM_cell");
        KineticLaw_addParameter(kl, param);

        sprintf(buf, "n_%d", parameterIndex);
        strcat (numer, buf);
        strcat (denom, buf);

        param = Parameter_createWith(buf, 1.0+floor(4*drand48()), "hill_coeff");
        KineticLaw_addParameter(kl, param);

        strcat (numer, ")");
        strcat (denom, ")");
        
        if((denom_sz - strlen(denom)) < NLT_SZ)
        {
          denom = realloc(denom, denom_sz + NLT_SZ);
          if(denom == NULL)
          {
            fprintf(stderr, "randomGeneralizedHill: realloc error, returning NULL Kinetic Law for %s\n", geneRegulated);
            return NULL;
          }
          denom_sz += NLT_SZ;
        }
        if((numer_sz - strlen(numer)) < NLT_SZ)
        {
          numer = realloc(numer, numer_sz + NLT_SZ);
          if(numer == NULL)
          {
            fprintf(stderr, "randomGeneralizedHill: realloc error, returning NULL Kinetic Law for %s\n", geneRegulated);
            return NULL;
          }
          numer_sz += NLT_SZ;
        }
#ifdef NON_LINEAR
        strcat (numer, insertNonLinearTerms(savP, sav_tfs));
        strcat (denom, insertNonLinearTerms(savP, sav_tfs));
#endif
      }
      else               /* repressor */
      {
        sprintf(buf, "%s", strstr(p, "P"));
        strcpy(savP, "-");
        strcat(savP, buf);
        msr = ModifierSpeciesReference_createWith(buf);
        Reaction_addModifier(react, msr);
        strcat (denom, buf);
        strcat (denom, "/");
        sprintf(buf, "K_%d", parameterIndex);
        strcat (denom, buf);
        strcat (denom, ", ");

        param = Parameter_createWith(buf, 0.5+drand48(), "microM_cell");
        KineticLaw_addParameter(kl, param);

        sprintf(buf, "n_%d", parameterIndex);
        strcat (denom, buf);

        param = Parameter_createWith(buf, 1.0+floor(4*drand48()), "hill_coeff");
        KineticLaw_addParameter(kl, param);
        strcat (denom, ")");
        
        if((denom_sz - strlen(denom)) < NLT_SZ)
        {
          denom = realloc(denom, denom_sz + NLT_SZ);
          if(denom == NULL)
          {
            fprintf(stderr, "randomGeneralizedHill: realloc error, returning NULL Kinetic Law for %s\n", geneRegulated);
            return NULL;
          }
          denom_sz += NLT_SZ;
        }
#ifdef NON_LINEAR
        strcat (denom, insertNonLinearTerms(savP, sav_tfs));
#endif
      }
      
      parameterIndex++;
    }
  }
  while(p);

  kineticLawString = (char *) malloc(strlen(numer)+strlen(denom)+8);
  if(kineticLawString==NULL)
  {
    fprintf(stderr, "randomGeneralizedHill: malloc error, returning NULL Kinetic Law for %s\n", geneRegulated);
    return NULL;
  }
  strcpy(kineticLawString, numer);
  strcat(kineticLawString, ")");
  strcat(kineticLawString, denom);
  strcat(kineticLawString, ")");
  if(kineticLawInfo)
    printf("Kinetic Law for %s = %s\n", geneRegulated, kineticLawString);

  free(denom);
  free(numer);
  free(sav_tfs);
  return kineticLawString;
}

char * insertNonLinearTerms(char *savP, char *sav_tfs)
{
  char buf[32], *pgp, *pt, *tfs;
  
  tfs = (char *) malloc(strlen(sav_tfs)+8);
  if(tfs == NULL)
  {
    fprintf(stderr, "insertNonLinearTerms: malloc error, returning NULL...\n");
    return NULL;
  }
  strcpy(tfs, sav_tfs);

  returnString[0] = 0x0;

  pgp = tfs;
  do
  {
    pt = strsep(&pgp, " ,;)");
    if(pt && strlen(pt) > 0)
    {
      if(strcmp(pt, savP))  /* different */
      {
        strcat(returnString, "*power(");
        sprintf(buf, "%s", strstr(pt, "P"));
        strcat(returnString, buf);
        strcat(returnString, "/");
        sprintf(buf, "K_%d", ++parameterIndex);
        strcat(returnString, buf);

        //param = Parameter_createWith(buf, 0.5*drand48(), "microM_cell"); /* 0.5 ~ 1.5 */
        param = Parameter_createWith(buf, 1+1000*drand48(), "microM_cell"); /* 1 ~ 1001 */
        KineticLaw_addParameter(kl, param);
        
        strcat(returnString, ", ");
        sprintf(buf, "n_%d", parameterIndex);
        strcat(returnString, buf);
        strcat(returnString, ")");

        //param = Parameter_createWith(buf, 1.0+floor(4*drand48()), "hill_coeff"); /* 1 - 4 */
        param = Parameter_createWith(buf, 0.2+4*drand48(), "hill_coeff"); /* 0.2 - 4.2 */
        KineticLaw_addParameter(kl, param);
      }
    }
  }
  while(pt && strlen(pt) > 0);
  
  free(tfs);
  return returnString;
}

char * explicitKineticLaw(char *geneRegulated, char *tfs, char *explicitFunction)
{
  int i, j;
  char a1[64], a2[64], buf[32], *p=0x0, proteins[64][32];
  KineticLaw_t  *dl;
  Reaction_t *degrad;

  /* degradation */
  dl = KineticLaw_create();
  degrad = Model_createReaction(model);
  strcpy(a1, "P");
  strcat(a1, strstr(geneRegulated, "G")+1);
  strcat(a1, "_degrad");
  Reaction_setId(degrad, a1);
  strcpy(a1, "P"); 
  strcat(a1, strstr(geneRegulated, "G")+1);
  strcat(a1, " degradation");
  Reaction_setName(degrad, a1);
  Reaction_setReversible(degrad, 0); /* balance of degradation is below */
  
  /* synthesis */
  strcpy(a1, "P"); 
  strcat(a1, strstr(geneRegulated, "G")+1);
  strcat(a1, "_synthesis");
  Reaction_setId(react, a1);

  strcpy(a1, "P"); 
  strcat(a1, strstr(geneRegulated, "G")+1);
  strcat(a1, " synthesis");
  Reaction_setName(react, a1);
  Reaction_setReversible(react, 0);

  strcpy(a1, "P"); 
  strcat(a1, strstr(geneRegulated, "G")+1);
  reactant = SpeciesReference_createWith(a1, 1.0, 1);
  Reaction_addProduct(react, reactant);
  Reaction_addReactant(react, SpeciesReference_createWith("devNull", 1.0, 1));
  
  /* balance of degradation */
  Reaction_addProduct(degrad, SpeciesReference_createWith("devNull", 1.0, 1));
  Reaction_addReactant(degrad, reactant);
  strcpy(a2, a1);
  KineticLaw_setFormula(dl, a2);
  Reaction_setKineticLaw(degrad, dl);
  
  if(!Model_getSpeciesById(model, a1)) /* has this species been created yet? */
  {
    species = Model_createSpecies(model);
    Species_setId(species, a1);
    Species_setName(species, a1);
    Species_setCompartment(species, sid);
    Species_setInitialConcentration(species, 1.0);
  }

  /* build a list of proteins (<=64) used in explicitFunction */
  p = explicitFunction;
  for(i=0; i<64; i++)
  {
    j = 0; proteins[i][0] = 0x0;
    if(p && (p=strstr(p, "P")))
    {
      j = 1;
      while(isdigit(p[j]))
      {
        proteins[i][j] = p[j++];
        proteins[i][j] = 0x0;
      }
      if(j>1) proteins[i][0] = 'P';
      p++;
    }
  }

  /* add reaction modifiers, and make sure that each modifier appears in explicitFunction */
  p = strtok(tfs, " ,;)");
  if(p)
  {
    sprintf(buf, "%s", strstr(p, "P"));
    /* is this modifier used in explicitFunction? */
    j = 0;
    for(i=0; i<64; i++)
      if(!strcmp(proteins[i], buf)) {j = 1; break;}

    if(!j)
    {
      fprintf(stderr, "explicitKineticLaw: protein %s unused in %s, returning NULL Kinetic Law for %s\n", buf, explicitFunction, geneRegulated);
      return NULL;
    }
    
    msr = ModifierSpeciesReference_createWith(buf);
    Reaction_addModifier(react, msr);
  }
  else
  {
    fprintf(stderr, "explicitKineticLaw: NULL tfs, returning NULL Kinetic Law for %s\n", geneRegulated);
    return NULL;
  }

  do
  {
    p = strtok(NULL, " ,;)");
    if(p)
    {
      sprintf(buf, "%s", strstr(p, "P"));
      /* is this modifier used in explicitFunction? */
      j = 0;
      for(i=0; i<64; i++)
        if(!strcmp(proteins[i], buf)) {j = 1; break;}
    
      if(!j)
      {
        fprintf(stderr, "explicitKineticLaw: protein %s unused in %s, returning NULL Kinetic Law for %s\n", buf, explicitFunction, geneRegulated);
        return NULL;
      }
      
      msr = ModifierSpeciesReference_createWith(buf);
      Reaction_addModifier(react, msr);
    }
  }
  while(p);
  
  if(kineticLawInfo)
    printf("Kinetic Law for %s = %s\n", geneRegulated, explicitFunction);
  
  return explicitFunction;
}

/* 
 write edge information into cytoBuf
 tfs (transcription factors string) format: ([+-]P[,;])* 
*/
void xgmmlXML(char *geneRegulated, char *tfs)
{
  int act;
  static int edgeId=1;

  p = strtok(tfs, " ,;)");
  if(p)
  {
    if((strlen(cytoBuf) + 1024) > cytoBufSz)
    {
      tmp = realloc(cytoBuf, cytoBufSz + 1024);
      if(!tmp)
      {
        fprintf(stderr, "nemo2sbml: realloc error for cytoBuf, unable to generate xgmml...");
        xgmml = 0;
      }
      else
      {
        cytoBuf = tmp;
        cytoBufSz += 1024;
      }
    }

    if(strstr(p, "+")) /* activator */
    {
      act = 1;
      sprintf(tmpCytoBuf, "  <edge id=\"%d\" source=\"G%s\" target=\"%s\" label=\"activation\">\n", 
              edgeId++, strstr(p, "P")+1, geneRegulated);
    }
    else               /* repressor */
    {
      act = 0;
      sprintf(tmpCytoBuf, "  <edge id=\"%d\" source=\"G%s\" target=\"%s\" label=\"repression\">\n", 
              edgeId++, strstr(p, "P")+1, geneRegulated);
    }
              
    strcat(cytoBuf, tmpCytoBuf);
    strcat(cytoBuf, "    <graphics>\n");
    strcat(cytoBuf, "      <att>\n");
    strcat(cytoBuf, "        <att name=\"sourceArrow\" value=\"0\"/>\n");
    
    if(act)
      strcat(cytoBuf, "        <att name=\"targetArrow\" value=\"3\"/>\n");
    else
      strcat(cytoBuf, "        <att name=\"targetArrow\" value=\"15\"/>\n");
      
    strcat(cytoBuf, "      </att>\n");
    strcat(cytoBuf, "    </graphics>\n");
    strcat(cytoBuf, "  </edge>\n");
  }
  
  do
  {
    p = strtok(NULL, " ,;)");
    if(p)
    {
      if((strlen(cytoBuf) + 1024) > cytoBufSz)
      {
        tmp = realloc(cytoBuf, cytoBufSz + 1024);
        if(!tmp)
        {
          fprintf(stderr, "nemo2sbml: realloc error for cytoBuf, unable to generate xgmml...");
          xgmml = 0;
        }
        else
        {
          cytoBuf = tmp;
          cytoBufSz += 1024;
        }
      }

      if(strstr(p, "+")) /* activator */
      {
        act = 1;
        sprintf(tmpCytoBuf, "  <edge id=\"%d\" source=\"G%s\" target=\"%s\" label=\"activation\">\n", 
                edgeId++, strstr(p, "P")+1, geneRegulated);
      }
      else               /* repressor */
      {
        act = 0;
        sprintf(tmpCytoBuf, "  <edge id=\"%d\" source=\"G%s\" target=\"%s\" label=\"repression\">\n", 
                edgeId++, strstr(p, "P")+1, geneRegulated);
      }
              
      strcat(cytoBuf, tmpCytoBuf);
      strcat(cytoBuf, "    <graphics>\n");
      strcat(cytoBuf, "      <att>\n");
      strcat(cytoBuf, "        <att name=\"sourceArrow\" value=\"0\"/>\n");
    
      if(act)
        strcat(cytoBuf, "        <att name=\"targetArrow\" value=\"3\"/>\n");
      else
        strcat(cytoBuf, "        <att name=\"targetArrow\" value=\"15\"/>\n");
      
      strcat(cytoBuf, "      </att>\n");
      strcat(cytoBuf, "    </graphics>\n");
      strcat(cytoBuf, "  </edge>\n");
    }
  }
  while(p);
}
