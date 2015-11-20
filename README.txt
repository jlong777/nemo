Changelog:

2015-05-25      James Long <jlong@alaska.edu>
        * nemo.y: removed the difference equation term from randomGeneralizedHill(),
        and removed references to 'tau'.
        Blunder detected: The Hill Function is *already* a rate equation, so there 
        is no need for a difference equation at all. What versions 1.5 and 1.6 are 
        doing is speeding up degradation by adding 'one' to the degradation 
        coefficient, and scaling the result by '1/tau'. 
        Also removed the difference equation for explicitKineticLaw(), be sure that 
        this is a rate equation if you explicitly specify one.

2012-05-24      James Long <jlong@alaska.edu>

        * nemo.y: removed the difference equation term from the degradation 
        reaction, moving it into a new reaction called "differenceEquationTerm".
        The explicitKineticLaw() routine is also now implemented as a difference 
        equation as above, was not done in ver 1.5, was unchanged from ver 1.4.
        Added a new option, "-t <double tau>" to set the time interval 'tau' 
        for the difference equations, default = 1.0 sec/femtoliter for 
        randomGeneralizedHill(), and sec/microMole for explicitKineticLaw().
        Added a check to make sure output file is different than input file.
        Changed "-h" option to output the default seed.


2012-05-22      James Long <jlong@alaska.edu>

        * nemo.y: modified the Hill Functions to represent *rate* 
        equations, expressed as a difference equation per unit time.
        This change produces proper behavior in the COPASI biochemical 
        simulator.
        
        * libsbml-2.3.5.tar.gz: modified source as documented below,
        and documented in libsbml's README.txt, modified tar ball is
        renamed libsbml_jlong-2.3.5.tar.gz


RANGE (RAndom Network GEnerator) generates a random network in the
NEMO (NEtwork MOtif) language. NEMO is recognized by a yacc grammar,
which outputs SBML that can be used as input to a biochemical 
simulator such as COPASI. See http://range.sf.net for examples of how
to use NEMO alone to create an explicit SBML model of your network.

In addition, a "-x" flag to nemo2sbml produces an XGMML (eXtensible Graph
Markup and Modeling Language) file for the network, allowing it to be viewed
in cytoscape (http://www.cytoscape.org/).

Code:

range.c     - source code to output a random network in the NEMO language
range-0.8.c - earlier version that makes different networks than current version
nemo.lex    - parser for the NEMO yacc grammar
nemo.y      - yacc file for NEMO
add_noise.r - R code to add noise to COPASI biochemical simulator output

INSTALL
=======
libsbml_jlong-2.3.5 included with this release, which is libsbml-2.3.5 with the 
following modifications applied to enable compatibility with gcc-4.4:

#include <cstdio>
#include <cstring>

added to Compartment.cpp, Event.cpp, FunctionDefinition.cpp, Parameter.cpp, 
Reaction.cpp, SimpleSpeciesReference.cpp, Species.cpp, and UnitDefinition.cpp

BinInputStream* SBMLSchemaInputSource::makeStream () const;

changed to

BinInputStream* makeStream () const;

in SBMLSchemaInputSource.h, line 79


yacc and lex also need to be installed on your system (or flex and bison, but 
the #defines for GENES and STR_SZ in nemo.y need to be reduced, otherwise the 
semantic stack in bison is too large. This reduction means that nemo2sbml can
handle networks with a max node size of about 3000).

1) gcc -o range range.c -lm
2) yacc -d nemo.y (or bison -y -d nemo.y)
3) lex nemo.lex   (or flex nemo.lex)
4) gcc -o nemo2sbml lex.yy.c y.tab.c -ll -lm -lsbml (may need -ly for yacc)
or gcc -o nemo2sbml lex.yy.c y.tab.c -lfl -lm -lsbml for flex/bison

usage: ./range <number of nodes in network (>=100, <=16,000)> | ./nemo2sbml
or if you have a file in the NEMO language do
       cat <file> | ./nemo2sbml
or
       ./nemo2sbml <file> [output_file_prefix]

where "output_file_prefix" will be used for the model id, and have ".xml"
appended for the output file name. "output_file_prefix" is one word. 

The above command with "output_file_prefix" omitted will produce a default
file "regulatoryNetwork_<number-of-desired-nodes-in-network>genes_0.xml"
which for 500 genes is regulatoryNetwork_500genes_0.xml

The xml file is then input to a biochemical simulator, such as COPASI.
For output exported from COPASI, use add_noise.r to simulate noisy data, 
see add_noise.r for details.

A -h to either range or nemo2sbml will list other options.


