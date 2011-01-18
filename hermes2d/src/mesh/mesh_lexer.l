%{
// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

// $Id$

#include <string.h>
#include <stdlib.h>
#include "mesh_lexer.h"

// modification of the declaration of yylex()
#define YY_DECL MToken* mesh_get_token()


static MToken simple_tokens[] =
{ 
  { MT_ERROR,   0, "???" }, //  0
  { MT_EOF,     0, "eof" }, //  1
  { MT_EQUAL,   0, "="   }, //  2
  { MT_SEMICOL, 0, ";"   }, //  3
  { MT_COMMA,   0, ","   }, //  4
  { MT_PLUS,    0, "+"   }, //  5
  { MT_MINUS,   0, "-"   }, //  6
  { MT_STAR,    0, "*"   }, //  7
  { MT_SLASH,   0, "/"   }, //  8
  { MT_BRA,     0, "("   }, //  9
  { MT_KET,     0, ")"   }, // 10
  { MT_POWER,   0, "^"   }, // 11
  { MT_BEGIN,   0, "{"   }, // 12
  { MT_END,     0, "}"   }  // 13
};

static MToken token;
int mesh_lexer_line_num;

extern void serror(const char* msg, ...);

%}
%%
  
<<EOF>>     { return simple_tokens+1;  }
\#.*$       { }
\n          { mesh_lexer_line_num++; }
[ \t\r]     { }
=           { return simple_tokens+2;  }
;           { return simple_tokens+3;  }
,           { return simple_tokens+4;  }
\+          { return simple_tokens+5;  }
-           { return simple_tokens+6;  }
\*          { return simple_tokens+7;  }
\/          { return simple_tokens+8;  }
\(          { return simple_tokens+9;  }
\)          { return simple_tokens+10; }
\^          { return simple_tokens+11; }
\{          { return simple_tokens+12; }
\}          { return simple_tokens+13; }

[a-zA-Z][a-zA-Z0-9_]* { token.type = MT_IDENT;
                        token.text = yytext;
                        token.value = 0.0;
                        return &token; }

[0-9]+(\.[0-9]*)?((e|E)(\+|-)?[0-9]+)?  { token.type = MT_NUMBER;
                                          sscanf(yytext, "%le", &(token.value));
                                          token.text = NULL;
                                          return &token; }

.           { serror("invalid character."); return NULL;  }

%%


void mesh_lexer_init(FILE* f)
{
  yyrestart(f);
  mesh_lexer_line_num = 1;
}


#if 0  // test code 
int main()
{
  MToken* t;
  FILE* f = fopen("../util/meshview/bracket.mesh", "r");
  mesh_lexer_init(f);
  while ((t = mesh_get_token())->type != MT_EOF)
    printf("{ %2d, %g, \"%s\" }\n", t->type, t->value, (t->text == NULL) ? "" : t->text);
  return 0;
}
#endif


int yywrap() { return 1; }
