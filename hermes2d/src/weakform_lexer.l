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
#include "weakform_lexer.h"

    
// modification of the declaration of yylex()
#define YY_DECL Token* wf_get_token()
  
  
const char* input_text;
  
// modification of the scanner input 
#define YY_INPUT(buf,result,max_size) \
  { \
    int c = *input_text++; \
    result = !c ? YY_NULL : (buf[0] = c, 1); \
  }

  
static Token simple_tokens[] =
{ 
  { T_ERROR, 0, "???",  NULL }, //  0
  { T_EOF,   0, "eof",  NULL }, //  1
  { T_VOL,   0, "vol",  NULL }, //  2
  { T_SURF,  0, "surf", NULL }, //  3
  { T_COMMA, 0, ",",    NULL }, //  4
  { T_COLON, 0, ":",    NULL }, //  5
  { T_X,     0, "x",    NULL }, //  6
  { T_Y,     0, "y",    NULL }, //  7
  { T_XX,    0, "xx",   NULL }, //  8
  { T_YY,    0, "yy",   NULL }, //  9
  { T_XY,    0, "xy",   NULL }, // 10
  { T_PLUS,  0, "+",    NULL }, // 11
  { T_MINUS, 0, "-",    NULL }, // 12
  { T_STAR,  0, "*",    NULL }, // 13
  { T_SLASH, 0, "/",    NULL }, // 14
  { T_BRA,   0, "(",    NULL }, // 15
  { T_KET,   0, ")",    NULL }, // 16
  { T_POWER, 0, "^",    NULL }, // 17
  { T_UNDER, 0, "_",    NULL }, // 18
  { T_IMAG,  0, "i",    NULL }  // 19
};


static char* strndup(char* str, int n)
{
  char* result;
  if (n == -1) n = strlen(str);
  result = new char[n+1];
  if (result == NULL) exit(1); // fixme
  strncpy(result, str, n);
  result[n] = 0;
  return result;
}


#define TABLE_SIZE 128 // must be a power of two
static Token* symbol_table[TABLE_SIZE];


static unsigned hash(char* str)
{
  unsigned hash = 0;
  while (*str)
    hash = (*str++ & 0x1f) ^ (hash << 4);
  return hash;
}


static Token* new_token(int type)
{
  Token* t = new Token;
  if (t == NULL) exit(1); // fixme
  memset(t, 0, sizeof(Token));
  t->type = type;
  return t;
}


static Token* find_symbol(char* name)
{
  unsigned index = hash(name) & (TABLE_SIZE-1);
  Token* t = symbol_table[index];
  while (t != NULL && strcmp(t->lexeme, name))
    t = t->next;
  if (t != NULL) return t;
  t = new_token(T_IDENT);
  t->lexeme = strndup(name, -1);
  t->next = symbol_table[index];
  return symbol_table[index] = t;
}


%}
%%
%{
  Token* t;
%}
  
<<EOF>>     { return simple_tokens+1;  }
[ \t\r\n]   { }
"vol"       { return simple_tokens+2;  }
"surf"      { return simple_tokens+3;  }
,           { return simple_tokens+4;  }
:           { return simple_tokens+5;  }
"x"         { return simple_tokens+6;  }
"y"         { return simple_tokens+7;  }
"xx"        { return simple_tokens+8;  }
"yy"        { return simple_tokens+9;  }
"xy"        { return simple_tokens+10; }
\+          { return simple_tokens+11; }
-           { return simple_tokens+12; }
\*          { return simple_tokens+13; }
\/          { return simple_tokens+14; }
\(          { return simple_tokens+15; }
\)          { return simple_tokens+16; }
\^          { return simple_tokens+17; }
\_          { return simple_tokens+18; }
"i"         { return simple_tokens+19; }

[a-zA-Z][a-zA-Z0-9]* { return find_symbol(yytext); }

[1-9][0-9]*(\.[0-9]*)?((e|E)(\+|-)?[0-9]+)?   { t = new_token(T_NUMBER);
                                                sscanf(yytext, "%le", &(t->value));
                                                return t; }

.           { return simple_tokens+0;  }

%%


void wf_lexer_init(const char* input)
{
  memset(symbol_table, 0, TABLE_SIZE * sizeof(Token*));
  input_text = input;
}


void wf_lexer_free()
{
  int i;
  Token *p1, *p2;
  for (i = 0; i < TABLE_SIZE; i++)
  {
    p1 = symbol_table[i];
    while (p1 != NULL)
    {
      p2 = p1;
      p1 = p1->next;
      free(p2);
    }
  }
}


#if 0 // test code 
int main()
{
  Token* t;
  wf_lexer_init("vol u,v:\" (u_x*v_x + u_y*v_y)/Re + u*v/tau + 15.6e-2*i + (xprev*u_x + yprev*u_y)*v - u^2.5E234");
  while ((t = wf_get_token())->type != T_EOF)
    printf("{ %2d, %g, \"%s\" }\n", t->type, t->value, (t->lexeme == NULL) ? "" : t->lexeme);
  wf_lexer_free();
  return 0;
}
#endif


int yywrap() { return 1; }
