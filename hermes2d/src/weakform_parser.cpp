// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
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

#include "h2d_common.h"
#include "weakform.h"
#include "weakform_lexer.h"
//#include "weakform_parser.h"

/*

  form    := type [ident "," ident ":"] expr
  type    := "vol" | "surf"
  expr    := term | term "+" expr | term "-" expr
  term    := power | power "*" term | power "/" term
  power   := factor | factor "^" expon
  expon   := number | "(" expr ")"
  factor  := number | ident | ident_partial | spvar | "(" expr ")"
  partial := "x" | "y" | "xx" | "yy" | "xy"
  spvar   := "x" | "y"

*/


static Token * token;


static void next_token(void)
{
  token = wf_get_token();
}


static volatile void error(char* message)
{
  printf("error: '%s' -- %s\n", token->lexeme, message);
  exit(1);
}


static void check_for(int type, char* message)
{
  if (token->type != type) error(message);
  next_token();
}


static Node* make_tree(Token* token, Node* left, Node* right, Node* cond)
{
  Node* n = (Node*) malloc(sizeof(Node));
  if (n == NULL) exit(1);
  n->token = token;
  n->left = left;
  n->right = right;
  n->cond = cond;
  return n;
}


NODE* factor(void)
{
  TOKEN* t = token;
  if (t->type != T_IDENT && t->type != T_INTEGER && t->type != T_STRING)
    error("expecting identifier, integer or string");
  next_token();
  return make_tree(t, NULL, NULL, NULL);
}


NODE* term(void)
{
  TOKEN* t;
  NODE *l, *r;
  l = factor();
  if (token->type == T_STAR || token->type == T_SLASH) {
    t = token;
    next_token();
    r = factor();
    return make_tree(t, l, r, NULL);
  }
  return l;
}


NODE* expression(void)
{
  TOKEN* t;
  NODE *l, *r;
  l = term();
  if (token->type == T_PLUS || token->type == T_MINUS) {
    t = token;
    next_token();
    r = term();
    return make_tree(t, l, r, NULL);
  }
  return l;
}


NODE* ifstatement(void)
{
  TOKEN* t = token;
  NODE *l, *r, *c;
  check_for(T_IF, "'if' expected");
  check_for(T_BRA, "'(' expected");
  c = expression();
  check_for(T_KET, "')' expected");
  l = goal();
  r = (token->type == T_ELSE) ? (next_token(), goal()) : NULL;
  return make_tree(t, l, r, c);
}


NODE* goal(void)
{
  NODE* n = (token->type == T_IF) ? ifstatement() : expression();
  check_for(T_SEMICOL, "semicolon expected");
  return n;
}


void wf_parser_init(const char* input)
{
  wf_lexer_init(input);
  token = get_token();
}


void wf_parser_free(void)
{
  wf_lexer_free();
}
