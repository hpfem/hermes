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

#include "h2d_common.h"
#include "mesh_lexer.h"
#include "mesh_parser.h"
#include <sstream>


////////////////////////////////////////////////////////////////////////////////////////////////////

static MToken* token; // current token from the lexer
static const char* H2D_PARSER_FILENAME;

inline void next_token()
  { token = mesh_get_token(); }

#define follows(what) (token->type == (what))


void serror(const char* msg, ...)
{
  char buffer[1000];
  va_list ap;
  va_start(ap, msg);
  vsprintf(buffer, msg, ap);
  va_end(ap);
  printf("ERROR: Mesh::load(const char*):\n"
         "%s, line %d: %s\n", H2D_PARSER_FILENAME, mesh_lexer_line_num, buffer);
  exit(1);
}


//// memory management /////////////////////////////////////////////////////////////////////////////

static const int H2D_PAGE_SIZE = 32*1024; // 32K
static std::vector<char*> pages;
static int H2D_LAST_PAGE_SIZE;

static void init_mem()
{
  H2D_LAST_PAGE_SIZE = H2D_PAGE_SIZE;
  pages.clear();
}

static void* new_block(int size)
{
  // this is a simple and fast memory allocation function
  char* p;
  if (H2D_LAST_PAGE_SIZE + size <= H2D_PAGE_SIZE)
    p = pages.back();
  else {
    p = (char*) malloc(H2D_PAGE_SIZE);
    H2D_LAST_PAGE_SIZE = 0;
    pages.push_back(p);
  }
  p += H2D_LAST_PAGE_SIZE;
  H2D_LAST_PAGE_SIZE += size;
  return p;
}

static void free_mem()
{
  // dump all blocks at once - no need to walk through the linked lists etc.
  for (unsigned int i = 0; i < pages.size(); i++)
    free(pages[i]);
}

static inline MItem* new_item()
{
  std::string* str= new std::string;
  MItem* it = (MItem*) new_block(sizeof(MItem));
  it->next = NULL; 
  it->marker = str;
  return it; 
}

static inline MSymbol* new_symbol()
  { return (MSymbol*) new_block(sizeof(MSymbol)); }

static inline char* new_string(const char* str)
  { char* p = (char*) new_block((strlen(str)+1 +3) & ~3);
    return strcpy(p, str); }


//// symbol table //////////////////////////////////////////////////////////////////////////////////

static const int H2D_TABLE_SIZE = 256;  // must be a power of two
MSymbol* symbol_table[H2D_TABLE_SIZE];

static inline unsigned hash(const char* str)
{
  unsigned hash = 0;
  while (*str)
    hash = (*str++ & 0x1f) ^ (hash << 4);
  return hash;
}


static MSymbol* get_symbol(const char* name)
{
  unsigned index = hash(name) &(H2D_TABLE_SIZE-1);
  MSymbol* s = symbol_table[index];
  while (s != NULL && strcmp(s->name, name))
    s = s->next_hash;
  if (s != NULL) return s;

  s = new_symbol();
  s->name = new_string(name);
  s->data = NULL;
  s->built_in = NULL;
  s->next_hash = symbol_table[index];
  return symbol_table[index] = s;
}


static MSymbol* peek_symbol(const char* name)
{
  unsigned index = hash(name) &(H2D_TABLE_SIZE-1);
  MSymbol* s = symbol_table[index];
  while (s != NULL && strcmp(s->name, name))
    s = s->next_hash;
  return s;
}


static void add_built_in(const char* name, void* fn, int narg)
{
  MSymbol* s = get_symbol(name);
  s->built_in = fn;
  s->data = (MItem*) narg;
}

static void add_built_in(const char* name, MSymbolFunc1* fn) { add_built_in(name, (void*)fn, 1); }
static void add_built_in(const char* name, MSymbolFunc2* fn) { add_built_in(name, (void*)fn, 2); }

static void add_const(const char* name, double val)
{
  MSymbol* s = get_symbol(name);
  s->data = new_item();
  s->data->n = -1;
  s->data->val = val;
}


static void init_symbols()
{
  memset(symbol_table, 0, sizeof(symbol_table));

  // add the standard math library functions
  add_built_in("acos", acos);
  add_built_in("asin", asin);
  add_built_in("atan", atan);
  add_built_in("atan2", atan2);
  add_built_in("cos", cos);
  add_built_in("cosh", cosh);
  add_built_in("sin", sin);
  add_built_in("sinh", sinh);
  add_built_in("tan", tan);
  add_built_in("tanh", tanh);
  add_built_in("exp", exp);
  add_built_in("log", log);
  add_built_in("log10", log10);
  add_built_in("exp2", exp2);
#ifdef HAVE_LOG2
  add_built_in("log2", log2);
#endif
  add_built_in("pow", pow);
  add_built_in("sqrt", sqrt);
  add_built_in("cbrt", cbrt);
  add_built_in("hypot", hypot);
  add_built_in("ceil", ceil);
  add_built_in("abs", fabs);
  add_built_in("fabs", fabs);
  add_built_in("floor", floor);
  add_built_in("mod", fmod);
  add_built_in("fmod", fmod);

  // add the constants pi and PI
  add_const("pi", M_PI);
  add_const("PI", M_PI);
}




//// parser ////////////////////////////////////////////////////////////////////////////////////////

/*  The following is a simple recursive-descent parser for the following grammar:

    mesh   := { assig }
    assig  := ident "=" item [";"]
    item   := expr | list
    list   := "{" item {"," item} "}"

    expr   := list_ident | term | term "+" term | term "-" term
    term   := power | power "*" power | power "/" power
    power  := factor | factor "^" factor
    factor := number | var_ident | funct | "(" expr ")"
    funct  := funct_ident "(" args ")"
    args   := expr | expr "," args                    */


static MItem* item();
static MItem* item_string_marker();
static MItem* expression();


static double eval_fn(MSymbol* sym)
{
  double args[4];
  unsigned narg = 0;

  if (!follows(MT_KET))
  {
    while (narg < 4)
    {
      MItem* exp = expression();
      if (exp->n >= 0) serror("invalid use of list.");
      args[narg++] = exp->val;

      if (follows(MT_KET)) break;
      if (!follows(MT_COMMA)) serror("',' expected.");
      next_token();
    }
  }
  next_token();

  if (narg < (unsigned long) sym->data) serror("too few arguments to '%s'.", sym->name);
  if (narg > (unsigned long) sym->data) serror("too many arguments to '%s'.", sym->name);

  if (narg == 1)
    return ((double (*)(double)) sym->built_in)(args[0]);
  else if (narg == 2)
    return ((double (*)(double, double)) sym->built_in)(args[0], args[1]);

  error("Internal error.");
  return 0.0;
}


static double factor()
{
  if (follows(MT_BRA))
  {
    next_token();
    MItem* it = expression();
    if (it->n >= 0)
      serror("invalid use of list.");
    if (!follows(MT_KET))
      serror("')' expected.");
    next_token();
    return it->val;
  }
  else if (follows(MT_NUMBER))
  {
    double val = token->value;
    next_token();
    return val;
  }
  else if (follows(MT_IDENT))
  {
    MSymbol* sym = peek_symbol(token->text);
    if (sym == NULL)
      serror("unknown identifier '%s'.", token->text);
    next_token();

    if (follows(MT_BRA))
    {
      if (sym->built_in == NULL)
        serror("unknown function '%s'.", sym->name);
      next_token();
      return eval_fn(sym);
    }
    else
    {
      if (sym->built_in != NULL)
        serror("invalid use of function.");
      if (sym->data->n >= 0)
        serror("invalid use of list.");
      return sym->data->val;
    }
  }
  else {
    serror("syntax error: '%s'", token->text);
    return 0.0;
  }
}


static double power()
{
  double base = factor();
  if (follows(MT_POWER))
  {
    next_token();
    return pow(base, factor());
  }
  return base;
}


static double term()
{
  double result = power();
  while (follows(MT_STAR) || follows(MT_SLASH))
  {
    if (follows(MT_STAR)) { next_token(); result *= power(); }
                     else { next_token(); result /= power(); }
  }
  return result;
}


static MItem* expression()
{
  // special case: expression is a list reference
  if (follows(MT_IDENT))
  {
    MSymbol* list = peek_symbol(token->text);
    if (list != NULL && list->built_in == NULL &&
        list->data != NULL && list->data->n >= 0)
    {
      next_token();
      // only "=", "}", ";" or "," must follow
      if (!follows(MT_EQUAL) && !follows(MT_END) &&
          !follows(MT_COMMA) && !follows(MT_SEMICOL))
        serror("invalid use of list.");
      return list->data;
    }
  }

  // handle unary + and -
  double unary = 1.0;
  if (follows(MT_PLUS) || follows(MT_MINUS))
  {
    if (follows(MT_MINUS)) unary = -1.0;
    next_token();
  }

  double result = unary * term();
  while (follows(MT_PLUS) || follows(MT_MINUS))
  {
    if (follows(MT_PLUS)) { next_token(); result += term(); }
                     else { next_token(); result -= term(); }
  }

  MItem* it = new_item();
  it->n = -1;
  it->val = result;
  return it;
}


static void list(MItem* parent)
{
  parent->list = NULL;
  parent->n = 0;
  if (follows(MT_END)) { next_token(); return; }

  parent->list = item();
  parent->n++;
  MItem* it = parent->list;
  while (follows(MT_COMMA))
  {
    next_token();
    it->next = item();
    parent->n++;
    it = it->next;
  };

  if (follows(MT_SEMICOL))
  {
    bool one_word_processed = false;
    while(!follows(MT_END))
    {
      next_token();
      if(follows(MT_END))
        break;
      if(one_word_processed++)
        parent->marker->append(" ");
      it->next = item_string_marker();
      if(token->type == MT_NUMBER) {
	      std::ostringstream sin;
	      sin << token->value;
        parent->marker->append(sin.str());
      }
      else
        parent->marker->append(token->text);

      it = it->next;
    }
  }

  if (!follows(MT_END)) serror("'}' expected.");
  next_token();
}


static MItem* item()
{
  if (follows(MT_BEGIN))
  {
    next_token();
    MItem* it = new_item();
    list(it);
    return it;
  }
  else
    return expression();
}

static MItem* item_string_marker()
{
  MItem* it = new_item();
  it->n = -2;
  it->next = NULL;
  return it;
}


static void assignment(bool debug)
{
  if (!follows(MT_IDENT)) serror("identifier expected.");
  MSymbol* sym = get_symbol(token->text);
  next_token();

  if (!follows(MT_EQUAL)) serror("'=' expected.");
  next_token();

  sym->data = item();

  // semicolon is optional
  if (follows(MT_SEMICOL)) next_token();

  // print numerical temporary variables
  if (debug)
    if (sym->data->n < 0) {
      info("%s = %.18g", sym->name, sym->data->val);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void mesh_parser_init(FILE* f, const char* filename)
{
  mesh_lexer_init(f);
  H2D_PARSER_FILENAME = filename;
  init_mem();
  init_symbols();
}


void mesh_parser_run(bool debug)
{
  next_token();
  while (!follows(MT_EOF))
    assignment(debug);
}


MSymbol* mesh_parser_find_symbol(const char* name)
{
  return peek_symbol(name);
}


void mesh_parser_free()
{
  free_mem();
}


// checks that "list" is a list of "n" doubles and retrieves them
bool mesh_parser_get_doubles(MItem* list, int n, ...)
{
  if (list->n != n) return false;
  MItem* it = list->list;
  va_list ap; va_start(ap, n);
  for (int i = 0; i < n; i++, it = it->next) {
    if (it->n >= 0) return false;
    *(va_arg(ap, double*)) = it->val;
  }
  va_end(ap);
  return true;
}


// checks that "list" is a list of "n" integers and retrieves them
bool mesh_parser_get_ints(MItem* list, int n, ...)
{
  if (list->n != n) return false;
  MItem* it = list->list;
  va_list ap; va_start(ap, n);
  for (int i = 0; i < n; i++, it = it->next) {
    if (it->n >= 0 || !HERMES_IS_INT(it->val)) return false;
    *(va_arg(ap, int*)) = (int) it->val;
  }
  va_end(ap);
  return true;
}

void mitem_drop_string_markers(MItem* mi)
{
  if(mi->n < 1)
    delete mi->marker;
  else {
    for(int i = 0; i < mi->n; i++) {
      MItem* mi_pass = &mi->list[0];
      for(int j = 0; j < i; j++)
        mi_pass = mi_pass->next;
      mitem_drop_string_markers(mi_pass);
    }
  }
}