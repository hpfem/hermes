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

#ifndef __H2D_MESH_PARSER_H
#define __H2D_MESH_PARSER_H


struct MItem
{
  int n; ///< -1..numeric value, otherwise list of length n
  union { double val; MItem* list; };
  MItem* next; ///< next item in the list this MItem is part of
};

typedef double (MSymbolFunc1)(double);
typedef double (MSymbolFunc2)(double, double);

struct MSymbol
{
  const char* name;   ///< symbol name
  MItem* data;        ///< pointer to symbol value (numeric or list)
  void* built_in;     ///< internal
  MSymbol* next_hash; ///< internal
};


void mesh_parser_init(FILE* f, const char* filename);
void mesh_parser_run(bool debug);
void mesh_parser_free();

MSymbol* mesh_parser_find_symbol(const char* name);

bool mesh_parser_get_doubles(MItem* list, int n, ...);
bool mesh_parser_get_ints(MItem* list, int n, ...);


#endif
