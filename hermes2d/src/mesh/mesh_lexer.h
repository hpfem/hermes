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

#ifndef __H2D_MESH_LEXER_H
#define __H2D_MESH_LEXER_H


enum MTType
{
  MT_ERROR, MT_EOF,
  MT_IDENT, MT_NUMBER,
  MT_EQUAL, MT_SEMICOL, MT_COMMA,
  MT_PLUS, MT_MINUS, MT_STAR, MT_SLASH,
  MT_BRA, MT_KET, MT_POWER,
  MT_BEGIN, MT_END
};


struct MToken
{
  MTType type;
  double value;
  const char* text;
};


void    mesh_lexer_init(FILE* f);
MToken* mesh_get_token();

extern int mesh_lexer_line_num;


#endif
