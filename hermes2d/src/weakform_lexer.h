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

#ifndef __H2D_WEAKFORM_LEXER_H
#define __H2D_WEAKFORM_LEXER_H


enum WFType
{
  T_ERROR, T_EOF,
  T_IDENT, T_NUMBER,
  T_VOL, T_SURF, T_COMMA, T_COLON,
  T_X, T_Y, T_XX, T_YY, T_XY,
  T_PLUS, T_MINUS, T_STAR, T_SLASH, T_BRA, T_KET,
  T_POWER, T_UNDER, T_IMAG
};


struct WFToken
{
  WFType type;
  double value;
  char*  lexeme;
  Token* next;
};


void   wf_lexer_init(const char* input);
Token* wf_get_token();
void   wf_lexer_free();


#endif
