#include "hermes2d.h"
#include "mesh_parser.h"


void print_item(MItem* it)
{
  if (it->n < 0)
  {
    printf("%g", it->val);
  }
  else
  {
    printf("{ ");
    MItem* it2 = it->list;
    for (int i = 0; i < it->n; i++, it2 = it2->next) {
      print_item(it2);
      if (i < it->n-1) printf(", ");
    }
    printf(" }");
  }
}


int main(int argc, char* argv[])
{
  const char* name = "/home/jakub/hermes2d/util/junk/bracket.mesh";
  FILE* f = fopen(name, "r");
  mesh_parser_init(f, name);
  mesh_parser_run();

  MSymbol* sym = mesh_parser_find_symbol("vertices");
  print_item(sym->data);
  printf("\n");

  mesh_parser_free();
  return 0;
}
