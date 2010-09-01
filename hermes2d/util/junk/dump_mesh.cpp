

void dump_mesh(Mesh* mesh)
{
  int i;
  printf("----------------------------------------------------------------\n\n");

  printf("nodes:\n");
  Node* n;
  for_all_nodes(n, mesh)
  {
    printf("id=%2d type=%s ref=%d bnd=%d p1=%2d p2=%2d ", n->id,
           (n->type == TYPE_VERTEX) ? "vert" : "edge", n->ref, n->bnd, n->p1, n->p2);

    if (n->type == TYPE_VERTEX)
      printf("x=%g y=%g\n", n->x, n->y);
    else
      printf("mrk=%d elem0=%d elem1=%d\n", n->marker,
             (n->elem[0] == NULL) ? -1 : n->elem[0]->id,
             (n->elem[1] == NULL) ? -1 : n->elem[1]->id);
  }

  printf("\nelements:\n");
  Element* e;
  for_all_elements(e, mesh)
  {
    printf("id=%2d nvert=%d act=%d ", e->id, e->nvert, e->active);
    if (e->active)
    {
      for (i = 0; i < e->nvert; i++)
        printf("v%d=%2d ", i, e->vn[i]->id);
      for (i = 0; i < e->nvert; i++)
        printf("e%d=%2d ", i, e->en[i]->id);
    }
    else
    {
      for (i = 0; i < 4; i++)
        printf("son%d=%2d ", i, (e->sons[i] == NULL) ? -1 : e->sons[i]->id);
    }
    printf("\n");
  }
  printf("\n");
}

