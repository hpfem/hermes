


void Filter::precalculate_sln(int i, int order, int mask)
{
  const Trf* top = sln[i]->get_top_transform();
  trf[i].m[0] = top->m[0] * ctm->m[0];
  trf[i].m[1] = top->m[1] * ctm->m[1];
  trf[i].t[0] = top->m[0] * ctm->t[0] + top->t[0];
  trf[i].t[1] = top->m[1] * ctm->t[1] + top->t[1];
  sln[i]->force_transform((unidata[i][element->id].idx << 3) + sub_idx + 1, trf + i);
  sln[i]->set_quad_order(order, item[i]);
}
