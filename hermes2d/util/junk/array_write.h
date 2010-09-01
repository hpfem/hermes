

template<class T>
void Array::save_raw(FILE* f)
{
  int itemsize = sizeof(T), pagesize = PAGE_SIZE;
  hermes2d_fwrite("H2DA\001\000\000\000", 1, 8, f);
  hermes2d_fwrite(&size, sizeof(int), 1, f);
  hermes2d_fwrite(&nitems, sizeof(int), 1, f);
  hermes2d_fwrite(&append_only, sizeof(int), 1, f);
  hermes2d_fwrite(&itemsize, sizeof(int), 1, f);
  hermes2d_fwrite(&pagesize, sizeof(int), 1, f);

  int npages = pages.size();
  hermes2d_fwrite(&npages, sizeof(int), 1, f);
  for (int i = 0; i < npages; i++)
    hermes2d_fwrite(pages[i], itemsize, PAGE_SIZE, f);
}


template<class T>
void Array::load_raw(FILE* f)
{
  struct { char magic[4]; int ver; } hdr;
  hermes2d_fread(&hdr, 8, 1, f);
  if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'A')
    error("Not a Hermes2D Array file.");
  if (hdr.ver > 1)
    error("Unsupported file version.", filename);

  hermes2d_fread(&size, sizeof(int), 1, f);
  hermes2d_fread(&nitems, sizeof(int), 1, f);
  hermes2d_fread(&append_only, sizeof(int), 1, f);

  int itemsize, pagesize;
  hermes2d_fread(&itemsize, sizeof(int), 1, f);
  hermes2d_fread(&pagesize, sizeof(int), 1, f);
  if (itemsize != sizeof(T))
    error("Item size mismatch.");
  if (pagesize != PAGE_SIZE)
    error("Page size mismatch.");

  int i, npages;
  hermes2d_fread(&npages, sizeof(int), 1, f);
  pages.clear();
  for (i = 0; i < npages; i++)
  {
    T* new_page = new T[PAGE_SIZE];
    pages.push_back(new_page);
    hermes2d_fread(new_page, itemsize, PAGE_SIZE, f);
  }

  unused.clear();
  for (i = 0; i < size; i++)
    if (!get_item(i).used)
      unused.push_back(i);
}
