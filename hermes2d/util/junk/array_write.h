

template<class T>
void Array::save_raw(FILE* f)
{
  int itemsize = sizeof(T), pagesize = PAGE_SIZE;
  hermes_fwrite("H2DA\001\000\000\000", 1, 8, f);
  hermes_fwrite(&size, sizeof(int), 1, f);
  hermes_fwrite(&nitems, sizeof(int), 1, f);
  hermes_fwrite(&append_only, sizeof(int), 1, f);
  hermes_fwrite(&itemsize, sizeof(int), 1, f);
  hermes_fwrite(&pagesize, sizeof(int), 1, f);

  int npages = pages.size();
  hermes_fwrite(&npages, sizeof(int), 1, f);
  for (int i = 0; i < npages; i++)
    hermes_fwrite(pages[i], itemsize, PAGE_SIZE, f);
}


template<class T>
void Array::load_raw(FILE* f)
{
  struct { char magic[4]; int ver; } hdr;
  hermes_fread(&hdr, 8, 1, f);
  if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'A')
    error("Not a Hermes2D Array file.");
  if (hdr.ver > 1)
    error("Unsupported file version.", filename);

  hermes_fread(&size, sizeof(int), 1, f);
  hermes_fread(&nitems, sizeof(int), 1, f);
  hermes_fread(&append_only, sizeof(int), 1, f);

  int itemsize, pagesize;
  hermes_fread(&itemsize, sizeof(int), 1, f);
  hermes_fread(&pagesize, sizeof(int), 1, f);
  if (itemsize != sizeof(T))
    error("Item size mismatch.");
  if (pagesize != PAGE_SIZE)
    error("Page size mismatch.");

  int i, npages;
  hermes_fread(&npages, sizeof(int), 1, f);
  pages.clear();
  for (i = 0; i < npages; i++)
  {
    T* new_page = new T[PAGE_SIZE];
    pages.push_back(new_page);
    hermes_fread(new_page, itemsize, PAGE_SIZE, f);
  }

  unused.clear();
  for (i = 0; i < size; i++)
    if (!get_item(i).used)
      unused.push_back(i);
}
