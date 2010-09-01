#include "hermes2d.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

#define MAX_BUFFER 1024 ///< A maximum lenght of the buffer

/// Dumps & compares a line.
int dump_compare_line(FILE* file, const char* line, const int line_inx) {
  if (file == NULL) {
    printf("%s", line);
    return ERROR_SUCCESS;
  }
  else {
    char buffer[MAX_BUFFER];
    if (fgets(buffer, MAX_BUFFER-1, file) == NULL) {
      error("unable to read line %d from the dump file", line_inx);
      return ERROR_FAILURE;
    }
    size_t len = strlen(buffer);
    if (buffer[len-1] == '\n')
      buffer[len-1] = '\0';
    if (strcmp(line, buffer) == 0)
      return ERROR_SUCCESS;
    else
      return ERROR_FAILURE;
  }
}

/// Compares the current line with a line in the dump file.
#define DUMP_CMP(__line) if (dump_compare_line(file, __line, line_cnt) == ERROR_FAILURE) goto quit; \
  line_cnt++; line_inx = 0
/// Prints the parameters to the line buffer.
#define DUMP_OUT(__line, ...) { char buffer[MAX_BUFFER]; \
  sprintf(buffer, __VA_ARGS__); \
  size_t sz_buf = strlen(buffer); \
  if ((sz_buf + line_inx) < MAX_BUFFER) { strcpy(__line + line_inx, buffer); line_inx += sz_buf; } \
  else { error("output line exceeds the maximum size of %d characters", MAX_BUFFER); } }

/// Compares dump with the mesh. If file_name_dump is NULL, it just prints the output.
int dump_compare(const Mesh &mesh, const char* file_name_dump) {
  int result = ERROR_FAILURE;

  //open source file if any
  FILE* file = NULL;
  if (file_name_dump != NULL) {
    file = fopen(file_name_dump, "rt");
    if (file == NULL) { error("unable to open the dump file \"%s\"", file_name_dump); goto quit; }
  }

  //dump/compare
  {
    char line[MAX_BUFFER];
    unsigned line_inx = 0;
    int line_cnt = 0;

    int ne = mesh.get_num_elements();
    DUMP_OUT(line, "Elements = %d", ne);
    DUMP_CMP(line);
    for (int eid = 0; eid < ne; eid++)
    {
      Element *e = mesh.get_element(eid);
      DUMP_OUT(line, " #%d:", e->id);
      for (unsigned iv = 0; iv < e->nvert; iv++)
      {
        DUMP_OUT(line, " %d", e->vn[iv]->id);
      }
      DUMP_OUT(line, " | %d", e->marker);
      DUMP_CMP(line);
    }

    int im = 0;
    DUMP_OUT(line, "Markers");
    DUMP_CMP(line);
    for (int eid = 0; eid < ne; eid++)
    {
      Element *e = mesh.get_element(eid);
      if (!e->active)
        continue;

      int nv = e->nvert;
      for (int iv = 0; iv < nv; iv++)
      {
        Node *nd = e->en[iv];
        if (nd->type == 1 && nd->bnd == 1)
        { // edge node
          DUMP_OUT(line, " %d, %d | %d", e->vn[iv]->id, e->vn[(iv + 1) % nv]->id, nd->marker);
          DUMP_CMP(line);
        }
      }
    }

    DUMP_OUT(line, "Curvilinears");
    DUMP_CMP(line);
    for (int eid = 0; eid < ne; eid++)
    {
      Element *e = mesh.get_element(eid);
      if (!e->active)
        continue;

      if (!e->is_curved())
        continue;

      int nv = e->nvert;
      for (int iv = 0; iv < nv; iv++)
      {
        Node *nd = e->en[iv];
        if (e->cm->nurbs[iv] == NULL)
          continue;  

        if (nd->type == 1 && nd->bnd == 1)
        { 
          DUMP_OUT(line, " %d, %d | %.6lf", e->vn[iv]->id, e->vn[(iv + 1) % nv]->id, e->cm->nurbs[iv]->angle);
          DUMP_CMP(line);
        }
      }
    }
  }

  result = ERROR_SUCCESS;

quit: //finish
  if (file != NULL)
    fclose(file);
  return result;
}

int main(int argc, char* argv[])
{
  int ret = ERROR_FAILURE;

  if (argc < 2)
  {
    printf("please input as this format: <mesh type> <meshfile> [meshfiledump] \n");
    return ERROR_FAILURE;
  }

  char *mtype = argv[1];
  char *file_name = argv[2];
  char *file_name_dump = NULL;
  if (argc > 2)
    file_name_dump = argv[3];

  // load the mesh file
  Mesh mesh;
  MeshLoader *mloader = NULL;
  if (strcmp(mtype, "exII") == 0) mloader = new ExodusIIReader();
  else if (strcmp(mtype, "h2d") == 0) mloader = new H2DReader();
  else if (strcmp(mtype, "h2d-str") == 0) {
    // load the file into a string
    FILE *file = fopen(file_name , "rb");
    error_if(file == NULL, "unable to open file '%s'", file_name);

    // obtain file size:
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    rewind(file);

    // allocate memory to contain the whole file:
    char *buffer = (char *) malloc (sizeof(char) * size);
    error_if(buffer == NULL, "memory error");

    // copy the file into the buffer:
    size_t result = fread(buffer, 1, size, file);
    error_if(result != size, "reading error");

    fclose(file);

    //
    H2DReader *hloader = new H2DReader();
    hloader->load_str(buffer, &mesh);
    ret = dump_compare(mesh, file_name_dump);
    delete hloader;
    free(buffer);
    return ret;
  }
  else if (strcmp(mtype, "h2d-old") == 0) {
    H2DReader *hloader = new H2DReader();
    hloader->load_old(file_name, &mesh);
    ret = dump_compare(mesh, file_name_dump);
    delete hloader;
    return ret;
  }
  else {
    error("unknown mesh loader type");
  }

  if (mloader->load(file_name, &mesh))
  {
    ret = dump_compare(mesh, file_name_dump);
  }
  else
  {
    error("failed");
    ret = ERROR_FAILURE;
  }

  delete mloader;

  return ret;
}

