#include "../h2d_common.h"
#include "../range.h"
#include "element_to_refine.h"

using namespace std;

ElementToRefine& ElementToRefine::operator=(const ElementToRefine& orig) {
  id = orig.id;
  comp = orig.comp;
  split = orig.split;
  copy_orders(p, orig.p);
  copy_orders(q, orig.q);
  return *this;
}

HERMES_API std::ostream& operator<<(std::ostream& stream, const ElementToRefine& elem_ref) {
  stream << "id:" << elem_ref.id << ";comp:" << elem_ref.comp << "; split:" << get_refin_str(elem_ref.split) << "; orders:[";
  int num_sons = elem_ref.get_num_sons();
  for(int i = 0; i < num_sons; i++) {
    if (i > 0)
      stream << " ";
    stream << Hermes2D::get_quad_order_str(elem_ref.p[i]);
  }
  stream << "]";
  return stream;
}

HERMES_API std::istream& operator>>(std::istream& stream, const TagChecker& checker) {
  stringstream tag;
  for(unsigned i = 0; i < checker.get_tag().size(); i++) {
    char tag_char;
    stream >> tag_char;
    tag << tag_char;
  }
  if (checker.get_tag().compare(tag.str()) != 0) {
    stringstream str;
    if (stream.eof())
      str << "Unexpected EOF found while reading tag '" << checker.get_tag() << "'";
    else
      str << "Expected '" << checker.get_tag() << "' but '" << tag.str() << "' found at offset " << (int)stream.tellg();
    throw runtime_error(str.str());
  }
  return stream;
};

const char* ElementToRefineStream::H2DER_START_TAG = "ERS";
const char* ElementToRefineStream::H2DER_BIN_TAG = "BIN";
const char* ElementToRefineStream::H2DER_VECTOR_TAG = "ERV";

ElementToRefineStream::ElementToRefineStream(const char* filename, std::ios_base::openmode mode) {
  uint16_t test_value = ((uint16_t)'M' << 8) | 'L';
  if (*((char*)&test_value) == 'L')
    little_endian = true;
  else
    little_endian = false;

  open(filename, mode);
}

void ElementToRefineStream::open(const char* filename, std::ios_base::openmode mode) {
  error_if((mode & ios_base::binary) == 0, "Only binary mode is supported.");
  error_if((mode & ios_base::in) == 0 && (mode & ios_base::out) == 0 && (mode & ios_base::app) == 0, "Only in, out, and append mode is supported.");

  //if append: check header and append
  if ((mode & ios_base::app) != 0) {
    //open for reading
    stream.open(filename, (mode & ~ios_base::app) | ios_base::in);

    //if fails: write new file
    if (!stream.is_open()) {
      stream.open(filename, (mode & ~ios_base::app) | ios_base::out);
      error_if(!stream.is_open(), "Unable to open the stream \"%s\" for writing.", filename);
      write_header(mode);
    }
    else {
      //read header
      bool header_ok = read_header(mode);
      stream.close();

      //if header was fine: open as append
      if (header_ok)
        stream.open(filename, mode);
    }
  }
  else {
    //open stream
    stream.open(filename, mode);

    //check header (read)
    if (stream.good()) {
      if ((mode & ios_base::in) != 0)
        read_header(mode);
      else if ((mode & ios_base::out) != 0)
        write_header(mode);
    }
  }
}

void ElementToRefineStream::write_header(std::ios_base::openmode mode) {
  assert_msg((mode & ios_base::binary) != 0, "Binary mode supported only.");

  //write header tag
  stream << H2DER_START_TAG << " " << H2DER_BIN_TAG << "\n";
}

bool ElementToRefineStream::read_header(std::ios_base::openmode mode) {
  assert_msg((mode & ios_base::binary) != 0, "Binary mode supported only.");

  //decode
  try {
    //read header tag
    stream >> TagChecker(H2DER_START_TAG) >> skipws >> TagChecker(H2DER_BIN_TAG) >> skipws;
    return true;
  }
  catch(std::runtime_error& err) {
    error("Invalid file tag or unsupported file format (%s)", err.what());
    return false;
  }
}

uint8_t ElementToRefineStream::get_byte_size(int value) {
  if (value == 0)
    return 1;
  else {
    int value_abs = abs(value);
    double byte_len = ceil((log2(value_abs) + 1) / 8);
    error_if(byte_len > 255, "Required calculation of byte size of %d but the size is larger than 256 bytes", value);
    return (uint8_t)byte_len;
  }
}

void ElementToRefineStream::write_bytes(const int integer, int num_bytes) {
  write_bytes(&integer, num_bytes);
}

void ElementToRefineStream::write_bytes(const void* integer_ptr, int num_bytes) {
  if (little_endian)
    stream.write((const char*)integer_ptr, num_bytes);
  else {
    const char* integer_ptr_c = (const char*)integer_ptr;
    for(unsigned i = sizeof(int) - 1; i > (sizeof(int) - num_bytes - 1); i--)
      stream << integer_ptr_c[i];
  }
}

int ElementToRefineStream::read_bytes(int num_bytes) {
  error_if((unsigned) num_bytes > sizeof(int), "Requested number of bytes (%d) exceedes size of integer (%d)", num_bytes, sizeof(int));
  int shift = 0, result = 0;
  uint8_t buffer;
  for(int i = 0; i < num_bytes; i++) {
    stream.read((char*)&buffer, 1);
    int component = buffer;
    if (shift > 0)
      component <<= shift;
    result |= component;
    shift += 8;
  }

  //expand sign
  if ((buffer & 0x80) != 0 && (unsigned) num_bytes < sizeof(int)) {
    int sign = -1 << (num_bytes * 8);
    result |= sign;
  }

  return result;
}

HERMES_API ElementToRefineStream& operator<<(ElementToRefineStream& stream, 
           const std::vector<ElementToRefine>& elem_refs) {
  //calculate range of values
  Range<int> range_id(0, 0), range_comp(0, 0), range_order(0, 0);
  vector<ElementToRefine>::const_iterator elem_ref = elem_refs.begin();
  while (elem_ref != elem_refs.end()) {
    range_id.enlarge_to_include(elem_ref->id);
    range_comp.enlarge_to_include(elem_ref->comp);

    const int num_sons = elem_ref->get_num_sons();
    for(int k = 0; k < num_sons; k++) {
      range_order.enlarge_to_include(H2D_GET_H_ORDER(elem_ref->p[k]));
      range_order.enlarge_to_include(H2D_GET_V_ORDER(elem_ref->p[k]));
    }

    //next element
    elem_ref++;
  }

  //calculate sizes
  uint8_t num_size = ElementToRefineStream::get_byte_size(elem_refs.size());
  uint8_t id_root_size = ElementToRefineStream::get_byte_size(range_id.lower());
  uint8_t id_offset_size = ElementToRefineStream::get_byte_size(range_id.upper() - range_id.lower());
  uint8_t comp_root_size = ElementToRefineStream::get_byte_size(range_comp.lower());
  uint8_t comp_offset_size = 0;
  if (range_comp.upper() != range_comp.lower())
    comp_offset_size = ElementToRefineStream::get_byte_size(range_comp.upper() - range_comp.lower());
  uint8_t order_size = ElementToRefineStream::get_byte_size(range_order.upper());

  //store header and sizes
  stream.stream << ElementToRefineStream::H2DER_VECTOR_TAG; //tag (3 bytes)
  stream.write_bytes(num_size, ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of number of elements
  stream.write_bytes(id_root_size, ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of ID root
  stream.write_bytes(id_offset_size, ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of ID offset
  stream.write_bytes(comp_root_size, ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of component root
  stream.write_bytes(comp_offset_size, ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of component offset; if zero, component is skipped in record
  stream.write_bytes(order_size, ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of order

  //store roots
  stream.write_bytes((int)elem_refs.size(), num_size); // number of elements
  int id_root = range_id.lower();
  stream.write_bytes(id_root, id_root_size); // ID root
  int comp_root = range_comp.lower();
  stream.write_bytes(comp_root, comp_root_size); // component root

  //store individual records
  elem_ref = elem_refs.begin();
  while (elem_ref != elem_refs.end()) {
    stream.write_bytes(elem_ref->id - id_root, id_offset_size); // ID offset
    if (comp_offset_size > 0)
      stream.write_bytes(elem_ref->comp - comp_root, comp_offset_size); // component offset
    stream.write_bytes(elem_ref->split, 1); // split type

    //orders
    const int num_sons = elem_ref->get_num_sons();
    for(int i = 0; i < num_sons; i++) {
      stream.write_bytes(H2D_GET_H_ORDER(elem_ref->p[i]), order_size);
      stream.write_bytes(H2D_GET_V_ORDER(elem_ref->p[i]), order_size);
    }

    //next element
    elem_ref++;
  }
  return stream;
}

HERMES_API ElementToRefineStream& operator>>(ElementToRefineStream& stream, std::vector<ElementToRefine>& elem_refs) {
  int pos = (int)stream.stream.tellg();
  //read tag
  try { stream.stream >> TagChecker(ElementToRefineStream::H2DER_VECTOR_TAG); }
  catch (runtime_error &err) { error("Unable to read record start tag (%s)", err.what()); };

  //read sizes
  int num_size = stream.read_bytes(ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of number of elements
  int id_root_size = stream.read_bytes(ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of ID root
  int id_offset_size = stream.read_bytes(ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of ID offset
  int comp_root_size = stream.read_bytes(ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of component root
  int comp_offset_size = stream.read_bytes(ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of component offset; if zero, component is skipped in record
  int order_size = stream.read_bytes(ElementToRefineStream::H2DER_SIZE_BYTESIZE); // byte length of order

  //read roots
  int num_elems = stream.read_bytes(num_size); // number of elements
  int id_root = stream.read_bytes(id_root_size); // ID root
  int comp_root = stream.read_bytes(comp_root_size); // component root

  //prepare space
  elem_refs.clear();
  elem_refs.reserve(num_elems);

  //read refinements
  for(int i = 0; i < num_elems; i++) {
    ElementToRefine elem_ref;
    elem_ref.id = id_root + stream.read_bytes(id_offset_size); // ID offset
    elem_ref.comp = comp_root;
    if (comp_offset_size > 0)
      elem_ref.comp += stream.read_bytes(comp_offset_size); // component offset
    elem_ref.split = stream.read_bytes(1); // split type

    //orders
    memset(elem_ref.p, 0, sizeof(elem_ref.p));
    memset(elem_ref.q, 0, sizeof(elem_ref.q));
    const int num_sons = elem_ref.get_num_sons();
    for(int k = 0; k < num_sons; k++) {
      int order_h = stream.read_bytes(order_size);
      int order_v = stream.read_bytes(order_size);
      elem_ref.p[k] = H2D_MAKE_QUAD_ORDER(order_h, order_v);
    }

    //store element
    elem_refs.push_back(elem_ref);
  }

  return stream;
}
