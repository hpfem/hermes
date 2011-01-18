#ifndef __H2D_AUTO_LOCAL_ARRAY_H
#define __H2D_AUTO_LOCAL_ARRAY_H

#include <stdlib.h>

/// \brief Allocation of a local array on the stack.
// An array of ordinal types or structures.
//#define AUTOLA_OR(__type, __name, __len) AutoLocalArray<__type> __name ((__type*)alloca((__len) * sizeof(__type)), __len)
#define AUTOLA_OR(__type, __name, __len) AutoLocalArray<__type> __name (new __type[__len], __len)
// An array of classes that has a default constructor that allocates memery.
#define AUTOLA_CL(__type, __name, __len) AutoLocalClassArray<__type> __name (new __type[__len], __len)
// An array of ordinal type or structures.
//#define AUTOLA2_OR(__type, __name, __len0, __len1) AutoLocalArray2<__type> __name ((__type*)alloca((__len0) * (__len1) * sizeof(__type)), __len0, __len1)
#define AUTOLA2_OR(__type, __name, __len0, __len1) AutoLocalArray2<__type> __name (new __type[(__len0) * (__len1)], __len0, __len1)

/// \brief A support for local arrays allocated on the stack. The class constains just information about the size.
template<class T>
class AutoLocalArray
{
protected:
  T* arr;
public:
  const size_t length;
  const size_t size;

  AutoLocalArray(T* arr, size_t length) : arr(arr), length(length), size(length*sizeof(T)) {};
  inline T * operator*() { return arr; };
  inline T & operator[](const int inx) { return arr[inx]; };

  inline operator T*() { return arr; };
  inline operator void*() { return arr; };

  virtual ~AutoLocalArray() { delete[] arr; arr = NULL; };
};

/// \brief A support for local arrays of class instances. The class constains just information about the size.
template<class T>
class AutoLocalClassArray : public AutoLocalArray<T>
{
public:
  AutoLocalClassArray(T* arr, size_t length) : AutoLocalArray<T>(arr, length) {};
  virtual ~AutoLocalClassArray() { delete[] this->arr; this->arr = NULL; }
};

/// \brief A support for local 2-dimensional arrays allocated on the stack. The class constains just information about the size.
template<class T>
class AutoLocalArray2
{
protected:
  T* arr;
public:
  const int length0;
  const int length1;
  const int size;

  AutoLocalArray2(T* arr, int length0, int length1) : arr(arr), length0(length0), length1(length1), size(length0*length1*sizeof(T)) {};
  inline T * operator*() { return arr; };
  inline T * operator[](const int inx) { return &arr[inx*length1]; };

  inline operator T*() { return arr; };
  inline operator void*() { return arr; };

  virtual ~AutoLocalArray2() { delete[] arr; arr = NULL; };
};



#endif
