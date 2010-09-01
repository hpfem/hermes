#ifndef __ASSOC_H
#define __ASSOC_H

// This may one day become a substitute for the slow Judy arrays.
// Nothing of it is currently used.


#include "common.h"


template<typename KeyType>
struct AssocKV
{
  int size, cap;
  struct { KeyType key; void* val; } kv[0];
};


const int H2D_DEFAULT_ASSOC_HASH_SIZE = 32;


template<typename KeyType>
void** assoc_ins(void** array, KeyType key, int hash_size = H2D_DEFAULT_ASSOC_HASH_SIZE)
{
  typedef AssocKV<KeyType> KV;
  KV** ht = (KV**) *array;

  if (ht == NULL)
  {
    // create a new hash table
    assert((hash_size & (hash_size-1)) == 0); // must be a power of 2
    ht = (KV**) (*array = malloc(hash_size * sizeof(KV*)));
    memset(ht, 0, hash_size * sizeof(KV*));
  }

  //KV** pp = ht + ((key * (KeyType) 0x29F5C7612A3B39CE) & (hash_size-1));
  KV** pp =  ht + ((key * 0x29F5C765) & (hash_size-1));
  KV*  kv = *pp;

  const int init_cap = 2;
  const int pair_size = sizeof(KeyType) + sizeof(void*);

  if (kv == NULL)
  {
    // create a new key-value table
    kv = *pp = (KV*) malloc(sizeof(KV) + init_cap*pair_size);
    kv->size = 1;
    kv->cap  = init_cap;
    kv->kv[0].key = key;
    kv->kv[0].val = NULL;
    return &(kv->kv[0].val);
  }
  else
  {
    // try finding an existing record in the key-value table
    int lo = 0, hi = kv->size-1, mid;
    while (1)
    {
      mid = (lo + hi) >> 1;
      KeyType midkey = kv->kv[mid].key;
      if (key < midkey)
        hi = mid-1;
      else if (key > midkey)
        lo = mid+1;
      else
        return &(kv->kv[mid].val); // found it!

      if (lo > hi)
      {
        // not found -- insert a new key-value pair
        if (key > midkey) mid++;

        if (kv->size >= kv->cap)
        {
          // reallocate the table if necessary
          kv->cap *= 2;
          kv = *pp = (KV*) realloc(kv, sizeof(KV) + kv->cap*pair_size);
        }

        // shift the rest and insert the new pair
        memmove(kv->kv + mid+1, kv->kv + mid, (kv->size - mid)*pair_size);
        kv->size++;
        kv->kv[mid].key = key;
        kv->kv[mid].val = NULL;
        return &(kv->kv[mid].val);
      }
    }
  }
}


template<typename KeyType>
void assoc_forall(void** array, void (*fn)(KeyType key, void* val), int hash_size = H2D_DEFAULT_ASSOC_HASH_SIZE)
{
  typedef AssocKV<KeyType> KV;
  KV** ht = (KV**) *array;
  if (ht == NULL) return;

  for (int i = 0; i < hash_size; i++)
  {
    KV* kv = ht[i];
    if (kv != NULL)
      for (int j = 0; j < kv->size; j++)
        fn(kv->kv[j].key, kv->kv[j].val);
  }
}


template<typename KeyType>
void assoc_free(void** array, int hash_size = H2D_DEFAULT_ASSOC_HASH_SIZE)
{
  typedef AssocKV<KeyType> KV;
  KV** ht = (KV**) *array;
  if (ht == NULL) return;

  for (int i = 0; i < hash_size; i++)
    if (ht[i] != NULL)
      free(ht[i]);

  free(ht);
  *ht = NULL;
}


template<typename KeyType>
void assoc_dump(void** array, int hash_size = H2D_DEFAULT_ASSOC_HASH_SIZE)
{
  typedef AssocKV<KeyType> KV;
  KV** ht = (KV**) *array;
  if (ht == NULL) return;

  for (int i = 0; i < hash_size; i++)
  {
    KV* kv = ht[i];
    printf("#%d: ", i);
    if (kv != NULL)
      for (int j = 0; j < kv->size; j++)
      {
        printf("{%d,%d}", kv->kv[j].key, (int) kv->kv[j].val);
        if (j < kv->size-1) printf(", ");
      }
    printf("\n");
  }
}



/*

void** l3_ins(void** array, int key)
{
  const int size = 16 * sizeof(void*);
  if (*array == NULL) { array = malloc(size); memset(array, 0, size); }
  if (key < 15) return *array + key;
  if (array[15] == NULL) { array[15] = malloc(size); memset(array, 0, size); }
  return array[15] + (key-15);
}

void** l3_next(void** array, int* key)
{
  void** result;
  for ( ; *key < 30; *key++)
    if (*(result = l3_ins(array, *key)) != NULL)
      return result;
  return NULL;
}

void** l3_first(void** array, int* key)
{
  *key = 0;
  return l3_next(array, key);
}

*/



#endif
