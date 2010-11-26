#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Testing templates of Judy arrays

bool testPrint(bool value, const char *msg, bool correct) {
	info("%s...", msg);
	if (value == correct) {
		info("OK.");
		return true;
	}
	else {
		info("failed.");
		return false;
	}
}

//
// Test Array template
//

int testArrayInt() {
	info("- Testing Array<int>-----");
	Array<int> int_array;
	bool r;
	unsigned int idx;

	// fill the array
	info("  * Filling the array with numbers.");
	for (int i = 0; i < 100; i++)
		int_array.set(i, i + 100);

	// test the values
	r = true;
	for (int i = 0; i < 100; i++)
		r &= int_array.get(i) == (i + 100);
	if (!testPrint(r, "  * Checking if all values were inserted correctly", true))
		return -1;

	// overwriting an existing item (checking memory leaks)
	int_array.set(0, 100);
	if (!testPrint(int_array[0] == 100, "  * Replacing existing value", true))
		return -1;

	// test [] operator
	int_array[1000] = 999;
	if (!testPrint(int_array[1000] == 999, "  * Setting value by [] operator", true))
		return -1;

	int_array[1000] = 1100;
	if (!testPrint(int_array[1000] == 1100, "  * Replacing existing value by [] operator", true))
		return -1;

	idx = int_array.by_count(101);
	if (!testPrint(int_array[1000] == int_array[idx], "  * Testing item 101th item", true))
		return -1;

	//
	info("  * Testing iterators.");
	r = true;
	for (unsigned int i = int_array.first(); i != INVALID_IDX; i = int_array.next(i))
		r &= int_array.get(i) == (int) (i + 100);
	if (!testPrint(r, "    - Forward iteration", true))
		return -1;

	r = true;
	for (unsigned int i = int_array.last(); i != INVALID_IDX; i = int_array.prev(i))
		r &= int_array.get(i) == (int) (i + 100);
	if (!testPrint(r, "    - Backward iteration", true))
		return -1;

	idx = int_array.add(2000);
	if (!testPrint(int_array[idx] == 2000, "  * Adding an item at the end of the array", true))
		return -1;

	// remove
	int_array.remove(idx);
	if (!testPrint(!int_array.exists(idx), "  * Removing an item from the array", true))
		return -1;

	// MEMORY
	unsigned int mem_used;

	// mem used
	mem_used = int_array.mem_used();
	if (!testPrint(mem_used > 2, "  * Used memory should be nonzero", true))
		return -1;

	// free memory
	int_array.remove_all();
	testPrint(true, "  * Freeing memory", true);

	mem_used = int_array.mem_used();
	if (!testPrint(mem_used == 0, "  * Used memory should be zero", true))
		return -1;

	return 0;
}

struct Point {
	int X, Y;

	Point() {
		X = 0;
		Y = 0;
	}

	Point(int x, int y) {
		X = x;
		Y = y;
	}

	bool operator==(Point &o) {
		return (X == o.X) && (Y == o.Y);
	}
};

int testArrayStruct() {
	info("- Testing Array<struct>-----");

	Array<Point> pt_array;
	bool r;
	unsigned int idx;

	// fill the array
	info("  * Filling the array with numbers.");
	for (int i = 0; i < 100; i++) {
		pt_array.set(i, Point(100 + i, 1000 + i));
	}

	// test the values
	r = true;
	for (int i = 0; i < 100; i++) {
		r &= (pt_array.get(i).X == i + 100) && (pt_array.get(i).Y == i + 1000);
	}
	if (!testPrint(r, "  * Checking if all values were inserted correctly", true))
		return -1;

	// overwriting an existing item (checking memory leaks)
	pt_array.set(0, Point(100, 1000));
	if (!testPrint((pt_array[0].X == 100) && (pt_array[0].Y == 1000), "  * Replacing existing value", true))
		return -1;

	// test [] operator
	pt_array[1000] = Point(199, 1999);
	if (!testPrint((pt_array[1000].X == 199) && (pt_array[1000].Y == 1999), "  * Setting value by [] operator", true))
		return -1;

	pt_array[1000] = Point(1100, 2000);
	if (!testPrint((pt_array[1000].X == 1100) && (pt_array[1000].Y == 2000), "  * Replacing existing value by [] operator", true))
		return -1;

	idx = pt_array.by_count(101);
	if (!testPrint(pt_array[1000] == pt_array[idx], "  * Testing item 101th item", true))
		return -1;

	//
	info("  * Testing iterators.");
	r = true;
	for (unsigned int i = pt_array.first(); i != INVALID_IDX; i = pt_array.next(i)) {
		r &= (pt_array.get(i).X == (int) (i + 100)) && (pt_array.get(i).Y == (int) (i + 1000));
	}
	if (!testPrint(r, "    - Forward iteration", true))
		return -1;

	r = true;
	for (unsigned int i = pt_array.last(); i != INVALID_IDX; i = pt_array.prev(i)) {
		r &= (pt_array.get(i).X == (int) (i + 100)) && (pt_array.get(i).Y == (int) (i + 1000));
	}
	if (!testPrint(r, "    - Backward iteration", true))
		return -1;

	idx = pt_array.add(Point(2000, 3000));
	if (!testPrint((pt_array[idx].X == 2000) && (pt_array[idx].Y == 3000), "  * Adding an item at the end of the array", true))
		return -1;

	// remove
	pt_array.remove(idx);
	if (!testPrint(!pt_array.exists(idx), "  * Removing an item from the array", true))
		return -1;

	// MEMORY
	unsigned int mem_used;

	// mem used
	mem_used = pt_array.mem_used();
	if (!testPrint(mem_used > 2, "  * Used memory should be nonzero", true))
		return -1;

	// free memory
	pt_array.remove_all();
	testPrint(true, "  * Freeing memory", true);

	mem_used = pt_array.mem_used();
	if (!testPrint(mem_used == 0, "  * Used memory should be zero", true))
		return -1;

	return 0;
}

// test array templates (with int and with struct)
int testArray() {
	int ret;

	// test the array of ints
	if ((ret = testArrayInt()) != 0)
		return ret;

	// test the array of structs
	if ((ret = testArrayStruct()) != 0)
		return ret;

	return 0;
}

//
// Test ArrayPtr template
//

int testArrayPtrStruct() {
	info("- Testing ArrayPtr<struct>-----");

	ArrayPtr<Point> ptr_array;
	bool r;
	unsigned int idx;

	// fill the array
	info("  * Filling the array with items.");
	for (int i = 0; i < 100; i++) {
		Point *pt = new Point(100 + i, 1000 + i);
		ptr_array.set(i, pt);
	}

	// test the values
	r = true;
	for (int i = 0; i < 100; i++) {
		Point *pt = ptr_array.get(i);
		r &= (pt->X == i + 100) && (pt->Y == i + 1000);
	}

	if (!testPrint(r, "  * Checking if all values were inserted correctly", true))
		return -1;

	// test [] operator
	ptr_array[1000] = new Point(999, 9998);
	if (!testPrint((ptr_array[1000]->X == 999) && (ptr_array[1000]->Y == 9998), "  * Setting value by [] operator", true))
		return -1;

	delete ptr_array.get(1000);
	ptr_array[1000] = new Point(1100, 2000);
	if (!testPrint((ptr_array[1000]->X == 1100) && (ptr_array[1000]->Y == 2000), "  * Replacing existing value by [] operator", true))
		return -1;

	idx = ptr_array.by_count(101);
	if (!testPrint(ptr_array[1000] == ptr_array[idx], "  * Testing item 101th item", true))
		return -1;

	//
	info("  * Testing iterators.");
	r = true;
	for (unsigned int i = ptr_array.first(); i != INVALID_IDX; i = ptr_array.next(i)) {
		Point *pt = ptr_array.get(i);
		r &= (pt->X == (int) (i + 100)) && (pt->Y == (int) (i + 1000));
	}
	if (!testPrint(r, "    - Forward iteration", true))
		return -1;

	r = true;
	for (unsigned int i = ptr_array.last(); i != INVALID_IDX; i = ptr_array.prev(i)) {
		Point *pt = ptr_array.get(i);
		r &= (pt->X == (int) (i + 100)) && (pt->Y == (int) (i + 1000));
	}
	if (!testPrint(r, "    - Backward iteration", true))
		return -1;

	//
	idx = ptr_array.add(new Point(2000, 3000));
	if (!testPrint((ptr_array[idx]->X == 2000) && (ptr_array[idx]->Y == 3000), "  * Adding an item at the end of the array", true))
		return -1;

	// remove
	delete ptr_array.get(idx);
	ptr_array.remove(idx);
	if (!testPrint(!ptr_array.exists(idx), "  * Removing an item from the array", true))
		return -1;

	// MEMORY
	unsigned int mem_used;

	// mem used
	mem_used = ptr_array.mem_used();
	if (!testPrint(mem_used > 2, "  * Used memory should be nonzero", true))
		return -1;

	// free memory
	for (unsigned int i = ptr_array.first(); i != INVALID_IDX; i = ptr_array.next(i))
		delete ptr_array.get(i);

	ptr_array.remove_all();
	testPrint(true, "  * Freeing memory", true);

	mem_used = ptr_array.mem_used();
	if (!testPrint(mem_used == 0, "  * Used memory should be zero", true))
		return -1;


	return 0;
}

// test Map
int testMap() {
	info("- Testing Map<struct>-----");
/*
	bool r;

	Map<char, Point> pt_map;
	unsigned int idx;

	// fill the array
	info("  * Filling the map with items.");
	pt_map.set("a", Point(100, 200));
	pt_map.set("b", Point(101, 201));
	pt_map.set("c", Point(102, 202));
	pt_map.set("d", Point(103, 203));

	// test the values
	r = true;
	Point pt;
	r &= (pt_map.lookup("a", pt) && pt.X == 100 && pt.Y == 200);
	r &= (pt_map.lookup("b", pt) && pt.X == 101 && pt.Y == 201);
	r &= (pt_map.lookup("c", pt) && pt.X == 102 && pt.Y == 202);
	r &= (pt_map.lookup("d", pt) && pt.X == 103 && pt.Y == 203);
	if (!testPrint(r, "  * Checking if all values were inserted correctly", true))
		return -1;

	//  iterators
	info("  * Testing iterators.");
	r = true;
	int i = 0;
	for (unsigned int iter = pt_map.first(); iter != INVALID_IDX; iter = pt_map.next(iter)) {
		Point pt = pt_map.get(iter);
		r &= (pt.X == i + 100) && (pt.Y == i + 200);
		i++;
	}
	if (!testPrint(r, "    - Forward iteration", true))
		return -1;

	r = true;
	i = 3;
	for (unsigned int iter = pt_map.last(); iter != INVALID_IDX; iter = pt_map.prev(iter)) {
		Point pt = pt_map.get(iter);
		r &= (pt.X == i + 100) && (pt.Y == i + 100);
		i--;
	}
	if (!testPrint(r, "    - Backward iteration", true))
		return -1;


	// non-existent item
	r = (pt_map.lookup("x", pt));
	if (!testPrint(r, "  * Checking non-existent item", false))
		return -1;

	// replace existing item
	pt_map.set("d", Point(104, 204));
	r = (pt_map.lookup("d", pt) && pt.X == 104 && pt.Y == 204);
	if (!testPrint(r, "  * Replacing existing item", true))
		return -1;

	// deleting item
	pt_map.remove("d");
	r = (pt_map.lookup("d", pt));
	if (!testPrint(r, "  * Deleting item", false))
		return -1;

	// count
	unsigned int count = pt_map.count();
	if (!testPrint(count == 3, "  * Number of items == 3", true))
		return -1;

	// free
	pt_map.remove_all();
	if (!testPrint(pt_map.is_empty(), "  * Empty map after removing all items", true))
		return -1;
*/
	return 0;
}

// test MapHSOrd
int testMapHSOrd() {
	info("- Testing MapHSOrd<struct>-----");

	bool r;

	MapHSOrd map;

	// fill the array
	info("  * Filling the map with items.");
	unsigned int k0[] = { 1, 2 };
	map.set(k0, 2, 100);
	unsigned int k1[] = { 1, 3 };
	map.set(k1, 2, 101);
	unsigned int k2[] = { 1, 4 };
	map.set(k2, 2, 102);
	unsigned int k3[] = { 1, 5 };
	map.set(k3, 2, 103);

	// test the values
	r = true;
	unsigned int w;
	r &= (map.lookup(k0, 2, w) && w == 100);
	r &= (map.lookup(k1, 2, w) && w == 101);
	r &= (map.lookup(k2, 2, w) && w == 102);
	r &= (map.lookup(k3, 2, w) && w == 103);
	if (!testPrint(r, "  * Checking if all values were inserted correctly", true))
		return -1;

	// non-existent item
	unsigned int nk[] = { 100, 100 };
	r = (map.lookup(nk, 2, w));
	if (!testPrint(r, "  * Checking non-existent item", false))
		return -1;

	// replace existing item
	map.set(k3, 2, 104);
	r = (map.lookup(k3, 2, w) && w == 104);
	if (!testPrint(r, "  * Replacing existing item", true))
		return -1;

	// deleting item
	map.remove(k3, 2);
	r = (map.lookup(k3, 2, w));
	if (!testPrint(r, "  * Deleting item", false))
		return -1;

	// free
	map.remove_all();

	return 0;
}

// test MapHS
int testMapHS() {
	info("- Testing MapHS<struct>-----");

	bool r;

	MapHS map;

	// fill the array
	info("  * Filling the map with items.");
	map.set((uint8_t *) "a", 1, 100);
	map.set((uint8_t *) "ab", 2, 101);
	map.set((uint8_t *) "abc", 3, 102);
	map.set((uint8_t *) "abcd", 4, 103);

	// test the values
	r = true;
	unsigned int w;
	r &= (map.lookup((uint8_t *) "a", 1, w) && w == 100);
	r &= (map.lookup((uint8_t *) "ab", 2, w) && w == 101);
	r &= (map.lookup((uint8_t *) "abc", 3, w) && w == 102);
	r &= (map.lookup((uint8_t *) "abcd", 4, w) && w == 103);
	if (!testPrint(r, "  * Checking if all values were inserted correctly", true))
		return -1;

	// non-existent item
	r = (map.lookup((uint8_t *) "xyz", 3, w));
	if (!testPrint(r, "  * Checking non-existent item", false))
		return -1;

	// replace existing item
	map.set((uint8_t *) "abcd", 4, 104);
	r = (map.lookup((uint8_t *) "abcd", 4, w) && w == 104);
	if (!testPrint(r, "  * Replacing existing item", true))
		return -1;

	// deleting item
	map.remove((uint8_t *) "abcd", 4);
	r = (map.lookup((uint8_t *) "abcd", 4, w));
	if (!testPrint(r, "  * Deleting item", false))
		return -1;

	// free
	map.remove_all();

	return 0;
}


// test bit array
int testBitArray() {
	BitArray bit_array;
	bool r;
	unsigned int count;
	unsigned int index, value;
	unsigned int mem_used;

	info("- Testing BitArray -----");

	// first set a value in the empty array
	r = bit_array.set(1);
	if (!testPrint(r, "  * Putting a non-existent value (1) into the array", true))
		return -1;

	// put the same value again (it should exist)
	r = bit_array.set(1);
	if (!testPrint(r, "  * Putting the value that already exist in the array", false))
		return -1;

	// test the value
	r = bit_array.is_set(1);
	if (!testPrint(r, "  * Testing if the value exist in the array", true))
		return -1;

	// set next value
	r = bit_array.set(3);
	if (!testPrint(r, "  * Putting another value (3) into the array", true))
		return -1;

	// test the count of items in the array
	count = bit_array.count();
	if (!testPrint(count == 2, "  * Number of values in the array", true))
		return -1;

	// n-th index
	value = bit_array.by_count(2);
	if (!testPrint(value == 3, "  * The value at position 2 should be 3", true))
		return -1;

	// put one more value into the array
	r = bit_array.set(2);
	if (!testPrint(r, "  * Putting one more value (2) into the array", true))
		return -1;

	// ITERATION

	// test the first index
	index = bit_array.first();
	if (!testPrint(index == 1, "  * The first index in the array should be one (1)", true))
		return -1;

	index = bit_array.next(index);
	if (!testPrint(index == 2, "  * The next index in the array should be two (2)", true))
		return -1;

	index = bit_array.last();
	if (!testPrint(index == 3, "  * The last index in the array should be two (3)", true))
		return -1;

	index = bit_array.prev(index);
	if (!testPrint(index == 2, "  * The previous index in the array should be two (2)", true))
		return -1;

	// COPY
	BitArray dup;
	r = dup.copy(&bit_array);
	if (!testPrint(r, "  * Making copy of the array", true))
		return -1;

	// test copied values
	value = bit_array.by_count(1);
	if (!testPrint(value == 1, "    - The value at position 1 should be 1", true))
		return -1;

	value = bit_array.by_count(2);
	if (!testPrint(value == 2, "    - The value at position 2 should be 3", true))
		return -1;

	value = bit_array.by_count(3);
	if (!testPrint(value == 3, "    - The value at position 3 should be 3", true))
		return -1;

	// MEMORY

	// mem used
	mem_used = bit_array.mem_used();
	if (!testPrint(mem_used > 2, "  * Used memory should be nonzero", true))
		return -1;

	// free memory
	bit_array.free();
	testPrint(true, "  * Freeing memory", true);

	mem_used = bit_array.mem_used();
	if (!testPrint(mem_used == 0, "  * Used memory should be zero", true))
		return -1;

	return 0;
}

int main() {
	int ret = 0;

	// test JudyV
	if ((ret = testArray()) != 0)
		return ret;

	if ((ret = testArrayPtrStruct()) != 0)
		return ret;

	// test Map
	if ((ret = testMap()) != 0)
		return ret;

	// test MapHS
	if ((ret = testMapHS()) != 0)
		return ret;

	// test MapHSOrd
	if ((ret = testMapHSOrd()) != 0)
		return ret;

	// test Judy1
	if ((ret = testBitArray()) != 0)
		return ret;

	return ret;
}
