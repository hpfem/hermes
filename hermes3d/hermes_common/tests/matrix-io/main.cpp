#include <iostream>
#include <stdexcept>

#include "matrix.h"
#include "matrixio.h"
#include "solvers.h"

#define EPS 1e-4
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_matrix_hb()
{
    CooMatrix Acoow(4);
    Acoow.add(0, 0, -1);
    Acoow.add(1, 1, -1);
    Acoow.add(2, 2, -1);
    Acoow.add(3, 3, -1);
    Acoow.add(0, 1, 2);
    Acoow.add(1, 0, 2);
    Acoow.add(1, 2, 2);
    Acoow.add(2, 1, 2);
    Acoow.add(2, 3, 2);
    Acoow.add(3, 2, 2);

    double resw[4] = {1., 1., 1., 1.};

    // coo
    // write
    write_hb_coo("/tmp/coo.rua", &Acoow, resw);
    // read
    CooMatrix *Acoor = read_hb_coo("/tmp/coo.rua");
    Acoow.print();
    Acoor->print();
    remove("/tmp/coo.rua");

    // csc
    CSCMatrix Acscw(&Acoow);
    // write
    write_hb_csc("/tmp/csc.rua", &Acscw, resw);
    // read
    CSCMatrix *Acscr = read_hb_csc("/tmp/csc.rua");
    Acscw.print();
    Acscr->print();
    remove("/tmp/csc.rua");

    // csc
    CSRMatrix Acsrw(&Acoow);
    // write
    write_hb_csr("/tmp/csr.rua", &Acsrw, resw);
    // read
    CSRMatrix *Acsrr = read_hb_csr("/tmp/csr.rua");
    Acsrw.print();
    Acsrr->print();
    remove("/tmp/csr.rua");
}
int main(int argc, char* argv[])
{
    long size;
    char *buf;
    char *ptr;

    size = pathconf(".", _PC_PATH_MAX);

    if ((buf = (char *)malloc((size_t)size)) != NULL)
        ptr = getcwd(buf, (size_t) size);

    std::cout << buf << std::endl;

    try {
        test_matrix_hb();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
