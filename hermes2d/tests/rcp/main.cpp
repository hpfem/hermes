#include "hermes2d.h"

using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcp;
using Teuchos::null;

int main(int argc, char* argv[])
{

    RCP<Mesh> m = rcp(new Mesh());
    Ptr<Mesh> p = m.ptr();
    if (m == null) return ERR_FAILURE;
    if (p == null) return ERR_FAILURE;

    return ERR_SUCCESS;
}
