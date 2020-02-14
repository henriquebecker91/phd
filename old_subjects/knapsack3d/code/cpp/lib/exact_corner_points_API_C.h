#include "altered_third_party/shared_typedefs.h" /* for {i,n,c}type */

struct space;
typedef struct space space;

space* new_space(itype l, itype w, itype h);
void del_space(space *s);

void pop(space *s);

void place(space *s, itype l, itype w, itype h, ctype cp_ix);

boolean can_be_placed(space *s, itype l, itype w, itype h, ctype cp_ix);

const ctype cpnb(space *s);
const itype* cpx(space *s);
const itype* cpy(space *s);
const itype* cpz(space *s);

/* Unnecessary as can_be_placed will check this for us.
const ntype* cp_clog(space *s);
const itype* cpmx(space *s);
const itype* cpmy(space *s);
const itype* cpmz(space *s);
*/

